/*
 * This file is part of g-mm-repair
 * <https://github.com/ftosoni/g-mm-repair>.
 * Copyright (c) 2026 Francesco Tosoni.
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 *     http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 */

#include <iostream>
#include <fstream>
#include <vector>
#include <random>
#include <cmath>
#include <algorithm>
#include <string>
#include <filesystem>

#include "grammar.h"
#include "build_schedule.h"

// Unified CUDA-free correctness harness -- what CI runs (see
// tests/run_integration_tests.py). It exercises the whole host pipeline
// (RePair grammar -> schedule -> streaming sweep) and checks the sweep
// bit-for-bit against independent oracles, so the paper's algorithm can be
// validated on any machine without a GPU. See main() for the two entry modes.

// Semiring op-sets (host, double-based -- mirrors gpu_engine.h but CUDA-free so CI
// can validate all three semirings without a GPU). Tropical uses 1e9 as +infinity.
struct SRPlus {
    static double zero() { return 0.0; }
    static double oplus(double a, double b) { return a + b; }
    static double leaf(double V, double x) { return V * x; }
    static double apply(double v, double c) { return c * v; }
    static constexpr const char* name = "plustimes";
};
struct SRBool {
    static double zero() { return 0.0; }
    static double oplus(double a, double b) { return (a != 0.0 || b != 0.0) ? 1.0 : 0.0; }
    static double leaf(double V, double x) { return (V != 0.0 && x != 0.0) ? 1.0 : 0.0; }
    static double apply(double v, double /*c*/) { return v; }
    static constexpr const char* name = "boolean";
};
struct SRTrop {
    static double zero() { return 1e9; }
    static double oplus(double a, double b) { return a < b ? a : b; }
    static double leaf(double V, double x) {
        long iv = (long)V, ix = (long)x;
        if (iv >= 1000000000L || ix >= 1000000000L) return 1e9;
        return (double)(iv + ix);
    }
    static double apply(double v, double /*c*/) { return v; }
    static constexpr const char* name = "tropical";
};

// The same proper-layered, double-buffered sweep the engine runs, on the CPU,
// parameterized by the semiring S.
template <class S>
static std::vector<double> cpu_sweep(const BuiltSchedule& sched,
                                     const std::vector<double>& V,
                                     const std::vector<float>& x,
                                     int rows, int cols) {
    int alpha = sched.alpha;
    std::vector<double> T(alpha, S::zero());
    for (int p = 1; p < alpha; ++p) {
        int t = sched.term_orig[p];
        T[p] = S::leaf(V[(t - 1) / cols], x[(t - 1) % cols]);
    }

    std::vector<double> y(rows, S::zero());
    for (const auto& te : sched.term_emit) {
        y[te.second] = S::oplus(y[te.second], T[te.first]);
    }

    auto cv = [&](int c, const std::vector<double>& prev) -> double {
        return (c >= 0) ? prev[c] : T[-c - 1];
    };

    std::vector<double> prev, cur;
    for (int k = 1; k <= sched.max_depth; ++k) {
        const auto& level = sched.layered_levels[k];
        cur.assign(level.size(), S::zero());
        for (std::size_t idx = 0; idx < level.size(); ++idx) {
            const LayeredRule& r = level[idx];
            cur[idx] = S::oplus(S::apply(cv(r.left, prev), r.coeff), cv(r.right, prev));
        }
        for (const auto& e : sched.emit_by_level[k]) {
            y[e.second] = S::oplus(y[e.second], cur[e.first]);
        }
        prev.swap(cur);
    }
    return y;
}

// Brute-force y = M x straight from the CSRV string, in the semiring S.
template <class S>
static std::vector<double> brute_force(const std::vector<int>& vc,
                                       const std::vector<double>& V,
                                       const std::vector<float>& x,
                                       int rows, int cols) {
    std::vector<double> y(rows, S::zero());
    int row = 0;
    for (int code : vc) {
        if (code == 0) { ++row; continue; }
        y[row] = S::oplus(y[row], S::leaf(V[(code - 1) / cols], x[(code - 1) % cols]));
    }
    return y;
}

// Helper to expand a symbol to its original terminals.
static void expand_symbol_local(const grammar::Grammar& g, int sym, std::vector<int>& out) {
    if (sym < g.alpha) { out.push_back(sym); return; }
    std::vector<std::pair<int,int>> st;
    st.push_back({sym, 1});
    while (!st.empty()) {
        auto [s, reps] = st.back();
        st.pop_back();
        if (s < g.alpha) {
            for (int r = 0; r < reps; ++r) out.push_back(s);
            continue;
        }
        int rule_idx = s - g.alpha;
        if (rule_idx < 0 || rule_idx >= (int)g.rules.size()) {
            std::cerr << "CRITICAL ERROR: Rule index " << rule_idx << " out of bounds of rules size " << g.rules.size() << " for symbol " << s << " (alpha=" << g.alpha << ")\n";
            exit(1);
        }
        const grammar::Rule& rule = g.rules[rule_idx];
        for (int r = 0; r < reps; ++r) {
            if (rule.kind == grammar::RUNLEN) {
                st.push_back({rule.left, rule.rep});
            } else {
                st.push_back({rule.right, 1});
                st.push_back({rule.left, 1});
            }
        }
    }
}

// Reconstruct the full CSRV vector from the compressed grammar.
static std::vector<int> reconstruct_csrv(const grammar::Grammar& g) {
    std::vector<int> out;
    for (int c : g.C) {
        if (c == 0) out.push_back(0);
        else expand_symbol_local(g, c, out);
    }
    return out;
}

// File-test of the RePair construction (the paper's path): load a real grammar
// mm-repair produced for <base_name>, lay it out, sweep on the CPU, and check
// against a brute-force y = M x reconstructed from the same grammar.
int run_repair_filetest(const std::string& base_name, int rows, int cols) {
    std::cout << "Running RePair-construction file-test for matrix base: " << base_name << "...\n";

    // 1. Load rules (.vc.R)
    std::string r_filename = base_name + ".vc.R";
    std::ifstream r_file(r_filename, std::ios::binary);
    if (!r_file.is_open()) {
        std::cerr << "Error: Could not open rules file: " << r_filename << "\n";
        return 1;
    }
    r_file.seekg(0, std::ios::end);
    size_t r_file_size = r_file.tellg();
    r_file.seekg(0, std::ios::beg);

    int alpha = 0;
    r_file.read(reinterpret_cast<char*>(&alpha), sizeof(int));
    size_t nt_num = (r_file_size - sizeof(int)) / (2 * sizeof(int));
    std::vector<int> raw_rules(nt_num * 2);
    r_file.read(reinterpret_cast<char*>(raw_rules.data()), nt_num * 2 * sizeof(int));
    r_file.close();

    // 2. Load top sequence (.vc.C)
    std::string c_filename = base_name + ".vc.C";
    std::ifstream c_file(c_filename, std::ios::binary);
    if (!c_file.is_open()) {
        std::cerr << "Error: Could not open C file: " << c_filename << "\n";
        return 1;
    }
    c_file.seekg(0, std::ios::end);
    size_t c_file_size = c_file.tellg();
    c_file.seekg(0, std::ios::beg);
    std::vector<int> C(c_file_size / sizeof(int));
    c_file.read(reinterpret_cast<char*>(C.data()), c_file_size);
    c_file.close();

    // 3. Load values (.val)
    std::string val_filename = base_name + ".val";
    std::ifstream val_file(val_filename, std::ios::binary);
    if (!val_file.is_open()) {
        std::cerr << "Error: Could not open val file: " << val_filename << "\n";
        return 1;
    }
    val_file.seekg(0, std::ios::end);
    size_t val_file_size = val_file.tellg();
    size_t num_vals = val_file_size / sizeof(double);
    val_file.seekg(0, std::ios::beg);
    std::vector<double> V(num_vals);
    val_file.read(reinterpret_cast<char*>(V.data()), val_file_size);
    val_file.close();

    // Build the grammar object
    grammar::Grammar g;
    g.alpha = alpha;
    g.C = C;
    g.rules.resize(nt_num);
    for (size_t i = 0; i < nt_num; ++i) {
        g.rules[i].kind = grammar::BINARY;
        g.rules[i].left = raw_rules[2 * i];
        g.rules[i].right = raw_rules[2 * i + 1];
        g.rules[i].rep = 1;
    }

    std::cout << "Loaded grammar: alpha=" << alpha << " NTs=" << nt_num << " |C|=" << C.size() << " values=" << num_vals << "\n";

    // Reconstruct CSRV to run brute force comparison
    std::vector<int> vc = reconstruct_csrv(g);
    std::cout << "Reconstructed CSRV length: " << vc.size() << "\n";

    // Build schedule once; the layout is semiring-independent.
    BuiltSchedule sched = build_schedule(g, rows, cols);

    // Cross-check the CPU sweep against the brute-force y = M x for ALL three semirings
    // (each with an appropriate random x): the GPU-free analog of the engine's
    // PlusTimes / Boolean / Tropical runs, so CI validates every semiring on any machine.
    std::mt19937 rng(456);
    bool all_ok = true;
    for (int which = 0; which < 3; ++which) {
        std::vector<float> x(cols);
        for (int j = 0; j < cols; ++j) {
            if (which == 1)      x[j] = (rng() % 2) ? 1.0f : 0.0f;                          // boolean
            else if (which == 2) x[j] = (rng() % 10 == 0) ? 1e9f : (float)(rng() % 100);    // tropical
            else                 x[j] = std::uniform_real_distribution<float>(-1.0f, 1.0f)(rng); // plustimes
        }
        std::vector<double> y_sweep, y_bf;
        const char* nm;
        if (which == 1)      { y_sweep = cpu_sweep<SRBool>(sched, V, x, rows, cols); y_bf = brute_force<SRBool>(vc, V, x, rows, cols); nm = SRBool::name; }
        else if (which == 2) { y_sweep = cpu_sweep<SRTrop>(sched, V, x, rows, cols); y_bf = brute_force<SRTrop>(vc, V, x, rows, cols); nm = SRTrop::name; }
        else                 { y_sweep = cpu_sweep<SRPlus>(sched, V, x, rows, cols); y_bf = brute_force<SRPlus>(vc, V, x, rows, cols); nm = SRPlus::name; }
        double max_abs = 0.0, max_rel = 0.0;
        for (int r = 0; r < rows; ++r) {
            double a = std::abs(y_sweep[r] - y_bf[r]);
            max_abs = std::max(max_abs, a);
            double denom = std::abs(y_bf[r]);
            if (denom > 1e-9) max_rel = std::max(max_rel, a / denom);
        }
        bool ok = max_abs < 1e-6;
        all_ok = all_ok && ok;
        std::cout << "Sweep vs brute force [" << nm << "]: max_abs=" << max_abs
                  << "  max_rel=" << max_rel << (ok ? "  PASS" : "  FAIL") << "\n";
    }
    std::cout << (all_ok ? "REPAIR FILE-TEST PASSED" : "REPAIR FILE-TEST FAILED") << "\n";
    return all_ok ? 0 : 1;
}

// As above, but validated against an external reference vector produced by the
// CPU baseline re32mm (the mm-repair remm binary) -- an independent oracle, not
// this file's own brute force.
int run_repair_filetest_with_reference(const std::string& base_name, int rows, int cols, const std::string& x_file_path, const std::string& y_ref_path) {
    std::cout << "Running RePair-construction file-test against reference vector...\n";

    // 1. Load rules (.vc.R)
    std::string r_filename = base_name + ".vc.R";
    std::ifstream r_file(r_filename, std::ios::binary);
    if (!r_file.is_open()) {
        std::cerr << "Error: Could not open rules file: " << r_filename << "\n";
        return 1;
    }
    r_file.seekg(0, std::ios::end);
    size_t r_file_size = r_file.tellg();
    r_file.seekg(0, std::ios::beg);

    int alpha = 0;
    r_file.read(reinterpret_cast<char*>(&alpha), sizeof(int));
    size_t nt_num = (r_file_size - sizeof(int)) / (2 * sizeof(int));
    std::vector<int> raw_rules(nt_num * 2);
    r_file.read(reinterpret_cast<char*>(raw_rules.data()), nt_num * 2 * sizeof(int));
    r_file.close();

    // 2. Load top sequence (.vc.C)
    std::string c_filename = base_name + ".vc.C";
    std::ifstream c_file(c_filename, std::ios::binary);
    if (!c_file.is_open()) {
        std::cerr << "Error: Could not open C file: " << c_filename << "\n";
        return 1;
    }
    c_file.seekg(0, std::ios::end);
    size_t c_file_size = c_file.tellg();
    c_file.seekg(0, std::ios::beg);
    std::vector<int> C(c_file_size / sizeof(int));
    c_file.read(reinterpret_cast<char*>(C.data()), c_file_size);
    c_file.close();

    // 3. Load values (.val)
    std::string val_filename = base_name + ".val";
    std::ifstream val_file(val_filename, std::ios::binary);
    if (!val_file.is_open()) {
        std::cerr << "Error: Could not open val file: " << val_filename << "\n";
        return 1;
    }
    val_file.seekg(0, std::ios::end);
    size_t val_file_size = val_file.tellg();
    size_t num_vals = val_file_size / sizeof(double);
    val_file.seekg(0, std::ios::beg);
    std::vector<double> V(num_vals);
    val_file.read(reinterpret_cast<char*>(V.data()), val_file_size);
    val_file.close();

    // Build the grammar object
    grammar::Grammar g;
    g.alpha = alpha;
    g.C = C;
    g.rules.resize(nt_num);
    for (size_t i = 0; i < nt_num; ++i) {
        g.rules[i].kind = grammar::BINARY;
        g.rules[i].left = raw_rules[2 * i];
        g.rules[i].right = raw_rules[2 * i + 1];
        g.rules[i].rep = 1;
    }

    std::cout << "Loaded grammar: alpha=" << alpha << " NTs=" << nt_num << " |C|=" << C.size() << " values=" << num_vals << "\n";

    // 4. Load input vector x (as doubles, convert to float)
    std::ifstream x_file(x_file_path, std::ios::binary);
    if (!x_file.is_open()) {
        std::cerr << "Error: Could not open input vector file: " << x_file_path << "\n";
        return 1;
    }
    std::vector<double> x_double(cols);
    x_file.read(reinterpret_cast<char*>(x_double.data()), cols * sizeof(double));
    x_file.close();
    std::vector<float> x(cols);
    for (int j = 0; j < cols; ++j) x[j] = static_cast<float>(x_double[j]);

    // 5. Load expected y vector (as doubles)
    std::ifstream y_file(y_ref_path, std::ios::binary);
    if (!y_file.is_open()) {
        std::cerr << "Error: Could not open expected y vector file: " << y_ref_path << "\n";
        return 1;
    }
    std::vector<double> y_ref(rows);
    y_file.read(reinterpret_cast<char*>(y_ref.data()), rows * sizeof(double));
    y_file.close();

    // Build schedule
    BuiltSchedule sched = build_schedule(g, rows, cols);

    // Compute sweep (re32mm is a PlusTimes CPU mat-vec)
    std::vector<double> y_sweep = cpu_sweep<SRPlus>(sched, V, x, rows, cols);

    double max_abs = 0.0, max_rel = 0.0;
    for (int r = 0; r < rows; ++r) {
        double a = std::abs(y_sweep[r] - y_ref[r]);
        max_abs = std::max(max_abs, a);
        double denom = std::abs(y_ref[r]);
        if (denom > 1e-9) max_rel = std::max(max_rel, a / denom);
    }
    std::cout << "Sweep vs re32mm CPU reference: max_abs=" << max_abs << "  max_rel=" << max_rel << "\n";

    bool ok = max_abs < 1e-4; // float precision difference might cause small drift
    std::cout << (ok ? "REPAIR REFERENCE TEST PASSED" : "REPAIR REFERENCE TEST FAILED") << "\n";
    return ok ? 0 : 1;
}

int run_unittest() {
    std::cout << "Running unit test covering all 4 kernel cases...\n";

    // Build the grammar object
    grammar::Grammar g;
    g.alpha = 5; // terminals 0, 1, 2, 3, 4
    
    // Rule 0 (BINARY): N_0 -> 1, 2 (T, T)
    // Rule 1 (BINARY): N_1 -> 3, N_0 (T, NT)
    // Rule 2 (RUNLEN): N_2 -> N_1^3 (RUNLEN)
    // Rule 3 (BINARY): N_3 -> N_2, N_0 (NT, NT)
    g.rules.resize(4);
    
    g.rules[0].kind = grammar::BINARY;
    g.rules[0].left = 1;
    g.rules[0].right = 2;
    g.rules[0].rep = 1;
    
    g.rules[1].kind = grammar::BINARY;
    g.rules[1].left = 3;
    g.rules[1].right = 5; // N_0 is (alpha + 0) = 5
    g.rules[1].rep = 1;
    
    g.rules[2].kind = grammar::RUNLEN;
    g.rules[2].left = 6; // N_1 is (alpha + 1) = 6
    g.rules[2].right = -1;
    g.rules[2].rep = 3;
    
    g.rules[3].kind = grammar::BINARY;
    g.rules[3].left = 7; // N_2 is (alpha + 2) = 7
    g.rules[3].right = 5; // N_0 is (alpha + 0) = 5
    g.rules[3].rep = 1;

    // Top sequence C has N_3 (alpha + 3 = 8) and row delimiter 0.
    g.C = {8, 0};

    int rows = 2;
    int cols = 3;
    
    std::vector<double> V = {1.5, 2.5};
    std::vector<float> x = {0.5f, -0.5f, 1.0f};

    BuiltSchedule sched = build_schedule(g, rows, cols);

    std::cout << "Unittest schedule details:\n";
    std::cout << "  max_depth: " << sched.max_depth << "\n";
    std::cout << "  max_width: " << sched.max_width << "\n";
    std::cout << "  num_layered: " << sched.num_layered << "\n";
    std::cout << "  num_passthrough: " << sched.num_passthrough << "\n";

    // Verify expectations
    if (sched.max_depth != 4) {
        std::cerr << "Assertion failed: max_depth should be 4, got " << sched.max_depth << "\n";
        return 1;
    }
    if (sched.num_passthrough != 5) {
        std::cerr << "Assertion failed: num_passthrough should be 5, got " << sched.num_passthrough << "\n";
        return 1;
    }

    // Reconstruct CSRV to run brute force comparison
    std::vector<int> vc = reconstruct_csrv(g);
    std::cout << "Reconstructed CSRV length: " << vc.size() << "\n";

    // Compute sweep and brute force
    std::vector<double> y_sweep = cpu_sweep<SRPlus>(sched, V, x, rows, cols);
    std::vector<double> y_bf    = brute_force<SRPlus>(vc, V, x, rows, cols);

    std::cout << "y_sweep: [" << y_sweep[0] << ", " << y_sweep[1] << "]\n";
    std::cout << "y_bf:    [" << y_bf[0] << ", " << y_bf[1] << "]\n";

    double max_abs = 0.0;
    for (int r = 0; r < rows; ++r) {
        max_abs = std::max(max_abs, std::abs(y_sweep[r] - y_bf[r]));
    }

    if (max_abs > 1e-6) {
        std::cerr << "Assertion failed: Sweep vs brute force max_abs = " << max_abs << " > 1e-6\n";
        return 1;
    }

    // Check expected values: y[0] should be 4.5, y[1] should be 0.0
    if (std::abs(y_sweep[0] - 4.5) > 1e-6) {
        std::cerr << "Assertion failed: y_sweep[0] should be 4.5, got " << y_sweep[0] << "\n";
        return 1;
    }
    if (std::abs(y_sweep[1] - 0.0) > 1e-6) {
        std::cerr << "Assertion failed: y_sweep[1] should be 0.0, got " << y_sweep[1] << "\n";
        return 1;
    }

    std::cout << "UNIT TEST PASSED SUCCESSFULLY!\n";
    return 0;
}

// Dispatch by argument count:
//   --unittest -> Run programmatic unit test covering all 4 kernel cases
//   base rows cols -> RePair file-test vs this file's brute force
//   + x_file y_file -> RePair file-test vs an external re32mm reference vector
int main(int argc, char** argv) {
    if (argc == 2 && std::string(argv[1]) == "--unittest") {
        return run_unittest();
    }
    if (argc < 4) {
        std::cerr << "Usage: " << argv[0] << " <matrix_base_name> <rows> <cols> [x_file y_file]" << std::endl;
        std::cerr << "       " << argv[0] << " --unittest" << std::endl;
        return 1;
    } else if (argc == 4) {
        // Run file-based RePair-construction test
        std::string base_name = argv[1];
        int rows = std::stoi(argv[2]);
        int cols = std::stoi(argv[3]);
        return run_repair_filetest(base_name, rows, cols);
    } else {
        // Run file-based RePair-construction test against re32mm reference
        std::string base_name = argv[1];
        int rows = std::stoi(argv[2]);
        int cols = std::stoi(argv[3]);
        std::string x_file = argv[4];
        std::string y_file = argv[5];
        return run_repair_filetest_with_reference(base_name, rows, cols, x_file, y_file);
    }
}
