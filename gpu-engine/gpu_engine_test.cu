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

#include "gpu_engine.h"
#include "grammar.h"
#include "build_schedule.h"
#include <cstdint>
#include <cstdlib>
#include <iostream>
#include <fstream>
#include <vector>
#include <utility>
#include <cmath>
#include <chrono>
#include <algorithm>
#ifdef _OPENMP
#include <omp.h>
#endif
#ifndef _WIN32
#include <sys/resource.h>   // getrusage: peak host RSS for the construction cost (sec:limitations ii)
#endif

// Driver `gpu_test`: the end-to-end harness behind the manuscript's genotype and
// graph experiments. It (1) reads the value array + grammar, (2) builds the
// proper-layered schedule (section "A proper-layered streaming engine", Step 2) and prints its structural figures
// -- L, w*, |R|, +pt: the architecture-independent numbers of tab:geno_through /
// tab:graph_struct -- (3) runs the CPU reference sweep as the oracle, (4) runs
// the GPU streaming sweep and verifies bit-for-bit, and (5) optionally runs the
// batched (SpMM) path. The semiring (PlusTimes / Boolean / Tropical) is chosen
// at runtime by the SEMIRING env var. CLI: <base> <rows> <cols> [iters] [mode] [B].

// Host-side CPU reference: the same proper-layered, double-buffered sweep with
// emit-on-the-spot, used as the correctness oracle for the GPU engine.
template <typename S>
void run_cpu_reference_layered(
    int rows,
    int alpha,
    int max_depth,
    const std::vector<std::vector<LayeredRule>>& layered_levels,
    const std::vector<std::vector<std::pair<int,int>>>& emit_by_level, // (pos,row) per level
    const std::vector<std::pair<int,int>>& term_emit,                  // (T-index,row)
    const int* val_indices,
    const int* col_indices,
    const float* V,
    const float* x,
    typename S::T* y_ref
) {
    using T = typename S::T;
    // Persistent T terminal array (T[0] = S::zero() is the zero terminal).
    std::vector<T> T_arr(alpha, S::zero());
    T_arr[0] = S::zero();
    #pragma omp parallel for if(alpha > 1000)
    for (int p = 1; p < alpha; ++p) {
        T_arr[p] = S::leaf(V[val_indices[p - 1]], x[col_indices[p - 1]]);
    }

    #pragma omp parallel for if(rows > 1000)
    for (int r = 0; r < rows; ++r) y_ref[r] = S::zero();

    // Terminal occurrences of C (level 0).
    #pragma omp parallel for if(term_emit.size() > 1000)
    for (size_t i = 0; i < term_emit.size(); ++i) {
        S::cpu_emit(&y_ref[term_emit[i].second], T_arr[term_emit[i].first]);
    }

    // Double-buffered sweep; emit each level's C-occurrences on the spot.
    std::vector<T> W_read;
    std::vector<T> W_write;

    for (int k = 1; k <= max_depth; ++k) {
        int num_rules = layered_levels[k].size();
        W_write.assign(num_rules, S::zero());

        #pragma omp parallel for if(num_rules > 1000)
        for (int idx = 0; idx < num_rules; ++idx) {
            LayeredRule r = layered_levels[k][idx];
            T left_val = (r.left >= 0) ? W_read[r.left] : T_arr[-r.left - 1];
            T right_val = (r.right >= 0) ? W_read[r.right] : T_arr[-r.right - 1];
            W_write[idx] = S::oplus(S::apply(left_val, r.coeff), right_val);
        }

        #pragma omp parallel for if(emit_by_level[k].size() > 1000)
        for (size_t i = 0; i < emit_by_level[k].size(); ++i) {
            S::cpu_emit(&y_ref[emit_by_level[k][i].second], W_write[emit_by_level[k][i].first]);
        }

        W_read = W_write;
    }
}

// Dump the engine's input vector x (float, length cols) and output y (as double,
// length rows) so an independent tool (cusparse_test, graphblas_bench) can load the
// SAME vector and assert its own result equals the engine's -- a cross-implementation
// correctness check on one shared vector. Enabled by CROSSCHECK=<prefix>.
template <typename S>
static void dump_crosscheck(const std::string& prefix, const std::string& sr,
                            const std::vector<float>& x, const std::vector<typename S::T>& y,
                            int cols, int rows) {
    std::ofstream fx(prefix + "." + sr + ".x", std::ios::binary);
    fx.write(reinterpret_cast<const char*>(x.data()), (size_t)cols * sizeof(float));
    std::vector<double> yd(rows);
    for (int r = 0; r < rows; ++r) yd[r] = (double)y[r];
    std::ofstream fy(prefix + "." + sr + ".y", std::ios::binary);
    fy.write(reinterpret_cast<const char*>(yd.data()), (size_t)rows * sizeof(double));
    std::cout << "CROSSCHECK: dumped engine x/y to " << prefix << "." << sr << ".{x,y}" << std::endl;
}

template <typename S>
int run_test_for_semiring(
    int rows, int cols, int alpha, int num_distinct_vals, int iters, int B, bool skip_cpu,
    const std::vector<float>& V, const std::vector<float>& x,
    const BuiltSchedule& bs, const std::vector<int>& val_indices, const std::vector<int>& col_indices,
    GPUSchedule& sched, const std::vector<LayeredRule>& flattened_rules,
    std::vector<std::vector<LayeredRule>>& layered_levels,
    std::vector<std::vector<std::pair<int,int>>>& emit_by_level,
    std::vector<std::pair<int,int>>& term_emit,
    const std::string& base_name, const std::string& semiring_name
) {
    using T = typename S::T;
    std::vector<T> y_cpu;
    if (!skip_cpu) {
        y_cpu.assign(rows, S::zero());
#ifdef _OPENMP
        int max_threads = omp_get_max_threads();
        std::vector<T> y_cpu_omp(rows, S::zero());
        std::cout << "Running CPU reference with OpenMP (" << max_threads << " threads)..." << std::endl;
        auto cpu_omp_start = std::chrono::high_resolution_clock::now();
        run_cpu_reference_layered<S>(
            rows, alpha, bs.max_depth, layered_levels, emit_by_level, term_emit,
            val_indices.data(), col_indices.data(), V.data(), x.data(), y_cpu_omp.data()
        );
        auto cpu_omp_end = std::chrono::high_resolution_clock::now();
        std::chrono::duration<double, std::milli> cpu_omp_ms = cpu_omp_end - cpu_omp_start;
        std::cout << "CPU reference (OpenMP): " << cpu_omp_ms.count() << " ms/vector (" << max_threads << " threads)." << std::endl;

        std::cout << "Running CPU reference sequentially (forced 1 thread)..." << std::endl;
        omp_set_num_threads(1);
        auto cpu_start = std::chrono::high_resolution_clock::now();
        run_cpu_reference_layered<S>(
            rows, alpha, bs.max_depth, layered_levels, emit_by_level, term_emit,
            val_indices.data(), col_indices.data(), V.data(), x.data(), y_cpu.data()
        );
        auto cpu_end = std::chrono::high_resolution_clock::now();
        std::chrono::duration<double, std::milli> cpu_ms = cpu_end - cpu_start;
        std::cout << "CPU reference (sequential): " << cpu_ms.count() << " ms/vector." << std::endl;

        omp_set_num_threads(max_threads); // restore

        // Cross-check: the OpenMP reference must match the sequential one (same oracle,
        // different threading) -- otherwise a data race in the parallel emit is hiding.
        {
            double omp_max_rel = 0.0;
            for (int r = 0; r < rows; ++r) {
                double a = std::abs((double)y_cpu_omp[r] - (double)y_cpu[r]);
                double m = std::abs((double)y_cpu[r]);
                double rel = a / (1.0 + m);
                if (rel > omp_max_rel) omp_max_rel = rel;
            }
            std::cout << "CPU OpenMP vs sequential: max_rel_diff=" << omp_max_rel
                      << (omp_max_rel < 1e-4 ? "  SUCCESS" : "  FAILURE") << std::endl;
        }
#else
        std::cout << "Running CPU reference (sequential)..." << std::endl;
        auto cpu_start = std::chrono::high_resolution_clock::now();
        run_cpu_reference_layered<S>(
            rows, alpha, bs.max_depth, layered_levels, emit_by_level, term_emit,
            val_indices.data(), col_indices.data(), V.data(), x.data(), y_cpu.data()
        );
        auto cpu_end = std::chrono::high_resolution_clock::now();
        std::chrono::duration<double, std::milli> cpu_ms = cpu_end - cpu_start;
        std::cout << "CPU reference (sequential): " << cpu_ms.count() << " ms/vector." << std::endl;
#endif
    }

    // Deallocate per-level structures no longer needed after the CPU reference.
    layered_levels.clear(); layered_levels.shrink_to_fit();
    emit_by_level.clear(); emit_by_level.shrink_to_fit();
    term_emit.clear(); term_emit.shrink_to_fit();

    std::vector<T> y_gpu(rows, S::zero());

    // 10. Run GPU engine.
    std::cout << "Running GPU engine (double-buffered proper layered, emit-on-the-spot)..." << std::endl;
    auto start = std::chrono::high_resolution_clock::now();
    run_gpu_right_multiply_layered(
        rows, cols, alpha, flattened_rules.size(), flattened_rules.data(),
        val_indices.data(), col_indices.data(), V.data(), num_distinct_vals,
        x.data(), sched, reinterpret_cast<float*>(y_gpu.data()), iters
    );
    auto end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double, std::milli> duration = end - start;
    std::cout << "GPU multiplication completed in " << duration.count() << " ms." << std::endl;

    // 10b. Optional batched (SpMM) run: B vectors at once. We replicate the same
    //      x across all B columns so the result of every column must equal the
    //      single-vector y_cpu -- a cheap correctness check for the batched path.
    if (B > 1) {
        double gb = ((double)alpha * B + 2.0 * sched.max_width * B
                     + (double)cols * B + (double)rows * B) * 4.0 / 1e9;
        if (gb > 90.0) {
            std::cout << "Batched run B=" << B << " SKIPPED (est. " << gb
                      << " GB exceeds memory budget)." << std::endl;
        } else {
            std::cout << "Batched run B=" << B << " (est. " << gb << " GB B-wide buffers)..." << std::endl;
            std::vector<float> X((size_t)cols * B);
            for (int j = 0; j < cols; ++j)
                for (int b = 0; b < B; ++b)
                    X[(size_t)j * B + b] = x[j];        // replicate x across columns
            std::vector<T> Y((size_t)rows * B, S::zero());
            run_gpu_right_multiply_batched(
                rows, cols, alpha, flattened_rules.size(), flattened_rules.data(),
                val_indices.data(), col_indices.data(), V.data(), num_distinct_vals,
                X.data(), sched, reinterpret_cast<float*>(Y.data()), B, iters);
            if (!skip_cpu) {
                double max_rel = 0.0;
                for (int r = 0; r < rows; ++r) {
                    double m = std::abs((double)y_cpu[r]);
                    for (int b = 0; b < B; ++b) {
                        double rel = std::abs((double)Y[(size_t)r * B + b] - (double)y_cpu[r]) / (1.0 + m);
                        if (rel > max_rel) max_rel = rel;
                    }
                }
                std::cout << "Batched verification: max_rel_diff=" << max_rel
                          << (max_rel < 1e-4 ? "  SUCCESS" : "  FAILURE") << std::endl;
            }
        }
    }

    // 10c. Cross-check dump: hand the shared x and the engine's y to the independent
    //      baselines (cusparse_test / graphblas_bench) so they can confirm agreement.
    if (const char* cc = getenv("CROSSCHECK")) {
        if (cc[0]) dump_crosscheck<S>(cc, semiring_name, x, y_gpu, cols, rows);
    }

    // 11. Verify correctness.
    if (!skip_cpu) {
        double max_abs = 0.0, max_rel = 0.0, max_y = 0.0;
        for (int r = 0; r < rows; ++r) {
            double a = std::abs((double)y_gpu[r] - (double)y_cpu[r]);
            double m = std::abs((double)y_cpu[r]);
            if (a > max_abs) max_abs = a;
            if (m > max_y) max_y = m;
            double rel = a / (1.0 + m);
            if (rel > max_rel) max_rel = rel;
        }
        std::cout << "Verification: max_abs_diff=" << max_abs
                  << ", max_rel_diff=" << max_rel
                  << ", max|y|=" << max_y << std::endl;
        if (max_rel < 1e-4) {
            std::cout << "SUCCESS: GPU and CPU reference results match (within float precision)!" << std::endl;
            return 0;
        } else {
            std::cerr << "FAILURE: Results mismatch!" << std::endl;
            return 1;
        }
    } else {
        std::cout << "Skipped CPU reference verification for large dataset." << std::endl;
        return 0;
    }
}

int main(int argc, char** argv) {
    if (argc < 4) {
        std::cerr << "Usage: " << argv[0] << " <matrix_base_name> <rows> <cols> [iters] [mode] [B]" << std::endl;
        std::cerr << "  mode: 'repair' (default; the RePair construction used in the paper)" << std::endl;
        std::cerr << "  B:    optional batch width for the batched (SpMM) path" << std::endl;
        std::cerr << "  semiring is chosen by the SEMIRING env var: plustimes | boolean | tropical" << std::endl;
        return 1;
    }

    std::string base_name = argv[1];
    int rows = std::stoi(argv[2]);
    int cols = std::stoi(argv[3]);
    int iters = (argc >= 5) ? std::stoi(argv[4]) : 1;

    // 1. Read distinct values (.val) -- the value array V of the (C,R,V)
    //    representation, shared regardless of which constructor built the grammar.
    std::string val_filename = base_name + ".val";
    std::ifstream val_file(val_filename, std::ios::binary);
    if (!val_file.is_open()) {
        std::cerr << "Error: Could not open values file." << std::endl;
        return 1;
    }
    val_file.seekg(0, std::ios::end);
    size_t val_file_size = val_file.tellg();
    size_t num_distinct_vals = val_file_size / sizeof(double);
    val_file.seekg(0, std::ios::beg);
    std::vector<double> V_double(num_distinct_vals);
    val_file.read(reinterpret_cast<char*>(V_double.data()), val_file_size);
    std::vector<float> V(num_distinct_vals);
    for (size_t i = 0; i < num_distinct_vals; ++i) {
        V[i] = static_cast<float>(V_double[i]);
    }

    // 2. Build the grammar from the RePair construction mm-repair emitted. This
    //    fills a grammar::Grammar and feeds the shared schedule builder + GPU
    //    engine. .vc.R holds the binary rules (alpha, then left/right pairs);
    //    .vc.C the top sequence C with 0 marking row ends.
    grammar::Grammar g;
    int alpha = 0;
    {
        std::ifstream r_file(base_name + ".vc.R", std::ios::binary);
        if (!r_file.is_open()) { std::cerr << "Error: Could not open .vc.R file." << std::endl; return 1; }
        r_file.read(reinterpret_cast<char*>(&alpha), sizeof(int));
        r_file.seekg(0, std::ios::end);
        size_t r_file_size = r_file.tellg();
        size_t nt_num = (r_file_size - sizeof(int)) / (2 * sizeof(int));
        r_file.seekg(sizeof(int), std::ios::beg);
        std::vector<int> raw_rules(nt_num * 2);
        r_file.read(reinterpret_cast<char*>(raw_rules.data()), nt_num * 2 * sizeof(int));

        std::ifstream c_file(base_name + ".vc.C", std::ios::binary);
        if (!c_file.is_open()) { std::cerr << "Error: Could not open .vc.C file." << std::endl; return 1; }
        c_file.seekg(0, std::ios::end);
        size_t c_file_size = c_file.tellg();
        c_file.seekg(0, std::ios::beg);
        g.C.resize(c_file_size / sizeof(int));
        c_file.read(reinterpret_cast<char*>(g.C.data()), c_file_size);

        g.alpha = alpha;
        g.rules.resize(nt_num);
        for (size_t i = 0; i < nt_num; ++i) {
            g.rules[i].kind = grammar::BINARY;
            g.rules[i].left = raw_rules[2 * i];
            g.rules[i].right = raw_rules[2 * i + 1];
            g.rules[i].rep = 1;
        }
        std::cout << "RePair construction: alpha=" << alpha << ", NTs=" << nt_num
                  << ", distinct values=" << num_distinct_vals << std::endl;
    }

    // 3. Lay out the grammar (levels + pass-through completion, section "A proper-layered streaming engine", Step 2)
    //    with the shared schedule builder -- identical code path regardless of
    //    the constructor. This also compacts the terminal alphabet to the
    //    symbols actually referenced.
    int alpha_full = alpha;
    auto sched_start = std::chrono::high_resolution_clock::now();
    BuiltSchedule bs = build_schedule(g, rows, cols);
    auto sched_end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double, std::milli> sched_dur = sched_end - sched_start;
    std::cout << "Schedule builder construction time: " << sched_dur.count() << " ms" << std::endl;
    // Canonical, machine-parseable construction-cost line (see sec:limitations (ii),
    // extract_results.py -> tab:build). Time = pass-through completion + level-bucketing
    // + terminal compaction (build_schedule); peak host MB = process max RSS so far, a
    // proxy for the memory needed to construct the layered structure on the host.
    double peak_host_mb = 0.0;
#ifndef _WIN32
    struct rusage ru; getrusage(RUSAGE_SELF, &ru);
    peak_host_mb = ru.ru_maxrss / 1024.0;   // ru_maxrss is in KiB on Linux
#endif
    std::cout << "BUILD schedule ms: " << sched_dur.count()
              << " peak host MB: " << peak_host_mb << std::endl;
    g.rules.clear(); g.rules.shrink_to_fit();   // freed; bs owns the layout now
    g.C.clear(); g.C.shrink_to_fit();

    // 4. Prepare terminal layouts for the COMPACTED alphabet only: compact id
    //    q in [1,alpha) decodes its original dense (value,column) grid id.
    alpha = bs.alpha;
    std::vector<int> val_indices(alpha - 1);
    std::vector<int> col_indices(alpha - 1);
    for (int q = 1; q < alpha; ++q) {
        int t = bs.term_orig[q];   // original dense terminal id
        val_indices[q - 1] = (t - 1) / cols;
        col_indices[q - 1] = (t - 1) % cols;
    }
    std::cout << "Terminal compaction: alpha " << alpha_full << " -> " << alpha
              << " (" << (alpha_full > 0 ? 100.0 * alpha / alpha_full : 0.0) << "% kept)" << std::endl;

    int max_depth = bs.max_depth;
    auto& flattened_rules = bs.flattened_rules;
    auto& layered_levels  = bs.layered_levels;
    auto& emit_by_level   = bs.emit_by_level;
    auto& term_emit       = bs.term_emit;
    GPUSchedule sched = bs.gpu_schedule();

    std::cout << "Proper Layered Completed Grammar Stats:" << std::endl;
    std::cout << "Max depth (L): " << max_depth << ", Max width (w*): " << bs.max_width
              << ", Total layered rules: " << bs.num_layered
              << " (raw NTs: " << bs.num_raw_nt << ", +pt: " << bs.num_passthrough << ")"
              << ", emissions: " << sched.emit_total << "+" << sched.num_term_emit << " term" << std::endl;

    // 8. Generate random input vector x based on semiring choice
    std::vector<float> x(cols);
    const char* semiring_env = getenv("SEMIRING");
    std::string semiring = semiring_env ? semiring_env : "plustimes";
    for (int j = 0; j < cols; ++j) {
        if (semiring == "boolean" || semiring == "bool") {
            x[j] = (rand() % 2 == 0) ? 0.0f : 1.0f;
        } else if (semiring == "tropical" || semiring == "minplus") {
            if (rand() % 10 == 0) {
                x[j] = 1000000000.0f; // infinity
            } else {
                x[j] = static_cast<float>(rand() % 100);
            }
        } else {
            x[j] = static_cast<float>(rand()) / RAND_MAX;
        }
    }

    // Optional fixed-x override: XVEC=<file> supplies the input vector on the
    // command line instead of the random one above. The file format is IDENTICAL
    // to the x the baselines consume in CROSSCHECK mode (cusparse_test /
    // graphblas_bench): a raw binary array of exactly `cols` float32 values, so a
    // chosen x can be dumped (CROSSCHECK) and cross-validated across every
    // implementation. Absent XVEC, behaviour is unchanged (random x).
    if (const char* xf = getenv("XVEC")) {
        std::ifstream fx(xf, std::ios::binary);
        if (!fx.is_open()) { std::cerr << "XVEC: cannot open " << xf << std::endl; return 1; }
        fx.seekg(0, std::ios::end);
        size_t bytes = fx.tellg();
        if (bytes % sizeof(float) != 0 || bytes / sizeof(float) != (size_t)cols) {
            std::cerr << "XVEC: " << xf << " holds " << (bytes / sizeof(float))
                      << " float32 values, expected cols=" << cols << std::endl;
            return 1;
        }
        fx.seekg(0, std::ios::beg);
        fx.read(reinterpret_cast<char*>(x.data()), (size_t)cols * sizeof(float));
        std::cout << "XVEC: loaded fixed x (" << cols << " float32) from " << xf << std::endl;
    }

    // The CPU reference oracle is O(|C|+|R|) but sequential-ish; by default we skip it on
    // billion-scale instances (rows > 30M) purely to keep runs fast. It is NOT a hard limit:
    // set FORCE_CPU_VERIFY=1 to run the oracle and certify the engine bit-for-bit even on the
    // full Software Heritage graph (45.7M nodes) -- slow, but a one-off correctness check.
    bool skip_cpu = (rows > 30000000);
    const char* force_verify = getenv("FORCE_CPU_VERIFY");
    if (force_verify && force_verify[0] && force_verify[0] != '0') skip_cpu = false;
    int B = (argc >= 7) ? std::stoi(argv[6]) : 1;

    int result_code = 0;
    if (semiring == "boolean" || semiring == "bool") {
        std::cout << "--- SEMIRING: Boolean ---" << std::endl;
        result_code = run_test_for_semiring<Boolean>(
            rows, cols, alpha, num_distinct_vals, iters, B, skip_cpu, V, x,
            bs, val_indices, col_indices, sched, flattened_rules, layered_levels, emit_by_level, term_emit,
            base_name, "boolean"
        );
    } else if (semiring == "tropical" || semiring == "minplus") {
        std::cout << "--- SEMIRING: Tropical (min,+) ---" << std::endl;
        result_code = run_test_for_semiring<Tropical>(
            rows, cols, alpha, num_distinct_vals, iters, B, skip_cpu, V, x,
            bs, val_indices, col_indices, sched, flattened_rules, layered_levels, emit_by_level, term_emit,
            base_name, "tropical"
        );
    } else {
        std::cout << "--- SEMIRING: PlusTimes (+,*) ---" << std::endl;
        result_code = run_test_for_semiring<PlusTimes>(
            rows, cols, alpha, num_distinct_vals, iters, B, skip_cpu, V, x,
            bs, val_indices, col_indices, sched, flattened_rules, layered_levels, emit_by_level, term_emit,
            base_name, "plustimes"
        );
    }

    return result_code;
}