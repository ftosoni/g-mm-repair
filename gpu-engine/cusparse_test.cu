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

#include "grammar.h"
#include <iostream>
#include <fstream>
#include <vector>
#include <chrono>
#include <cmath>
#include <cstdint>
#include <cuda_runtime.h>
#include <cusparse.h>
#include <nvml.h>

#define CHECK_CUDA(func) \
    { \
        cudaError_t status = (func); \
        if (status != cudaSuccess) { \
            std::cerr << "CUDA Error at line " << __LINE__ << ": " << cudaGetErrorString(status) << std::endl; \
            return 1; \
        } \
    }

#define CHECK_CUSPARSE(func) \
    { \
        cusparseStatus_t status = (func); \
        if (status != CUSPARSE_STATUS_SUCCESS) { \
            std::cerr << "cuSPARSE Error at line " << __LINE__ << ": " << status << std::endl; \
            return 1; \
        } \
    }

// Recursively expand a grammar symbol to reconstruct original terminals
void expand_symbol(int s, int alpha, const std::vector<grammar::Rule>& rules,
                   int cols, const std::vector<float>& V,
                   std::vector<int>& col_indices, std::vector<float>& values) {
    if (s < alpha) {
        if (s > 0) {
            int col = (s - 1) % cols;
            int val_idx = (s - 1) / cols;
            col_indices.push_back(col);
            values.push_back(V[val_idx]);
        }
    } else {
        int rule_idx = s - alpha;
        const auto& r = rules[rule_idx];
        if (r.kind == grammar::BINARY) {
            expand_symbol(r.left, alpha, rules, cols, V, col_indices, values);
            expand_symbol(r.right, alpha, rules, cols, V, col_indices, values);
        } else if (r.kind == grammar::RUNLEN) {
            for (int t = 0; t < r.rep; ++t) {
                expand_symbol(r.left, alpha, rules, cols, V, col_indices, values);
            }
        }
    }
}

int main(int argc, char** argv) {
    if (argc < 4) {
        std::cerr << "Usage: " << argv[0] << " <matrix_base_name> <rows> <cols> [iters] [mode]" << std::endl;
        std::cerr << "  mode: 'repair' (default)" << std::endl;
        return 1;
    }

    std::string base_name = argv[1];
    int rows = std::stoi(argv[2]);
    int cols = std::stoi(argv[3]);
    int iters = (argc >= 5) ? std::stoi(argv[4]) : 100;
    std::string mode = (argc >= 6) ? std::string(argv[5]) : std::string("repair");
    if (mode == "recompress" || mode == "partII" || mode == "II") {
        std::cerr << "Error: Recompression mode (Part II) is no longer supported/needed." << std::endl;
        return 1;
    }

    // 1. Read distinct values (.val)
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

    // 2. Build the grammar representation
    grammar::Grammar g;
    int alpha = 0;
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

    // 3. Decompress grammar to reconstruct original CSR matrix
    std::vector<int> h_row_ptr;
    std::vector<int> h_col_indices;
    std::vector<float> h_values;

    h_row_ptr.reserve(rows + 1);
    h_row_ptr.push_back(0);

    for (int s : g.C) {
        if (s == 0) {
            h_row_ptr.push_back(h_col_indices.size());
        } else {
            expand_symbol(s, g.alpha, g.rules, cols, V, h_col_indices, h_values);
        }
    }

    // Safety padding for row pointers if the string was cut off
    while (h_row_ptr.size() <= (size_t)rows) {
        h_row_ptr.push_back(h_col_indices.size());
    }

    int nnz = h_values.size();
    std::cout << "Reconstructed CSR Matrix Stats:" << std::endl;
    std::cout << "Rows: " << rows << ", Cols: " << cols << ", NNZ: " << nnz << std::endl;

    // 4. Input vector x. In CROSSCHECK mode we load the SAME vector the engine used
    //    (dumped by gpu_test as "<prefix>.<sr>.x") so this independent CSR
    //    reconstruction can confirm the engine's output bit-for-bit / within float
    //    precision. Otherwise we generate a random vector for the timing baseline.
    const char* cc = getenv("CROSSCHECK");
    const char* sr_env = getenv("SEMIRING");
    std::string sr = sr_env ? sr_env : "plustimes";
    if (sr == "bool") sr = "boolean";
    if (sr == "minplus") sr = "tropical";
    std::vector<float> h_x(cols);
    if (cc && cc[0]) {
        std::string xf = std::string(cc) + "." + sr + ".x";
        std::ifstream fx(xf, std::ios::binary);
        if (!fx.is_open()) { std::cerr << "CROSSCHECK: cannot open " << xf << std::endl; return 1; }
        fx.read(reinterpret_cast<char*>(h_x.data()), (size_t)cols * sizeof(float));
        std::cout << "CROSSCHECK: loaded engine x from " << xf << std::endl;
    } else {
        for (int j = 0; j < cols; ++j) h_x[j] = static_cast<float>(rand()) / RAND_MAX;
    }

    // Loader for the engine's dumped y (double, length rows).
    auto load_engine_y = [&](std::vector<double>& ey) -> bool {
        std::string yf = std::string(cc) + "." + sr + ".y";
        std::ifstream fy(yf, std::ios::binary);
        if (!fy.is_open()) { std::cerr << "CROSSCHECK: cannot open " << yf << std::endl; return false; }
        ey.resize(rows);
        fy.read(reinterpret_cast<char*>(ey.data()), (size_t)rows * sizeof(double));
        return true;
    };

    // Tropical crosscheck: cuSPARSE has no min-plus kernel, so we validate the engine
    // against an independent CPU min-plus SpMV over the SAME reconstructed CSR (same M,
    // same 1e9 infinity convention as gpu_engine.h) and stop before the +x machinery.
    if (cc && cc[0] && sr == "tropical") {
        const double INF = 1e9;
        std::vector<double> y_trop(rows, INF);
        for (int r = 0; r < rows; ++r) {
            double best = INF;
            for (int idx = h_row_ptr[r]; idx < h_row_ptr[r + 1]; ++idx) {
                double mv = (double)h_values[idx];
                double xv = (double)h_x[h_col_indices[idx]];
                double cand = (mv >= INF || xv >= INF) ? INF : (mv + xv);
                if (cand < best) best = cand;
            }
            y_trop[r] = best;
        }
        std::vector<double> ey;
        if (!load_engine_y(ey)) return 1;
        double maxd = 0.0; long mism = 0;
        for (int r = 0; r < rows; ++r) {
            double d = std::abs(y_trop[r] - ey[r]);
            if (d > maxd) maxd = d;
            if (d > 0.5) ++mism;
        }
        std::cout << "CROSSCHECK cpu_csr_minplus vs engine: max_abs_diff=" << maxd
                  << ", mismatches=" << mism
                  << (mism == 0 ? "  SUCCESS" : "  FAILURE") << std::endl;
        return 0;
    }

    std::vector<float> h_y(rows, 0.0f);
    std::vector<float> h_y_ref(rows, 0.0f);

    // Run CPU sequential reference for verification
    std::cout << "Running CPU reference SpMV..." << std::endl;
    for (int r = 0; r < rows; ++r) {
        double sum = 0.0;
        int row_start = h_row_ptr[r];
        int row_end = h_row_ptr[r + 1];
        for (int idx = row_start; idx < row_end; ++idx) {
            sum += (double)h_values[idx] * (double)h_x[h_col_indices[idx]];
        }
        h_y_ref[r] = (float)sum;
    }

    // 5. Device memory allocations
    int* d_row_ptr = nullptr;
    int* d_col_indices = nullptr;
    float* d_values = nullptr;
    float* d_x = nullptr;
    float* d_y = nullptr;

    // Baseline free memory before any of our allocations (apples-to-apples with
    // the grammar engine: both use unified memory on the GB10 coherent address space).
    size_t free0 = 0, total0 = 0;
    CHECK_CUDA(cudaMemGetInfo(&free0, &total0));

    // Unified memory (cudaMallocManaged) so the CSR baseline is measured under the
    // same allocator and coherent address space as the grammar engine.
    CHECK_CUDA(cudaMallocManaged(&d_row_ptr, h_row_ptr.size() * sizeof(int)));
    CHECK_CUDA(cudaMallocManaged(&d_col_indices, h_col_indices.size() * sizeof(int)));
    CHECK_CUDA(cudaMallocManaged(&d_values, h_values.size() * sizeof(float)));
    CHECK_CUDA(cudaMallocManaged(&d_x, h_x.size() * sizeof(float)));
    CHECK_CUDA(cudaMallocManaged(&d_y, h_y.size() * sizeof(float)));

    CHECK_CUDA(cudaMemcpy(d_row_ptr, h_row_ptr.data(), h_row_ptr.size() * sizeof(int), cudaMemcpyHostToDevice));
    CHECK_CUDA(cudaMemcpy(d_col_indices, h_col_indices.data(), h_col_indices.size() * sizeof(int), cudaMemcpyHostToDevice));
    CHECK_CUDA(cudaMemcpy(d_values, h_values.data(), h_values.size() * sizeof(float), cudaMemcpyHostToDevice));
    CHECK_CUDA(cudaMemcpy(d_x, h_x.data(), h_x.size() * sizeof(float), cudaMemcpyHostToDevice));
    CHECK_CUDA(cudaMemset(d_y, 0, h_y.size() * sizeof(float)));

    // 6. cuSPARSE Setup
    cusparseHandle_t handle;
    CHECK_CUSPARSE(cusparseCreate(&handle));

    cusparseSpMatDescr_t matA;
    CHECK_CUSPARSE(cusparseCreateCsr(&matA, rows, cols, nnz,
                                     d_row_ptr, d_col_indices, d_values,
                                     CUSPARSE_INDEX_32I, CUSPARSE_INDEX_32I,
                                     CUSPARSE_INDEX_BASE_ZERO, CUDA_R_32F));

    cusparseDnVecDescr_t vecX, vecY;
    CHECK_CUSPARSE(cusparseCreateDnVec(&vecX, cols, d_x, CUDA_R_32F));
    CHECK_CUSPARSE(cusparseCreateDnVec(&vecY, rows, d_y, CUDA_R_32F));

    float alpha_param = 1.0f;
    float beta_param = 0.0f;
    size_t bufferSize = 0;
    CHECK_CUSPARSE(cusparseSpMV_bufferSize(
        handle, CUSPARSE_OPERATION_NON_TRANSPOSE,
        &alpha_param, matA, vecX, &beta_param, vecY, CUDA_R_32F,
        CUSPARSE_SPMV_ALG_DEFAULT, &bufferSize));

    void* dBuffer = nullptr;
    CHECK_CUDA(cudaMallocManaged(&dBuffer, bufferSize));

    // Warm up
    CHECK_CUSPARSE(cusparseSpMV(
        handle, CUSPARSE_OPERATION_NON_TRANSPOSE,
        &alpha_param, matA, vecX, &beta_param, vecY, CUDA_R_32F,
        CUSPARSE_SPMV_ALG_DEFAULT, dBuffer));
    CHECK_CUDA(cudaDeviceSynchronize());

    // --- Memory accounting. Analytic CSR footprint: row pointers + (value,index)
    // per non-zero + the dense vectors + the cuSPARSE workspace. Measured peak:
    // resident device pages after the warm-up, minus the pre-allocation baseline.
    size_t mem_analytic_bytes =
          (size_t)h_row_ptr.size() * sizeof(int)        // d_row_ptr
        + (size_t)h_col_indices.size() * sizeof(int)    // d_col_indices
        + (size_t)h_values.size() * sizeof(float)       // d_values
        + (size_t)h_x.size() * sizeof(float)            // d_x
        + (size_t)h_y.size() * sizeof(float)            // d_y
        + (size_t)bufferSize;                           // cuSPARSE workspace
    size_t free1 = 0, total1 = 0;
    CHECK_CUDA(cudaMemGetInfo(&free1, &total1));
    size_t mem_peak_bytes = (free0 > free1) ? (free0 - free1) : 0;

    // Verify GPU results
    CHECK_CUDA(cudaMemcpy(h_y.data(), d_y, h_y.size() * sizeof(float), cudaMemcpyDeviceToHost));
    float max_abs_diff = 0.0f;
    float max_rel_diff = 0.0f;
    for (int r = 0; r < rows; ++r) {
        float diff = std::abs(h_y[r] - h_y_ref[r]);
        if (diff > max_abs_diff) max_abs_diff = diff;
        float rel = diff / (1.0f + std::abs(h_y_ref[r]));
        if (rel > max_rel_diff) max_rel_diff = rel;
    }
    std::cout << "Verification: max absolute diff = " << max_abs_diff 
              << ", max relative diff = " << max_rel_diff << std::endl;
    if (max_rel_diff > 1e-4f) {
        std::cerr << "WARNING: Verification mismatch!" << std::endl;
    } else {
        std::cout << "SUCCESS: Verification passed!" << std::endl;
    }

    // CROSSCHECK: assert both this CSR reconstruction (CPU SpMV) and the cuSPARSE GPU
    // kernel agree with the engine's dumped y on the shared vector. For Boolean the
    // engine emits 0/1 (OR-AND); the +x count agrees iff (count != 0) matches. Then stop.
    if (cc && cc[0]) {
        std::vector<double> ey;
        if (load_engine_y(ey)) {
            if (sr == "boolean") {
                long mism_cpu = 0, mism_gpu = 0;
                for (int r = 0; r < rows; ++r) {
                    int eb = (ey[r] != 0.0) ? 1 : 0;
                    if (((h_y_ref[r] != 0.0f) ? 1 : 0) != eb) ++mism_cpu;
                    if (((h_y[r]     != 0.0f) ? 1 : 0) != eb) ++mism_gpu;
                }
                std::cout << "CROSSCHECK cpu_csr_spmv(bool) vs engine: mismatches=" << mism_cpu
                          << (mism_cpu == 0 ? "  SUCCESS" : "  FAILURE") << std::endl;
                std::cout << "CROSSCHECK cusparse(bool) vs engine: mismatches=" << mism_gpu
                          << (mism_gpu == 0 ? "  SUCCESS" : "  FAILURE") << std::endl;
            } else {
                double rc = 0.0, rg = 0.0;
                for (int r = 0; r < rows; ++r) {
                    double m = std::abs(ey[r]);
                    double dc = std::abs((double)h_y_ref[r] - ey[r]) / (1.0 + m);
                    double dg = std::abs((double)h_y[r]     - ey[r]) / (1.0 + m);
                    if (dc > rc) rc = dc;
                    if (dg > rg) rg = dg;
                }
                std::cout << "CROSSCHECK cpu_csr_spmv vs engine: max_rel_diff=" << rc
                          << (rc < 1e-4 ? "  SUCCESS" : "  FAILURE") << std::endl;
                std::cout << "CROSSCHECK cusparse vs engine: max_rel_diff=" << rg
                          << (rg < 1e-4 ? "  SUCCESS" : "  FAILURE") << std::endl;
            }
        }
        // Skip the timing/energy/SpMM machinery in crosscheck mode.
        CHECK_CUDA(cudaFree(dBuffer));
        CHECK_CUSPARSE(cusparseDestroySpMat(matA));
        CHECK_CUSPARSE(cusparseDestroyDnVec(vecX));
        CHECK_CUSPARSE(cusparseDestroyDnVec(vecY));
        CHECK_CUSPARSE(cusparseDestroy(handle));
        CHECK_CUDA(cudaFree(d_row_ptr));
        CHECK_CUDA(cudaFree(d_col_indices));
        CHECK_CUDA(cudaFree(d_values));
        CHECK_CUDA(cudaFree(d_x));
        CHECK_CUDA(cudaFree(d_y));
        return 0;
    }

    // Benchmark loop
    std::cout << "Running " << iters << " iterations of cuSPARSE SpMV..." << std::endl;
    auto start = std::chrono::high_resolution_clock::now();
    for (int i = 0; i < iters; ++i) {
        CHECK_CUSPARSE(cusparseSpMV(
            handle, CUSPARSE_OPERATION_NON_TRANSPOSE,
            &alpha_param, matA, vecX, &beta_param, vecY, CUDA_R_32F,
            CUSPARSE_SPMV_ALG_DEFAULT, dBuffer));
    }
    CHECK_CUDA(cudaDeviceSynchronize());
    auto end = std::chrono::high_resolution_clock::now();

    std::chrono::duration<double, std::milli> duration = end - start;
    double avg_ms = duration.count() / iters;
    std::cout << "cuSPARSE SpMV completed. Average time: " << avg_ms << " ms/vector" << std::endl;

    // --- Energy per vector via the NVML energy counter, over a sustained loop
    // (>=1.5 s), AFTER the timed loop so the ms/vector above is unaffected.
    double mj_per_vector = -1.0, avg_power_mw = -1.0;
    {
        if (nvmlInit() == NVML_SUCCESS) {
            nvmlDevice_t nvd;
            if (nvmlDeviceGetHandleByIndex(0, &nvd) == NVML_SUCCESS) {
                unsigned long long e0 = 0, e1 = 0;
                nvmlReturn_t r0 = nvmlDeviceGetTotalEnergyConsumption(nvd, &e0);
                auto et0 = std::chrono::high_resolution_clock::now();
                long ecount = 0;
                double elapsed_s = 0.0;
                do {
                    for (int rep = 0; rep < 10; ++rep) {
                        cusparseSpMV(handle, CUSPARSE_OPERATION_NON_TRANSPOSE,
                                     &alpha_param, matA, vecX, &beta_param, vecY, CUDA_R_32F,
                                     CUSPARSE_SPMV_ALG_DEFAULT, dBuffer);
                        ++ecount;
                    }
                    CHECK_CUDA(cudaDeviceSynchronize());
                    elapsed_s = std::chrono::duration<double>(std::chrono::high_resolution_clock::now() - et0).count();
                } while (elapsed_s < 1.5);
                nvmlReturn_t r1 = nvmlDeviceGetTotalEnergyConsumption(nvd, &e1);
                if (r0 == NVML_SUCCESS && r1 == NVML_SUCCESS && ecount > 0) {
                    mj_per_vector = (double)(e1 - e0) / (double)ecount;
                    avg_power_mw  = (double)(e1 - e0) / elapsed_s; // mJ/s == mW
                }
            }
            nvmlShutdown();
        }
    }

    std::cout << "MEM analytic bytes: " << mem_analytic_bytes << std::endl;
    std::cout << "MEM peak bytes: " << mem_peak_bytes << std::endl;
    std::cout << "ENERGY mJ/vector: " << mj_per_vector << std::endl;
    std::cout << "ENERGY avg power mW: " << avg_power_mw << std::endl;

    // --- cuSPARSE SpMM (Y = M X over B vectors), apples-to-apples with the
    //     batched grammar engine: X is cols x B and Y is rows x B, row-major
    //     (ld = B), i.e. vectors innermost (X[j*B+b]), the same layout the
    //     engine uses. We replicate x across all B columns so each column of Y
    //     must equal the single-vector reference -- a cheap correctness check.
    int B = (argc >= 7) ? std::stoi(argv[6]) : 1;
    if (B > 1) {
        float* d_Xb = nullptr; float* d_Yb = nullptr;
        CHECK_CUDA(cudaMallocManaged(&d_Xb, (size_t)cols * B * sizeof(float)));
        CHECK_CUDA(cudaMallocManaged(&d_Yb, (size_t)rows * B * sizeof(float)));
        for (int j = 0; j < cols; ++j)
            for (int b = 0; b < B; ++b)
                d_Xb[(size_t)j * B + b] = h_x[j];
        CHECK_CUDA(cudaMemset(d_Yb, 0, (size_t)rows * B * sizeof(float)));

        cusparseDnMatDescr_t matX, matY;
        CHECK_CUSPARSE(cusparseCreateDnMat(&matX, cols, B, B, d_Xb, CUDA_R_32F, CUSPARSE_ORDER_ROW));
        CHECK_CUSPARSE(cusparseCreateDnMat(&matY, rows, B, B, d_Yb, CUDA_R_32F, CUSPARSE_ORDER_ROW));

        struct AlgEntry { const char* name; cusparseSpMMAlg_t alg; };
        AlgEntry algs[] = {
            {"default", CUSPARSE_SPMM_ALG_DEFAULT},
            {"csr_alg2", CUSPARSE_SPMM_CSR_ALG2},
            {"csr_alg3", CUSPARSE_SPMM_CSR_ALG3},
        };
        for (AlgEntry& AE : algs) {
            size_t bufMM = 0;
            cusparseStatus_t st = cusparseSpMM_bufferSize(
                handle, CUSPARSE_OPERATION_NON_TRANSPOSE, CUSPARSE_OPERATION_NON_TRANSPOSE,
                &alpha_param, matA, matX, &beta_param, matY, CUDA_R_32F,
                AE.alg, &bufMM);
            if (st != CUSPARSE_STATUS_SUCCESS) {
                std::cout << "cuSPARSE SpMM (B=" << B << ", " << AE.name << "): UNSUPPORTED" << std::endl;
                continue;
            }
            void* dBufMM = nullptr;
            CHECK_CUDA(cudaMallocManaged(&dBufMM, bufMM ? bufMM : 1));
            CHECK_CUDA(cudaMemset(d_Yb, 0, (size_t)rows * B * sizeof(float)));
            st = cusparseSpMM(
                handle, CUSPARSE_OPERATION_NON_TRANSPOSE, CUSPARSE_OPERATION_NON_TRANSPOSE,
                &alpha_param, matA, matX, &beta_param, matY, CUDA_R_32F,
                AE.alg, dBufMM);
            if (st != CUSPARSE_STATUS_SUCCESS) {
                std::cout << "cuSPARSE SpMM (B=" << B << ", " << AE.name << "): UNSUPPORTED_RUN" << std::endl;
                CHECK_CUDA(cudaFree(dBufMM));
                continue;
            }
            CHECK_CUDA(cudaDeviceSynchronize());
            float mm_rel = 0.0f;
            int check_rows = rows < 2000 ? rows : 2000;
            for (int r = 0; r < check_rows; ++r)
                for (int b = 0; b < B; ++b) {
                    float rel = std::abs(d_Yb[(size_t)r * B + b] - h_y_ref[r]) / (1.0f + std::abs(h_y_ref[r]));
                    if (rel > mm_rel) mm_rel = rel;
                }
            auto smm0 = std::chrono::high_resolution_clock::now();
            for (int i = 0; i < iters; ++i) {
                CHECK_CUSPARSE(cusparseSpMM(
                    handle, CUSPARSE_OPERATION_NON_TRANSPOSE, CUSPARSE_OPERATION_NON_TRANSPOSE,
                    &alpha_param, matA, matX, &beta_param, matY, CUDA_R_32F,
                    AE.alg, dBufMM));
            }
            CHECK_CUDA(cudaDeviceSynchronize());
            auto smm1 = std::chrono::high_resolution_clock::now();
            double mm_batch = std::chrono::duration<double, std::milli>(smm1 - smm0).count() / iters;
            std::cout << "cuSPARSE SpMM (B=" << B << ", " << AE.name << "): " << mm_batch
                      << " ms/batch, " << (mm_batch / B) << " ms/vector  rel=" << mm_rel
                      << (mm_rel < 1e-4f ? " SUCCESS" : " FAILURE") << std::endl;
            CHECK_CUDA(cudaFree(dBufMM));
        }

        CHECK_CUSPARSE(cusparseDestroyDnMat(matX));
        CHECK_CUSPARSE(cusparseDestroyDnMat(matY));
        CHECK_CUDA(cudaFree(d_Xb));
        CHECK_CUDA(cudaFree(d_Yb));
    }

    // 7. Cleanup
    CHECK_CUDA(cudaFree(dBuffer));
    CHECK_CUSPARSE(cusparseDestroySpMat(matA));
    CHECK_CUSPARSE(cusparseDestroyDnVec(vecX));
    CHECK_CUSPARSE(cusparseDestroyDnVec(vecY));
    CHECK_CUSPARSE(cusparseDestroy(handle));

    CHECK_CUDA(cudaFree(d_row_ptr));
    CHECK_CUDA(cudaFree(d_col_indices));
    CHECK_CUDA(cudaFree(d_values));
    CHECK_CUDA(cudaFree(d_x));
    CHECK_CUDA(cudaFree(d_y));

    return 0;
}
