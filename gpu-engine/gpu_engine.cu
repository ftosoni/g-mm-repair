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
#include <iostream>
#include <vector>
#include <chrono>
#include <nvml.h>

// The CUDA kernels of the streaming, double-buffered, level-synchronous sweep
// (manuscript sec:partI-gpu, Listing 1). One kernel launch per level k reads the
// previous frontier and writes the next; the host swaps the two buffers so they
// alternate read-only/write-only, keeping the live set at O(max_k w_k) rather
// than O(|R|) (lem:liveness). Every node update is one branch-free fused
// multiply-add (init_T_terminals / eval_level_layered) and each root of the top
// sequence C is atomically added into y the instant its level is computed and
// then freed (emit_level / emit_terminals -- "emit on the spot"). All four
// kernels are templated on a semiring policy from gpu_engine.h. The *_batched
// twins below do the same B vectors at once (SpMM, Y = M X), vectors innermost.
// Buffers are cudaMallocManaged: on the GB10's unified coherent memory there is
// no explicit host<->device copy for the frontiers (sec:setup).

#define CUDA_CHECK(call) \
    do { \
        cudaError_t err = call; \
        if (err != cudaSuccess) { \
            std::cerr << "CUDA Error: " << cudaGetErrorString(err) \
                      << " at " << __FILE__ << ":" << __LINE__ << std::endl; \
            exit(1); \
        } \
    } while (0)

template <typename T>
__global__ void fill_array_kernel(T* arr, T val, long size) {
    long idx = (long)blockIdx.x * blockDim.x + threadIdx.x;
    if (idx < size) {
        arr[idx] = val;
    }
}

template <typename T>
void fill_array(T* arr, T val, long size) {
    if (size <= 0) return;
    int threads = 256;
    long blocks = (size + threads - 1) / threads;
    fill_array_kernel<<<blocks, threads>>>(arr, val, size);
}

// Kernel to pre-populate terminal values in T (persistent input, never freed).
// T[0] = S::zero() is the "zero terminal" used as the second child of unary nodes.
template <typename S>
__global__ void init_T_terminals(
    typename S::T* T,
    const float* V,
    const float* x,
    const int* val_indices,
    const int* col_indices,
    int alpha
) {
    int p = blockIdx.x * blockDim.x + threadIdx.x;
    if (p == 0) {
        T[0] = S::zero();
    } else if (p < alpha) {
        int val_idx = val_indices[p - 1];
        int col_idx = col_indices[p - 1];
        T[p] = S::leaf(V[val_idx], x[col_idx]);
    }
}

// Level evaluation kernel (double-buffered, proper layered).
// One branch-free fused multiply-add covers binary / pass-through / run-length.
template <typename S>
__global__ void eval_level_layered(
    const LayeredRule* rules,
    const typename S::T* W_read,
    typename S::T* W_write,
    const typename S::T* T,
    int num_rules
) {
    int idx = blockIdx.x * blockDim.x + threadIdx.x;
    if (idx < num_rules) {
        LayeredRule r = rules[idx];
        typename S::T left_val = (r.left >= 0) ? W_read[r.left] : T[-r.left - 1];
        typename S::T right_val = (r.right >= 0) ? W_read[r.right] : T[-r.right - 1];
        W_write[idx] = S::oplus(S::apply(left_val, r.coeff), right_val);
    }
}

// Emit-on-the-spot: add the just-computed level-k frontier values to their rows
// of y. Multiple occurrences may target the same row, so the add is atomic.
template <typename S>
__global__ void emit_level(
    const int* emit_pos,
    const int* emit_row,
    const typename S::T* W_cur,
    typename S::T* y,
    int num_emit
) {
    int e = blockIdx.x * blockDim.x + threadIdx.x;
    if (e < num_emit) {
        S::emit(&y[emit_row[e]], W_cur[emit_pos[e]]);
    }
}

// Terminal occurrences of C (level 0): add T[sym] to the corresponding row of y.
template <typename S>
__global__ void emit_terminals(
    const int* term_sym,
    const int* term_row,
    const typename S::T* T,
    typename S::T* y,
    int num_term
) {
    int e = blockIdx.x * blockDim.x + threadIdx.x;
    if (e < num_term) {
        S::emit(&y[term_row[e]], T[term_sym[e]]);
    }
}

template <typename S>
void run_gpu_right_multiply_layered_templated(
    int rows,
    int cols,
    int alpha,
    int total_layered_rules,
    const LayeredRule* h_rules,
    const int* h_val_indices,
    const int* h_col_indices,
    const float* h_V,
    int num_distinct_vals,
    const float* h_x,
    const GPUSchedule& h_sched,
    typename S::T* h_y,
    int iters
) {
    using T = typename S::T;

    std::cout << "[GPU Setup Debug] rows=" << rows << ", cols=" << cols << ", alpha=" << alpha
              << ", total_rules=" << total_layered_rules << ", distinct_vals=" << num_distinct_vals
              << ", sched_max_depth=" << h_sched.max_depth << ", sched_max_width=" << h_sched.max_width
              << ", emit_total=" << h_sched.emit_total << ", num_term_emit=" << h_sched.num_term_emit
              << std::endl;

    // 1. Allocate device memory
    T* d_T = nullptr;
    LayeredRule* d_rules = nullptr;
    int* d_val_indices = nullptr;
    int* d_col_indices = nullptr;
    float* d_V = nullptr;
    float* d_x = nullptr;
    T* d_y = nullptr;
    T* d_W1 = nullptr;
    T* d_W2 = nullptr;
    int* d_emit_pos = nullptr;
    int* d_emit_row = nullptr;
    int* d_term_sym = nullptr;
    int* d_term_row = nullptr;

    size_t mem_analytic_bytes =
          (size_t)alpha * sizeof(T)                           // d_T
        + (size_t)total_layered_rules * sizeof(LayeredRule)   // d_rules
        + (size_t)(alpha - 1) * sizeof(int)                   // d_val_indices
        + (size_t)(alpha - 1) * sizeof(int)                   // d_col_indices
        + (size_t)num_distinct_vals * sizeof(float)           // d_V
        + (size_t)cols * sizeof(float)                        // d_x
        + (size_t)rows * sizeof(T)                            // d_y
        + (size_t)h_sched.max_width * sizeof(T) * 2           // d_W1 + d_W2
        + (size_t)h_sched.emit_total * sizeof(int) * 2        // d_emit_pos + d_emit_row
        + (size_t)h_sched.num_term_emit * sizeof(int) * 2;    // d_term_sym + d_term_row

    size_t free0 = 0, total0 = 0;
    CUDA_CHECK(cudaMemGetInfo(&free0, &total0));

    CUDA_CHECK(cudaMallocManaged(&d_T, alpha * sizeof(T)));
    CUDA_CHECK(cudaMallocManaged(&d_rules, total_layered_rules * sizeof(LayeredRule)));
    CUDA_CHECK(cudaMallocManaged(&d_val_indices, (alpha - 1) * sizeof(int)));
    CUDA_CHECK(cudaMallocManaged(&d_col_indices, (alpha - 1) * sizeof(int)));
    CUDA_CHECK(cudaMallocManaged(&d_V, num_distinct_vals * sizeof(float)));
    CUDA_CHECK(cudaMallocManaged(&d_x, cols * sizeof(float)));
    CUDA_CHECK(cudaMallocManaged(&d_y, rows * sizeof(T)));

    CUDA_CHECK(cudaMallocManaged(&d_W1, h_sched.max_width * sizeof(T)));
    CUDA_CHECK(cudaMallocManaged(&d_W2, h_sched.max_width * sizeof(T)));

    if (h_sched.emit_total > 0) {
        CUDA_CHECK(cudaMallocManaged(&d_emit_pos, h_sched.emit_total * sizeof(int)));
        CUDA_CHECK(cudaMallocManaged(&d_emit_row, h_sched.emit_total * sizeof(int)));
    }
    if (h_sched.num_term_emit > 0) {
        CUDA_CHECK(cudaMallocManaged(&d_term_sym, h_sched.num_term_emit * sizeof(int)));
        CUDA_CHECK(cudaMallocManaged(&d_term_row, h_sched.num_term_emit * sizeof(int)));
    }

    // 2. Copy data host -> device
    CUDA_CHECK(cudaMemcpy(d_rules, h_rules, total_layered_rules * sizeof(LayeredRule), cudaMemcpyHostToDevice));
    CUDA_CHECK(cudaMemcpy(d_val_indices, h_val_indices, (alpha - 1) * sizeof(int), cudaMemcpyHostToDevice));
    CUDA_CHECK(cudaMemcpy(d_col_indices, h_col_indices, (alpha - 1) * sizeof(int), cudaMemcpyHostToDevice));
    CUDA_CHECK(cudaMemcpy(d_V, h_V, num_distinct_vals * sizeof(float), cudaMemcpyHostToDevice));
    CUDA_CHECK(cudaMemcpy(d_x, h_x, cols * sizeof(float), cudaMemcpyHostToDevice));
    if (h_sched.emit_total > 0) {
        CUDA_CHECK(cudaMemcpy(d_emit_pos, h_sched.emit_pos, h_sched.emit_total * sizeof(int), cudaMemcpyHostToDevice));
        CUDA_CHECK(cudaMemcpy(d_emit_row, h_sched.emit_row, h_sched.emit_total * sizeof(int), cudaMemcpyHostToDevice));
    }
    if (h_sched.num_term_emit > 0) {
        CUDA_CHECK(cudaMemcpy(d_term_sym, h_sched.term_emit_sym, h_sched.num_term_emit * sizeof(int), cudaMemcpyHostToDevice));
        CUDA_CHECK(cudaMemcpy(d_term_row, h_sched.term_emit_row, h_sched.num_term_emit * sizeof(int), cudaMemcpyHostToDevice));
    }

    int threads_per_block = 256;
    int blocks = (alpha + threads_per_block - 1) / threads_per_block;

    // 3. Time the per-vector kernel work over `iters` input vectors.
    cudaEvent_t start_event, stop_event;
    CUDA_CHECK(cudaEventCreate(&start_event));
    CUDA_CHECK(cudaEventCreate(&stop_event));
    CUDA_CHECK(cudaEventRecord(start_event, 0));

    for (int it = 0; it < iters; ++it) {
        fill_array<T>(d_y, S::zero(), rows);

        // Terminals into T, then the terminal occurrences of C into y.
        init_T_terminals<S><<<blocks, threads_per_block>>>(d_T, d_V, d_x, d_val_indices, d_col_indices, alpha);
        if (h_sched.num_term_emit > 0) {
            int tb = (h_sched.num_term_emit + threads_per_block - 1) / threads_per_block;
            emit_terminals<S><<<tb, threads_per_block>>>(d_term_sym, d_term_row, d_T, d_y, h_sched.num_term_emit);
        }

        // Double-buffered level-by-level sweep with emit-on-the-spot.
        T* d_W_read = d_W1;
        T* d_W_write = d_W2;
        fill_array<T>(d_W_read, S::zero(), h_sched.max_width);
        fill_array<T>(d_W_write, S::zero(), h_sched.max_width);

        for (int k = 1; k <= h_sched.max_depth; ++k) {
            int start = h_sched.level_offsets[k - 1];
            int end = h_sched.level_offsets[k];
            int num_rules = end - start;
            if (num_rules > 0) {
                int eval_blocks = (num_rules + threads_per_block - 1) / threads_per_block;
                eval_level_layered<S><<<eval_blocks, threads_per_block>>>(
                    d_rules + start, d_W_read, d_W_write, d_T, num_rules
                );
            }
            int estart = h_sched.emit_offsets[k - 1];
            int eend = h_sched.emit_offsets[k];
            int num_emit = eend - estart;
            if (num_emit > 0) {
                int emit_blocks = (num_emit + threads_per_block - 1) / threads_per_block;
                emit_level<S><<<emit_blocks, threads_per_block>>>(
                    d_emit_pos + estart, d_emit_row + estart, d_W_write, d_y, num_emit
                );
            }
            T* temp = d_W_read;
            d_W_read = d_W_write;
            d_W_write = temp;
        }
    }

    CUDA_CHECK(cudaEventRecord(stop_event, 0));
    CUDA_CHECK(cudaEventSynchronize(stop_event));

    float total_ms = 0;
    CUDA_CHECK(cudaEventElapsedTime(&total_ms, start_event, stop_event));
    std::cout << "GPU kernel: " << total_ms << " ms total over " << iters
              << " vectors = " << (total_ms / iters) << " ms/vector." << std::endl;

    size_t free1 = 0, total1 = 0;
    CUDA_CHECK(cudaMemGetInfo(&free1, &total1));
    size_t mem_peak_bytes = (free0 > free1) ? (free0 - free1) : 0;

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
                        fill_array<T>(d_y, S::zero(), rows);
                        init_T_terminals<S><<<blocks, threads_per_block>>>(d_T, d_V, d_x, d_val_indices, d_col_indices, alpha);
                        if (h_sched.num_term_emit > 0) {
                            int tb = (h_sched.num_term_emit + threads_per_block - 1) / threads_per_block;
                            emit_terminals<S><<<tb, threads_per_block>>>(d_term_sym, d_term_row, d_T, d_y, h_sched.num_term_emit);
                        }
                        T* rd = d_W1; T* wr = d_W2;
                        fill_array<T>(rd, S::zero(), h_sched.max_width);
                        fill_array<T>(wr, S::zero(), h_sched.max_width);
                        for (int k = 1; k <= h_sched.max_depth; ++k) {
                            int s = h_sched.level_offsets[k - 1], e = h_sched.level_offsets[k];
                            int nr = e - s;
                            if (nr > 0) {
                                int eb = (nr + threads_per_block - 1) / threads_per_block;
                                eval_level_layered<S><<<eb, threads_per_block>>>(d_rules + s, rd, wr, d_T, nr);
                            }
                            int es = h_sched.emit_offsets[k - 1], ee = h_sched.emit_offsets[k];
                            int ne = ee - es;
                            if (ne > 0) {
                                int eb = (ne + threads_per_block - 1) / threads_per_block;
                                emit_level<S><<<eb, threads_per_block>>>(d_emit_pos + es, d_emit_row + es, wr, d_y, ne);
                            }
                            T* tmp = rd; rd = wr; wr = tmp;
                        }
                        ++ecount;
                    }
                    CUDA_CHECK(cudaDeviceSynchronize());
                    elapsed_s = std::chrono::duration<double>(std::chrono::high_resolution_clock::now() - et0).count();
                } while (elapsed_s < 1.5);
                nvmlReturn_t r1 = nvmlDeviceGetTotalEnergyConsumption(nvd, &e1);
                if (r0 == NVML_SUCCESS && r1 == NVML_SUCCESS && ecount > 0) {
                    mj_per_vector = (double)(e1 - e0) / (double)ecount;
                    avg_power_mw  = (double)(e1 - e0) / elapsed_s;
                }
            }
            nvmlShutdown();
        }
    }

    std::cout << "MEM analytic bytes: " << mem_analytic_bytes << std::endl;
    std::cout << "MEM peak bytes: " << mem_peak_bytes << std::endl;
    std::cout << "ENERGY mJ/vector: " << mj_per_vector << std::endl;
    std::cout << "ENERGY avg power mW: " << avg_power_mw << std::endl;

    CUDA_CHECK(cudaEventDestroy(start_event));
    CUDA_CHECK(cudaEventDestroy(stop_event));

    // 5. Copy output back to host
    CUDA_CHECK(cudaMemcpy(h_y, d_y, rows * sizeof(T), cudaMemcpyDeviceToHost));

    // 6. Free device allocations
    CUDA_CHECK(cudaFree(d_T));
    CUDA_CHECK(cudaFree(d_rules));
    CUDA_CHECK(cudaFree(d_val_indices));
    CUDA_CHECK(cudaFree(d_col_indices));
    CUDA_CHECK(cudaFree(d_V));
    CUDA_CHECK(cudaFree(d_x));
    CUDA_CHECK(cudaFree(d_y));
    CUDA_CHECK(cudaFree(d_W1));
    CUDA_CHECK(cudaFree(d_W2));
    if (d_emit_pos) CUDA_CHECK(cudaFree(d_emit_pos));
    if (d_emit_row) CUDA_CHECK(cudaFree(d_emit_row));
    if (d_term_sym) CUDA_CHECK(cudaFree(d_term_sym));
    if (d_term_row) CUDA_CHECK(cudaFree(d_term_row));
}

extern "C" void run_gpu_right_multiply_layered(
    int rows,
    int cols,
    int alpha,
    int total_layered_rules,
    const LayeredRule* h_rules,
    const int* h_val_indices,
    const int* h_col_indices,
    const float* h_V,
    int num_distinct_vals,
    const float* h_x,
    const GPUSchedule& h_sched,
    float* h_y,
    int iters
) {
    const char* semiring_env = getenv("SEMIRING");
    std::string semiring = semiring_env ? semiring_env : "plustimes";
    if (semiring == "boolean" || semiring == "bool") {
        run_gpu_right_multiply_layered_templated<Boolean>(
            rows, cols, alpha, total_layered_rules, h_rules, h_val_indices, h_col_indices,
            h_V, num_distinct_vals, h_x, h_sched, reinterpret_cast<int*>(h_y), iters
        );
    } else if (semiring == "tropical" || semiring == "minplus") {
        run_gpu_right_multiply_layered_templated<Tropical>(
            rows, cols, alpha, total_layered_rules, h_rules, h_val_indices, h_col_indices,
            h_V, num_distinct_vals, h_x, h_sched, reinterpret_cast<int*>(h_y), iters
        );
    } else {
        run_gpu_right_multiply_layered_templated<PlusTimes>(
            rows, cols, alpha, total_layered_rules, h_rules, h_val_indices, h_col_indices,
            h_V, num_distinct_vals, h_x, h_sched, h_y, iters
        );
    }
}

// =====================================================================
// Batched (SpMM) variant: B vectors at once.
// =====================================================================

template <typename S>
__global__ void init_T_batched(
    typename S::T* T, const float* V, const float* x,
    const int* val_indices, const int* col_indices, int alpha, int B
) {
    long tid = (long)blockIdx.x * blockDim.x + threadIdx.x;
    long total = (long)alpha * B;
    if (tid >= total) return;
    int p = (int)(tid / B);
    int b = (int)(tid % B);
    if (p == 0) { T[b] = S::zero(); return; }
    int v = val_indices[p - 1];
    int c = col_indices[p - 1];
    T[(long)p * B + b] = S::leaf(V[v], x[(long)c * B + b]);
}

template <typename S>
__global__ void eval_level_batched(
    const LayeredRule* rules, const typename S::T* W_read, typename S::T* W_write,
    const typename S::T* T, int num_rules, int B
) {
    long tid = (long)blockIdx.x * blockDim.x + threadIdx.x;
    long total = (long)num_rules * B;
    if (tid >= total) return;
    int idx = (int)(tid / B);
    int b = (int)(tid % B);
    LayeredRule r = rules[idx];
    typename S::T lv = (r.left  >= 0) ? W_read[(long)r.left  * B + b] : T[(long)(-r.left  - 1) * B + b];
    typename S::T rv = (r.right >= 0) ? W_read[(long)r.right * B + b] : T[(long)(-r.right - 1) * B + b];
    W_write[(long)idx * B + b] = S::oplus(S::apply(lv, r.coeff), rv);
}

template <typename S>
__global__ void emit_level_batched(
    const int* emit_pos, const int* emit_row, const typename S::T* W_cur,
    typename S::T* y, int num_emit, int B
) {
    long tid = (long)blockIdx.x * blockDim.x + threadIdx.x;
    long total = (long)num_emit * B;
    if (tid >= total) return;
    int e = (int)(tid / B);
    int b = (int)(tid % B);
    S::emit(&y[(long)emit_row[e] * B + b], W_cur[(long)emit_pos[e] * B + b]);
}

template <typename S>
__global__ void emit_terminals_batched(
    const int* term_sym, const int* term_row, const typename S::T* T,
    typename S::T* y, int num_term, int B
) {
    long tid = (long)blockIdx.x * blockDim.x + threadIdx.x;
    long total = (long)num_term * B;
    if (tid >= total) return;
    int e = (int)(tid / B);
    int b = (int)(tid % B);
    S::emit(&y[(long)term_row[e] * B + b], T[(long)term_sym[e] * B + b]);
}

template <typename S>
void run_gpu_right_multiply_batched_templated(
    int rows, int cols, int alpha, int total_layered_rules,
    const LayeredRule* h_rules, const int* h_val_indices, const int* h_col_indices,
    const float* h_V, int num_distinct_vals, const float* h_X,
    const GPUSchedule& h_sched, typename S::T* h_Y, int B, int iters
) {
    using T = typename S::T;

    T* d_T = nullptr; LayeredRule* d_rules = nullptr;
    int* d_val_indices = nullptr; int* d_col_indices = nullptr;
    float* d_V = nullptr; float* d_X = nullptr; T* d_Y = nullptr;
    T* d_W1 = nullptr; T* d_W2 = nullptr;
    int* d_emit_pos = nullptr; int* d_emit_row = nullptr;
    int* d_term_sym = nullptr; int* d_term_row = nullptr;

    CUDA_CHECK(cudaMallocManaged(&d_T, (long)alpha * B * sizeof(T)));
    CUDA_CHECK(cudaMallocManaged(&d_rules, total_layered_rules * sizeof(LayeredRule)));
    CUDA_CHECK(cudaMallocManaged(&d_val_indices, (alpha - 1) * sizeof(int)));
    CUDA_CHECK(cudaMallocManaged(&d_col_indices, (alpha - 1) * sizeof(int)));
    CUDA_CHECK(cudaMallocManaged(&d_V, num_distinct_vals * sizeof(float)));
    CUDA_CHECK(cudaMallocManaged(&d_X, (long)cols * B * sizeof(float)));
    CUDA_CHECK(cudaMallocManaged(&d_Y, (long)rows * B * sizeof(T)));
    CUDA_CHECK(cudaMallocManaged(&d_W1, (long)h_sched.max_width * B * sizeof(T)));
    CUDA_CHECK(cudaMallocManaged(&d_W2, (long)h_sched.max_width * B * sizeof(T)));
    if (h_sched.emit_total > 0) {
        CUDA_CHECK(cudaMallocManaged(&d_emit_pos, h_sched.emit_total * sizeof(int)));
        CUDA_CHECK(cudaMallocManaged(&d_emit_row, h_sched.emit_total * sizeof(int)));
    }
    if (h_sched.num_term_emit > 0) {
        CUDA_CHECK(cudaMallocManaged(&d_term_sym, h_sched.num_term_emit * sizeof(int)));
        CUDA_CHECK(cudaMallocManaged(&d_term_row, h_sched.num_term_emit * sizeof(int)));
    }

    CUDA_CHECK(cudaMemcpy(d_rules, h_rules, total_layered_rules * sizeof(LayeredRule), cudaMemcpyHostToDevice));
    CUDA_CHECK(cudaMemcpy(d_val_indices, h_val_indices, (alpha - 1) * sizeof(int), cudaMemcpyHostToDevice));
    CUDA_CHECK(cudaMemcpy(d_col_indices, h_col_indices, (alpha - 1) * sizeof(int), cudaMemcpyHostToDevice));
    CUDA_CHECK(cudaMemcpy(d_V, h_V, num_distinct_vals * sizeof(float), cudaMemcpyHostToDevice));
    CUDA_CHECK(cudaMemcpy(d_X, h_X, (long)cols * B * sizeof(float), cudaMemcpyHostToDevice));
    if (h_sched.emit_total > 0) {
        CUDA_CHECK(cudaMemcpy(d_emit_pos, h_sched.emit_pos, h_sched.emit_total * sizeof(int), cudaMemcpyHostToDevice));
        CUDA_CHECK(cudaMemcpy(d_emit_row, h_sched.emit_row, h_sched.emit_total * sizeof(int), cudaMemcpyHostToDevice));
    }
    if (h_sched.num_term_emit > 0) {
        CUDA_CHECK(cudaMemcpy(d_term_sym, h_sched.term_emit_sym, h_sched.num_term_emit * sizeof(int), cudaMemcpyHostToDevice));
        CUDA_CHECK(cudaMemcpy(d_term_row, h_sched.term_emit_row, h_sched.num_term_emit * sizeof(int), cudaMemcpyHostToDevice));
    }

    const int TPB = 256;
    auto grid = [&](long work) { return (int)((work + TPB - 1) / TPB); };

    cudaEvent_t start_event, stop_event;
    CUDA_CHECK(cudaEventCreate(&start_event));
    CUDA_CHECK(cudaEventCreate(&stop_event));
    CUDA_CHECK(cudaEventRecord(start_event, 0));

    for (int it = 0; it < iters; ++it) {
        fill_array<T>(d_Y, S::zero(), (long)rows * B);
        init_T_batched<S><<<grid((long)alpha * B), TPB>>>(d_T, d_V, d_X, d_val_indices, d_col_indices, alpha, B);
        if (h_sched.num_term_emit > 0)
            emit_terminals_batched<S><<<grid((long)h_sched.num_term_emit * B), TPB>>>(d_term_sym, d_term_row, d_T, d_Y, h_sched.num_term_emit, B);

        T* rd = d_W1; T* wr = d_W2;
        fill_array<T>(rd, S::zero(), (long)h_sched.max_width * B);
        fill_array<T>(wr, S::zero(), (long)h_sched.max_width * B);

        for (int k = 1; k <= h_sched.max_depth; ++k) {
            int s = h_sched.level_offsets[k - 1], e = h_sched.level_offsets[k];
            int nr = e - s;
            if (nr > 0)
                eval_level_batched<S><<<grid((long)nr * B), TPB>>>(d_rules + s, rd, wr, d_T, nr, B);
            int es = h_sched.emit_offsets[k - 1], ee = h_sched.emit_offsets[k];
            int ne = ee - es;
            if (ne > 0)
                emit_level_batched<S><<<grid((long)ne * B), TPB>>>(d_emit_pos + es, d_emit_row + es, wr, d_Y, ne, B);
            T* tmp = rd; rd = wr; wr = tmp;
        }
    }

    CUDA_CHECK(cudaEventRecord(stop_event, 0));
    CUDA_CHECK(cudaEventSynchronize(stop_event));
    float total_ms = 0;
    CUDA_CHECK(cudaEventElapsedTime(&total_ms, start_event, stop_event));
    double per_batch = total_ms / iters;
    double per_vec = per_batch / B;
    std::cout << "GPU batched (B=" << B << "): " << per_batch << " ms/batch, "
              << per_vec << " ms/vector over " << iters << " batches." << std::endl;
    CUDA_CHECK(cudaEventDestroy(start_event));
    CUDA_CHECK(cudaEventDestroy(stop_event));

    CUDA_CHECK(cudaMemcpy(h_Y, d_Y, (long)rows * B * sizeof(T), cudaMemcpyDeviceToHost));

    CUDA_CHECK(cudaFree(d_T)); CUDA_CHECK(cudaFree(d_rules));
    CUDA_CHECK(cudaFree(d_val_indices)); CUDA_CHECK(cudaFree(d_col_indices));
    CUDA_CHECK(cudaFree(d_V)); CUDA_CHECK(cudaFree(d_X)); CUDA_CHECK(cudaFree(d_Y));
    CUDA_CHECK(cudaFree(d_W1)); CUDA_CHECK(cudaFree(d_W2));
    if (d_emit_pos) CUDA_CHECK(cudaFree(d_emit_pos));
    if (d_emit_row) CUDA_CHECK(cudaFree(d_emit_row));
    if (d_term_sym) CUDA_CHECK(cudaFree(d_term_sym));
    if (d_term_row) CUDA_CHECK(cudaFree(d_term_row));
}

extern "C" void run_gpu_right_multiply_batched(
    int rows, int cols, int alpha, int total_layered_rules,
    const LayeredRule* h_rules, const int* h_val_indices, const int* h_col_indices,
    const float* h_V, int num_distinct_vals, const float* h_X,
    const GPUSchedule& h_sched, float* h_Y, int B, int iters
) {
    const char* semiring_env = getenv("SEMIRING");
    std::string semiring = semiring_env ? semiring_env : "plustimes";
    if (semiring == "boolean" || semiring == "bool") {
        run_gpu_right_multiply_batched_templated<Boolean>(
            rows, cols, alpha, total_layered_rules, h_rules, h_val_indices, h_col_indices,
            h_V, num_distinct_vals, h_X, h_sched, reinterpret_cast<int*>(h_Y), B, iters
        );
    } else if (semiring == "tropical" || semiring == "minplus") {
        run_gpu_right_multiply_batched_templated<Tropical>(
            rows, cols, alpha, total_layered_rules, h_rules, h_val_indices, h_col_indices,
            h_V, num_distinct_vals, h_X, h_sched, reinterpret_cast<int*>(h_Y), B, iters
        );
    } else {
        run_gpu_right_multiply_batched_templated<PlusTimes>(
            rows, cols, alpha, total_layered_rules, h_rules, h_val_indices, h_col_indices,
            h_V, num_distinct_vals, h_X, h_sched, h_Y, B, iters
        );
    }
}