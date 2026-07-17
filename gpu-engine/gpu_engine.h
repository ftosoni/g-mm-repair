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

#ifndef GPU_ENGINE_H
#define GPU_ENGINE_H

#include <cuda_runtime.h>
#include "layered_types.h"   // LayeredRule, GPUSchedule (CUDA-free, shared)

// Semiring policy structs -- the "concise leaf/combine/emit policy" of the
// manuscript's monoid-homomorphism engine (section "Beyond (+,x): a monoid-homomorphism engine"). The kernels below are
// templated on one of these and selected at runtime by the SEMIRING env var, so
// the *same* schedule and sweep evaluate any monoid over the grammar; only these
// primitives change. Correctness needs only that oplus be associative (a
// monoid); a skipped/zero entry must equal zero() (the oplus-identity).
// Each policy provides:
//   T          value type carried through the sweep
//   zero()     additive identity of oplus (also the value of the zero terminal
//              T[0] and the initial y); see gotcha (1): Tropical needs a fill,
//              not memset
//   oplus(a,b) the associative combine (the level update's reduction)
//   leaf(V,x)  the leaf/multiplicative map applied at terminals: eval of one
//              matrix entry
//   apply(v,c) the coeff power used by binary (c=1) and run-length (c=t) rules;
//              a plain scalar multiply for (+,x), a no-op for idempotent oplus
//   emit/cpu_emit  the atomic scatter of a finished root into y (emit-on-the-spot)
// Mapping (manuscript section "Beyond (+,x): a monoid-homomorphism engine"):
//   PlusTimes -> SpMV;  Boolean -> BFS frontier/reachability;  Tropical ->
//   Bellman-Ford relaxation (SSSP; APSP when batched).

struct PlusTimes {
    using T = float;
    static __device__ __host__ T zero() { return 0.0f; }
    static __device__ __host__ T oplus(T a, T b) { return a + b; }
    static __device__ __host__ T leaf(float Vv, float xv) { return Vv * xv; }
    static __device__ __host__ T apply(T v, float c) { return c * v; }
    static __device__ void emit(T* y, T v) { atomicAdd(y, v); }
    static void cpu_emit(T* y, T v) {
        #pragma omp atomic
        *y += v;
    }
};

struct Boolean {
    using T = int;
    static __device__ __host__ T zero() { return 0; }
    static __device__ __host__ T oplus(T a, T b) { return a | b; }
    static __device__ __host__ T leaf(float Vv, float xv) { return (Vv != 0.0f && xv != 0.0f) ? 1 : 0; }
    static __device__ __host__ T apply(T v, float c) { return v; }
    static __device__ void emit(T* y, T v) { atomicOr(y, v); }
    static void cpu_emit(T* y, T v) {
        #pragma omp atomic
        *y |= v;
    }
};

// Tropical (min,+). Uses INTEGER distances so atomicMin is native (gotcha 1),
// with 1e9 as the +inf sentinel (cannot memset +inf, so zero() returns it and a
// fill kernel initialises buffers). A stored weight-0 edge is a real edge; a
// skipped entry is "no edge" = +inf (gotcha 2), so graphs are encoded as each
// source's out-edges <weight,target>.
struct Tropical {
    using T = int;
    static __device__ __host__ T zero() { return 1000000000; }
    static __device__ __host__ T oplus(T a, T b) { return a < b ? a : b; }
    static __device__ __host__ T leaf(float Vv, float xv) {
        int iv = (int)Vv;
        int ix = (int)xv;
        if (iv >= 1000000000 || ix >= 1000000000) return 1000000000;
        return iv + ix;
    }
    static __device__ __host__ T apply(T v, float c) { return v; }
    static __device__ void emit(T* y, T v) { atomicMin(y, v); }
    static void cpu_emit(T* y, T v) {
        #pragma omp critical(tropical_emit)
        {
            if (v < *y) {
                *y = v;
            }
        }
    }
};

// Main GPU execution function (right multiplication y = M x).
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
);

// Batched right multiplication Y = M X over B vectors at once (SpMM). The
// grammar structure is read once per level and applied B-wide, amortising the
// per-vector traversal overhead. Vectors are the innermost dimension: X is
// column-major-by-vector with x[j*B + b], Y is y[r*B + b]. Returns the average
// time per (batch) and per (vector) via stdout. h_Y must hold rows*B floats.
extern "C" void run_gpu_right_multiply_batched(
    int rows,
    int cols,
    int alpha,
    int total_layered_rules,
    const LayeredRule* h_rules,
    const int* h_val_indices,
    const int* h_col_indices,
    const float* h_V,
    int num_distinct_vals,
    const float* h_X,        // cols * B
    const GPUSchedule& h_sched,
    float* h_Y,              // rows * B
    int B,
    int iters
);

#endif // GPU_ENGINE_H