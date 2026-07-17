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

#ifndef LAYERED_TYPES_H
#define LAYERED_TYPES_H

// Plain (CUDA-free) definitions of the proper-layered schedule that the paper's
// streaming, double-buffered sweep consumes (manuscript sec:partI, Listing 1).
// gpu_engine.h includes this; the host-side schedule builder (build_schedule.*)
// also includes it, so it compiles with a plain C++ host compiler -- no CUDA
// toolchain required for the grammar analysis the paper describes as "computed
// once on the host".

// Layered rule for the proper-layered, double-buffered sweep. This is the
// branch-free unified node update of the manuscript (sec:partI-gpu, eq. for
// cur = coeff * cv(left) + cv(right)); expressed here for the (+,x) semiring,
// with the generic (leaf/combine) version templated in gpu_engine.h.
// Each node computes  W[idx] = coeff * cv(left) + cv(right), where
//   cv(c) = W_read[c]      if c >= 0   (a slot of the previous frontier)
//         = T[-c - 1]      if c <  0   (a persistent terminal; T[0] = 0)
// Encoding of the three kinds (a "zero" child is the terminal index for T[0]):
//   binary       N -> A B : left=A, right=B,    coeff=1
//   pass-through N -> M    : left=M, right=zero, coeff=1   (identity node added
//                                                           by completion)
//   run-length   N -> B^t  : left=B, right=zero, coeff=t   (eq:rleval)
// Pass-through nodes are inserted by pass-through completion (sec:partI-layer)
// and are the only structural overhead the paper reports (+pt). The run-length
// kind is a reserved capability the sweep can evaluate; the RePair path that the
// paper evaluates produces binary and pass-through rules only.
struct LayeredRule {
    int left;
    int right;
    float coeff;
};

// Level-synchronous schedule plus the emit-on-the-spot lists.
//
// Rules of level k occupy [level_offsets[k-1], level_offsets[k]) in the
// flattened rule array (k = 1..max_depth).
//
// A symbol occurring in the top sequence C is emitted into y at the level where
// it is computed and then freed -- it is NOT carried to the top. Non-terminal
// emissions of level k occupy [emit_offsets[k-1], emit_offsets[k]) and add
// frontier_k[emit_pos[e]] to y[emit_row[e]]. Terminal occurrences of C add
// T[term_emit_sym[e]] to y[term_emit_row[e]] once, before the sweep.
struct GPUSchedule {
    int max_depth;
    int max_width;
    int* level_offsets;      // size max_depth + 1
    int* level_widths;       // size max_depth + 1

    int* emit_offsets;       // size max_depth + 1
    int* emit_pos;           // size emit_total
    int* emit_row;           // size emit_total
    int  emit_total;

    int* term_emit_sym;      // size num_term_emit (index into T)
    int* term_emit_row;      // size num_term_emit
    int  num_term_emit;
};

#endif // LAYERED_TYPES_H