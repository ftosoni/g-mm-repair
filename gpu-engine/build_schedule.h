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

#ifndef BUILD_SCHEDULE_H
#define BUILD_SCHEDULE_H

#include <vector>
#include <utility>
#include "grammar.h"
#include "layered_types.h"

// Turns the RePair grammar into the proper-layered, double-buffered schedule the
// GPU engine consumes. This is "Step 2: pass-through completion" of the
// manuscript (section "A proper-layered streaming engine"): the single host-side step the paper describes
// as an offline, amortized preprocessing cost. It is CUDA-free and runs on the
// host. Levels are computed from the rule structure (the manuscript's level equation,
// extended for run-length rules: lvl(N->B^t) = 1 + lvl(B)); proper layering
// (the "Proper layering" definition) is then enforced by inserting identity pass-through nodes on
// every nonterminal edge that skips a level. The number of inserted nodes is
// the +pt inflation reported in the paper's structural tables (the genotype and
// Wikidata ones); it is small on a natively layered grammar and grows with
// the edge-span profile of a RePair grammar.
struct BuiltSchedule {
    int alpha = 0;
    int rows = 0;
    int cols = 0;
    int max_depth = 0;
    int max_width = 0;
    std::size_t num_raw_nt = 0;        // grammar non-terminals (before completion)
    std::size_t num_layered = 0;       // after completion (raw + pass-throughs)
    std::size_t num_passthrough = 0;   // num_layered - num_raw_nt  (+pt)

    // GPU-facing flattened arrays (owned by this struct).
    std::vector<LayeredRule> flattened_rules;
    std::vector<int> level_offsets;    // size max_depth + 1
    std::vector<int> level_widths;     // size max_depth + 1
    std::vector<int> emit_pos;
    std::vector<int> emit_row;
    std::vector<int> emit_offsets;     // size max_depth + 1
    std::vector<int> term_emit_sym;
    std::vector<int> term_emit_row;

    // Terminal compaction: after layout, only terminals actually referenced are
    // kept, remapped to a contiguous [0,alpha) range (alpha is updated to the
    // compacted count U). term_orig[q] is the ORIGINAL dense (value,column)
    // terminal id of compact id q (term_orig[0] = 0, the zero terminal), so the
    // caller can rebuild val/col indices for the U kept terminals only.
    std::vector<int> term_orig;

    // CPU-reference nested structures (indexed by level 1..max_depth; [0] empty).
    std::vector<std::vector<LayeredRule>> layered_levels;
    std::vector<std::vector<std::pair<int,int>>> emit_by_level; // (pos, row)
    std::vector<std::pair<int,int>> term_emit;                  // (T-index, row)

    // View over the flattened arrays for the GPU entry point. The returned
    // schedule points into this struct, which must outlive it.
    GPUSchedule gpu_schedule();
};

// rows/cols are the matrix dimensions (cols needed only to keep the schedule's
// row bookkeeping consistent; terminal decoding happens in the engine).
BuiltSchedule build_schedule(const grammar::Grammar& g, int rows, int cols);

#endif // BUILD_SCHEDULE_H