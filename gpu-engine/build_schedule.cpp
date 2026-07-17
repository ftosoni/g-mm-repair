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

#include "build_schedule.h"
#include <algorithm>

using grammar::Grammar;
using grammar::Rule;
using grammar::RuleKind;

GPUSchedule BuiltSchedule::gpu_schedule() {
    GPUSchedule s;
    s.max_depth     = max_depth;
    s.max_width     = max_width;
    s.level_offsets = level_offsets.data();
    s.level_widths  = level_widths.data();
    s.emit_offsets  = emit_offsets.data();
    s.emit_pos      = emit_pos.data();
    s.emit_row      = emit_row.data();
    s.emit_total    = static_cast<int>(emit_pos.size());
    s.term_emit_sym = term_emit_sym.data();
    s.term_emit_row = term_emit_row.data();
    s.num_term_emit = static_cast<int>(term_emit_sym.size());
    return s;
}

BuiltSchedule build_schedule(const Grammar& g, int rows, int cols) {
    BuiltSchedule out;
    out.alpha = g.alpha;
    out.rows = rows;
    out.cols = cols;

    const int alpha = g.alpha;
    const std::size_t nt_num = g.rules.size();
    out.num_raw_nt = nt_num;

    // Level of a child symbol: terminals (incl. the $ delimiter 0) sit at 0.
    auto child_lvl = [&](int sym, const std::vector<int>& lvl) -> int {
        return (sym >= alpha) ? lvl[sym - alpha] : 0;
    };

    // --- (3) Levels of the raw rules, eq. (3) extended for run-length. ---
    std::vector<int> lvl(nt_num, 0);
    int max_depth = 0;
    for (std::size_t i = 0; i < nt_num; ++i) {
        const Rule& r = g.rules[i];
        int l = child_lvl(r.left, lvl);
        if (r.kind == grammar::RUNLEN) {
            lvl[i] = 1 + l;                       // lvl(N->B^t) = 1 + lvl(B)
        } else {
            int rt = child_lvl(r.right, lvl);
            lvl[i] = 1 + std::max(l, rt);         // lvl(N->AB) = 1 + max
        }
        if (lvl[i] > max_depth) max_depth = lvl[i];
    }
    out.max_depth = max_depth;

    // --- (4) Highest rule level at which each symbol is referenced. ---
    // (C references do NOT force a symbol to the top: it is emitted at its own
    //  level instead, so only rule references govern pass-through length.)
    std::vector<int> max_ref_lvl(nt_num, 0);
    for (std::size_t i = 0; i < nt_num; ++i) max_ref_lvl[i] = lvl[i];
    for (std::size_t i = 0; i < nt_num; ++i) {
        const Rule& r = g.rules[i];
        if (r.left >= alpha)
            max_ref_lvl[r.left - alpha] = std::max(max_ref_lvl[r.left - alpha], lvl[i]);
        if (r.kind == grammar::BINARY && r.right >= alpha)
            max_ref_lvl[r.right - alpha] = std::max(max_ref_lvl[r.right - alpha], lvl[i]);
    }

    // --- (5) Build layered levels with pass-through completion. ---
    const int stride = max_depth + 1;
    std::vector<int> rep(nt_num * static_cast<std::size_t>(stride), -1);
    out.layered_levels.assign(max_depth + 1, {});
    std::vector<std::vector<int>> level_orig_symbol(max_depth + 1);

    auto ref_in_prev = [&](int sym, int k) -> int {
        // Reference to a child as seen in frontier k: terminal -> -sym-1 (T),
        // non-terminal -> its slot index in level k.
        return (sym < alpha) ? (-sym - 1)
                             : rep[(sym - alpha) * static_cast<std::size_t>(stride) + k];
    };

    for (int k = 1; k <= max_depth; ++k) {
        // Original non-terminals whose level is exactly k.
        for (std::size_t i = 0; i < nt_num; ++i) {
            if (lvl[i] != k) continue;
            const Rule& r = g.rules[i];
            LayeredRule lr;
            lr.left = ref_in_prev(r.left, k - 1);
            if (r.kind == grammar::RUNLEN) {
                lr.right = -1;                       // zero terminal T[0] = 0
                lr.coeff = static_cast<float>(r.rep);// eval(N) = t * eval(B)
            } else {
                lr.right = ref_in_prev(r.right, k - 1);
                lr.coeff = 1.0f;
            }
            int idx = static_cast<int>(out.layered_levels[k].size());
            out.layered_levels[k].push_back(lr);
            rep[i * static_cast<std::size_t>(stride) + k] = idx;
            level_orig_symbol[k].push_back(static_cast<int>(i));
        }

        // Pass-throughs carrying level-(k-1) values still referenced above.
        if (k > 1) {
            int num_prev = static_cast<int>(out.layered_levels[k - 1].size());
            for (int j = 0; j < num_prev; ++j) {
                int orig = level_orig_symbol[k - 1][j];
                if (max_ref_lvl[orig] > k - 1) {
                    LayeredRule lr;
                    lr.left = j;        // slot in previous frontier
                    lr.right = -1;      // zero terminal
                    lr.coeff = 1.0f;
                    int idx = static_cast<int>(out.layered_levels[k].size());
                    out.layered_levels[k].push_back(lr);
                    rep[orig * static_cast<std::size_t>(stride) + k] = idx;
                    level_orig_symbol[k].push_back(orig);
                }
            }
        }
    }

    // --- (6) Emit-on-the-spot lists from the top sequence C. ---
    out.emit_by_level.assign(max_depth + 1, {});
    {
        int cur_row = 0;
        for (std::size_t j = 0; j < g.C.size(); ++j) {
            int c = g.C[j];
            if (c == 0) {                 // row delimiter $
                cur_row++;
            } else if (c < alpha) {       // terminal occurrence
                out.term_emit.push_back(std::make_pair(c, cur_row));
            } else {                      // non-terminal: emit at its own level
                int N = c - alpha;
                int k = lvl[N];
                int pos = rep[N * static_cast<std::size_t>(stride) + k];
                out.emit_by_level[k].push_back(std::make_pair(pos, cur_row));
            }
        }
    }

    // --- (7) Flatten rules + emissions; collect stats. ---
    out.level_offsets.push_back(0);
    out.level_widths.push_back(0);
    int max_width = 0;
    for (int k = 1; k <= max_depth; ++k) {
        out.flattened_rules.insert(out.flattened_rules.end(),
                                   out.layered_levels[k].begin(),
                                   out.layered_levels[k].end());
        out.level_offsets.push_back(static_cast<int>(out.flattened_rules.size()));
        int w = static_cast<int>(out.layered_levels[k].size());
        out.level_widths.push_back(w);
        if (w > max_width) max_width = w;
    }
    out.max_width = max_width;
    out.num_layered = out.flattened_rules.size();
    out.num_passthrough = out.num_layered - out.num_raw_nt;

    out.emit_offsets.assign(max_depth + 1, 0);
    for (int k = 1; k <= max_depth; ++k) {
        for (const auto& pr : out.emit_by_level[k]) {
            out.emit_pos.push_back(pr.first);
            out.emit_row.push_back(pr.second);
        }
        out.emit_offsets[k] = static_cast<int>(out.emit_pos.size());
    }
    for (const auto& pr : out.term_emit) {
        out.term_emit_sym.push_back(pr.first);
        out.term_emit_row.push_back(pr.second);
    }

    // --- (8) Terminal compaction: keep only terminals actually referenced. ---
    // The terminal id space is the dense (value,column) grid of size alpha, but
    // only the distinct pairs that occur are ever read (as negative children of
    // rules or as term_emit symbols). Remap the referenced ids to a contiguous
    // [0,U) range so the engine's init_T touches U terminals per vector instead
    // of the full grid -- a large win when values are near-unique and the grid is
    // sparse (continuous-valued matrices). Terminal 0 (the zero terminal, T[0]=0)
    // always keeps compact id 0.
    {
        std::vector<int> remap(alpha, -1);
        remap[0] = 0;                                  // zero terminal stays id 0
        auto mark = [&](int child) {                   // flag a terminal as used
            if (child < 0) { int t = -child - 1; if (remap[t] < 0) remap[t] = 0; }
        };
        for (const auto& lr : out.flattened_rules) { mark(lr.left); mark(lr.right); }
        for (int sym : out.term_emit_sym) { if (remap[sym] < 0) remap[sym] = 0; }

        out.term_orig.clear();
        out.term_orig.push_back(0);                    // compact id 0 = zero terminal
        for (int t = 1; t < alpha; ++t) {
            if (remap[t] == 0) {                       // referenced -> assign next id
                remap[t] = static_cast<int>(out.term_orig.size());
                out.term_orig.push_back(t);
            }
        }
        int new_alpha = static_cast<int>(out.term_orig.size());

        auto remap_child = [&](int& child) {
            if (child < 0) { int t = -child - 1; child = -remap[t] - 1; }
        };
        for (auto& lr : out.flattened_rules) { remap_child(lr.left); remap_child(lr.right); }
        for (auto& lvlvec : out.layered_levels)
            for (auto& lr : lvlvec) { remap_child(lr.left); remap_child(lr.right); }
        for (auto& s : out.term_emit_sym) s = remap[s];
        for (auto& pr : out.term_emit) pr.first = remap[pr.first];

        out.alpha = new_alpha;                         // compacted terminal count
    }

    return out;
}