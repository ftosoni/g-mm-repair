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

#ifndef GRAMMAR_H
#define GRAMMAR_H

#include <vector>
#include <cstdint>

// A grammar over the CSRV string (the (C, R, V) representation of the
// manuscript, sec:background), built from RePair (compression-first): rules are
// read from mm-repair's .vc.R, all of kind BINARY. This is the construction the
// manuscript evaluates (section "A proper-layered streaming engine", "Step 1: RePair (unchanged)"), and it is fed
// into build_schedule().
//
// Symbol encoding (identical to mm-repair's .vc / .vc.R convention):
//   a child/symbol value  s  is a *terminal*    if  0 <= s < alpha,
//                            a *non-terminal* rule index (s - alpha) otherwise.
//   Terminal 0 is the row delimiter "$"; terminal p in [1,alpha) decodes to
//   value index (p-1)/cols and column (p-1)%cols (see the GPU engine).
//
// A rule with index i is the non-terminal whose symbol value is (alpha + i).
// Numbering is topological: every child of rule i has a strictly smaller level
// (RePair's index order i<j-whenever-N_i-occurs-in-N_j guarantees this). This is
// what lets build_schedule() assign levels in one pass.
namespace grammar {

enum RuleKind : int {
    BINARY = 0,   // N -> A B          (left = A, right = B)
    RUNLEN = 1    // N -> B^t , t >= 2 (left = B, right unused, rep = t).
                  // Reserved run-length kind the layered engine can evaluate; the
                  // RePair path used in the paper produces BINARY rules only.
};

struct Rule {
    int kind;     // RuleKind
    int left;     // first child (terminal or alpha+idx)
    int right;    // BINARY: second child; RUNLEN: ignored
    int rep;      // RUNLEN: repetition count t (>=2); BINARY: 1
};

struct Grammar {
    int alpha = 0;            // terminal/non-terminal boundary
    std::vector<Rule> rules;  // rules[i] defines non-terminal (alpha + i)
    std::vector<int>  C;      // top sequence: symbols with 0 marking row ends
};

} // namespace grammar

#endif // GRAMMAR_H