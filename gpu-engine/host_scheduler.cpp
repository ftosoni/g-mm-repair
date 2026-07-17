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
#include <algorithm>
#include <string>
#include <filesystem>

// Standalone, CUDA-free inspector for a RePair grammar (.vc.R): reads the raw
// rules, computes each rule's level by the manuscript's level equation (before any pass-through
// completion), and prints the depth L and the per-level width distribution.
// Handy for a quick look at a grammar's shape; the numbers the paper reports
// (L, w* AFTER completion, +pt) come from gpu_test / build_schedule, not here.
int main(int argc, char** argv) {
    if (argc < 2) {
        std::cerr << "Usage: " << argv[0] << " <matrix_base_name>" << std::endl;
        std::cerr << "Example: " << argv[0] << " mm-repair/data/geno22" << std::endl;
        return 1;
    }

    std::string base_name = argv[1];
    std::string rules_filename = base_name + ".vc.R";

    std::ifstream file(rules_filename, std::ios::binary | std::ios::ate);
    if (!file.is_open()) {
        std::cerr << "Error: Could not open rules file: " << rules_filename << std::endl;
        return 1;
    }

    std::streamsize file_size = file.tellg();
    file.seekg(0, std::ios::beg);

    if (file_size < static_cast<std::streamsize>(sizeof(int))) {
        std::cerr << "Error: Rules file is too small." << std::endl;
        return 1;
    }

    int alpha = 0;
    if (!file.read(reinterpret_cast<char*>(&alpha), sizeof(int))) {
        std::cerr << "Error: Could not read alphabet size (Alpha)." << std::endl;
        return 1;
    }

    size_t rules_bytes = file_size - sizeof(int);
    size_t nt_num = rules_bytes / (2 * sizeof(int));

    std::cout << "Alphabet size (Alpha): " << alpha << std::endl;
    std::cout << "Number of non-terminals (NTnum): " << nt_num << std::endl;

    std::vector<int> rules(nt_num * 2);
    if (!file.read(reinterpret_cast<char*>(rules.data()), rules_bytes)) {
        std::cerr << "Error: Could not read rules array." << std::endl;
        return 1;
    }

    std::vector<int> lvl(nt_num, 0);
    int max_depth = 0;

    for (size_t i = 0; i < nt_num; ++i) {
        int left = rules[2 * i];
        int right = rules[2 * i + 1];

        int left_lvl = (left >= alpha) ? lvl[left - alpha] : 0;
        int right_lvl = (right >= alpha) ? lvl[right - alpha] : 0;

        lvl[i] = 1 + std::max(left_lvl, right_lvl);
        if (lvl[i] > max_depth) {
            max_depth = lvl[i];
        }
    }

    std::cout << "Grammar depth (L): " << max_depth << std::endl;

    std::vector<size_t> level_widths(max_depth + 1, 0);
    for (size_t i = 0; i < nt_num; ++i) {
        level_widths[lvl[i]]++;
    }

    std::cout << "\nLevel Width Distribution:" << std::endl;
    std::cout << "Level\tWidth" << std::endl;
    for (int k = 1; k <= max_depth; ++k) {
        std::cout << k << "\t" << level_widths[k] << std::endl;
    }

    return 0;
}