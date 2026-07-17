#!/usr/bin/env python3

# This file is part of g-mm-repair
# <https://github.com/ftosoni/g-mm-repair>.
# Copyright (c) 2026 Francesco Tosoni.
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#     http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

"""Compare two binary float64 vectors and print a single CROSSCHECK PASS/FAIL line.

Used by the reproduce.sh `crosscheck` driver to compare a baseline's output vector
(e.g. mm-repair's re32mm) against the engine's dumped y on the same shared input.

Usage: crosscheck_compare.py <label> <a.f64> <b_engine.f64> [tol]
"""
import sys
import numpy as np


def main():
    label, fa, fb = sys.argv[1], sys.argv[2], sys.argv[3]
    tol = float(sys.argv[4]) if len(sys.argv) > 4 else 1e-4
    a = np.fromfile(fa, dtype=np.float64)
    b = np.fromfile(fb, dtype=np.float64)
    if len(a) != len(b):
        print(f"CROSSCHECK {label} vs engine: LENGTH MISMATCH n1={len(a)} n2={len(b)}  FAILURE")
        return
    rel = float((np.abs(a - b) / (1.0 + np.abs(b))).max())
    print(f"CROSSCHECK {label} vs engine: max_rel_diff={rel:.3e}  "
          f"{'SUCCESS' if rel < tol else 'FAILURE'} (n={len(a)})")


if __name__ == "__main__":
    main()
