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

"""SuiteSparse:GraphBLAS (python-graphblas) baseline for the graph relations.

Reads a textual sparse edge list "row col" (as produced by `process_wikidata.py sparse`),
builds the boolean adjacency matrix, and times one *semiring-native* mat-vec (GrB_mxv) —
the direct CPU analog of the engine's single-vector graph sweep (tab:graph):
  - Boolean : lor_land semiring  (one BFS-frontier / reachability step)
  - Tropical: min_plus semiring  (one Bellman-Ford relaxation step)

Reports both single-threaded and 20-threaded times (to match the manuscript's CPU seq /
OpenMP-20 references). SuiteSparse:GraphBLAS is the semiring-native baseline named in the
methodology (there is no vendor dense kernel for these semirings).

Usage: graphblas_bench.py <sparse_file> <rows> <cols> <bool|tropical> [iters]
"""
import sys, os, time
import numpy as np
import graphblas as gb
from graphblas import Matrix, Vector, semiring, dtypes


def main():
    if len(sys.argv) < 5:
        print("Usage: graphblas_bench.py <sparse_file> <rows> <cols> <bool|tropical> [iters]")
        sys.exit(1)
    path = sys.argv[1]
    rows, cols = int(sys.argv[2]), int(sys.argv[3])
    sr = sys.argv[4]
    iters = int(sys.argv[5]) if len(sys.argv) > 5 else 50

    # fast load of the "row col" text edge list
    ij = np.array(open(path).read().split(), dtype=np.int64).reshape(-1, 2)
    r, c = ij[:, 0], ij[:, 1]
    nnz = len(r)
    rows = max(rows, int(r.max()) + 1)   # be robust to the exact remapped extent
    cols = max(cols, int(c.max()) + 1)

    # CROSSCHECK mode: load the engine's shared x and y (dumped by gpu_test), run the
    # semiring-native mat-vec on the SAME vector, and assert the result matches the
    # engine's -- GraphBLAS is thus certified against the engine on one input vector.
    cc = os.environ.get("CROSSCHECK")
    if cc:
        sr_full = {"bool": "boolean", "tropical": "tropical"}.get(sr, sr)
        xf = np.fromfile(f"{cc}.{sr_full}.x", dtype=np.float32)
        ey = np.fromfile(f"{cc}.{sr_full}.y", dtype=np.float64)
        ncols, nrows = len(xf), len(ey)
        if sr == "bool":
            A = Matrix.from_coo(r, c, True, nrows=nrows, ncols=ncols, dtype=dtypes.BOOL)
            xv = Vector.from_dense((xf != 0).astype(bool))
            yd = semiring.lor_land(A @ xv).new().to_dense(fill_value=False).astype(np.int64)
            ei = (ey != 0).astype(np.int64)
            op = "lor_land"
        else:
            A = Matrix.from_coo(r, c, 1, nrows=nrows, ncols=ncols, dtype=dtypes.INT64)
            xv = Vector.from_dense(xf.astype(np.int64))
            yd = np.minimum(semiring.min_plus(A @ xv).new().to_dense(fill_value=10**9).astype(np.int64), 10**9)
            ei = np.minimum(ey.astype(np.int64), 10**9)
            op = "min_plus"
        mism = int((yd != ei).sum())
        print(f"CROSSCHECK graphblas {op} vs engine: mismatches={mism}  "
              f"{'SUCCESS' if mism == 0 else 'FAILURE'} (rows={nrows} cols={ncols})")
        return

    if sr == "bool":
        A = Matrix.from_coo(r, c, True, nrows=rows, ncols=cols, dtype=dtypes.BOOL)
        x = Vector.from_dense(np.ones(cols, dtype=bool))
        S = semiring.lor_land
    elif sr == "tropical":
        A = Matrix.from_coo(r, c, 1, nrows=rows, ncols=cols, dtype=dtypes.INT64)
        x = Vector.from_dense(np.zeros(cols, dtype=np.int64))
        S = semiring.min_plus
    else:
        print(f"bad semiring: {sr}"); sys.exit(1)

    for nt in (1, 20):
        gb.ss.config["nthreads"] = nt
        S(A @ x).new()                       # warmup (also triggers JIT/analysis)
        t = time.perf_counter()
        for _ in range(iters):
            S(A @ x).new()
        dt = (time.perf_counter() - t) / iters * 1000.0
        print(f"GraphBLAS {sr} nthreads={nt:>2}: {dt:.4f} ms/vector  "
              f"(rows={rows} cols={cols} nnz={nnz})")


if __name__ == "__main__":
    main()
