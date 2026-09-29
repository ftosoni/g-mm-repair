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

"""cugraph (RAPIDS) GPU baseline for the graph relations.

Reads a textual sparse edge list "row col" (from `process_wikidata.py sparse`) and times the
canonical GPU graph traversals as end-to-end baselines for the semiring sweeps:
  - bool     -> cugraph.bfs   (full reachability from a source; Boolean analog)
  - tropical -> cugraph.sssp  (single-source shortest paths; Tropical analog)

NOTE: unlike the engine's *single* semiring mat-vec (one BFS frontier / one relaxation), cuGraph
BFS/SSSP run the *whole* traversal to convergence — so this is an end-to-end reference, not a
per-mat-vec one. Objects are offset by `rows` so subjects and objects are disjoint node ids
(a proper bipartite graph).

Usage: cugraph_bench.py <sparse_file> <rows> <cols> <bool|tropical> [iters]
"""
import sys, time
import numpy as np
import cudf, cugraph


def main():
    if len(sys.argv) < 5:
        print("Usage: cugraph_bench.py <sparse_file> <rows> <cols> <bool|tropical> [iters]")
        sys.exit(1)
    path = sys.argv[1]
    rows = int(sys.argv[2])
    sr = sys.argv[4]
    iters = int(sys.argv[5]) if len(sys.argv) > 5 else 50

    ij = np.array(open(path).read().split(), dtype=np.int64).reshape(-1, 2)
    src = ij[:, 0]
    dst = ij[:, 1] + rows                      # disjoint object ids -> bipartite graph
    nnz = len(src)

    df = cudf.DataFrame({"src": src, "dst": dst})
    G = cugraph.Graph(directed=True)
    if sr == "bool":
        G.from_cudf_edgelist(df, source="src", destination="dst")
        op = lambda: cugraph.bfs(G, start=0)
        label = "BFS (full reachability)"
    elif sr == "tropical":
        df["weight"] = np.ones(nnz, dtype=np.float32)
        G.from_cudf_edgelist(df, source="src", destination="dst", edge_attr="weight")
        op = lambda: cugraph.sssp(G, source=0)
        label = "SSSP (single-source shortest paths)"
    else:
        print(f"bad semiring: {sr}"); sys.exit(1)

    op()                                       # warmup
    t = time.perf_counter()
    for _ in range(iters):
        op()
    dt = (time.perf_counter() - t) / iters * 1000.0
    print(f"cuGraph {sr} [{label}]: {dt:.4f} ms/traversal  (nodes={rows}+cols, nnz={nnz})")


if __name__ == "__main__":
    main()
