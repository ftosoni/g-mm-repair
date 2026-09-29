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

# Build the knowledge-graph relation matrices used in the manuscript's semiring
# experiments (sec:semiring; REPRODUCIBILITY.md sec. D). From a triples file <s p o> plus its predicate
# dictionary <dat.P>, four modes:
#   analyze -> list predicate IDs + counts (pick the relations to build)
#   build   -> DENSE int32 matrix for ONE relation  (small relations only, e.g. YAGO)
#   sparse  -> textual "row col" edge list per relation, feeding `matrepair --bool`
#              WITHOUT ever materializing the dense matrix (needed at Wikidata scale)
#   square  -> like sparse but remaps subjects+objects into ONE shared node space,
#              for endomorphic relations where reachability/SSSP fixpoint y=M(x)x applies
# All outputs use the same dense-int32 / edge-list formats the genotype matrices use,
# so the CSRV/RePair pipeline and the GPU engine are entirely unchanged.
#
# Graph datasets (YAGO2s, Wikidata) from Zenodo: https://zenodo.org/record/7254968
# Reference: Diego Arroyuelo, Adrián Gómez-Brandón, Aidan Hogan, Gonzalo Navarro, & Javiel Rojas-Ledesma. (2022). Datasets of Time- and Space-Efficient Regular Path Queries.
import sys
from collections import defaultdict
import numpy as np

def analyze_yago(dat_file, p_file):
    # Read predicate names
    pred_names = {}
    with open(p_file, 'r', encoding='utf-8', errors='ignore') as f:
        for line in f:
            parts = line.strip().split(None, 1)
            if len(parts) == 2:
                pred_names[int(parts[0])] = parts[1]

    # Count triples per predicate
    pred_counts = defaultdict(int)
    pred_subjects = defaultdict(set)
    pred_objects = defaultdict(set)

    print("Reading triples file...")
    with open(dat_file, 'r') as f:
        for line in f:
            parts = line.strip().split()
            if len(parts) == 3:
                s, p, o = int(parts[0]), int(parts[1]), int(parts[2])
                pred_counts[p] += 1
                pred_subjects[p].add(s)
                pred_objects[p].add(o)

    print("\nPredicate Statistics:")
    print(f"{'ID':<5} | {'Predicate Name':<50} | {'Triples':<12} | {'Unique S':<12} | {'Unique O':<12} | {'Dense Size':<15}")
    print("-" * 115)
    for p in sorted(pred_counts.keys()):
        name = pred_names.get(p, f"Unknown_{p}")
        num_triples = pred_counts[p]
        num_s = len(pred_subjects[p])
        num_o = len(pred_objects[p])
        dense_size_bytes = num_s * num_o * 4
        dense_size_mb = dense_size_bytes / (1024 * 1024)
        print(f"{p:<5} | {name:<50} | {num_triples:<12} | {num_s:<12} | {num_o:<12} | {dense_size_mb:.2f} MB")

def build_predicate_matrix(dat_file, target_p, out_file):
    """DENSE builder (int32 rows x cols). Only feasible for small relations
    (rows*cols must fit in RAM/disk). For large graphs use build_sparse()."""
    print(f"Building matrix for predicate {target_p}...")
    triples = []
    subjects = set()
    objects = set()

    with open(dat_file, 'r') as f:
        for line in f:
            parts = line.strip().split()
            if len(parts) == 3:
                s, p, o = int(parts[0]), int(parts[1]), int(parts[2])
                if p == target_p:
                    triples.append((s, o))
                    subjects.add(s)
                    objects.add(o)

    sorted_s = sorted(list(subjects))
    sorted_o = sorted(list(objects))

    s_map = {s: i for i, s in enumerate(sorted_s)}
    o_map = {o: j for j, o in enumerate(sorted_o)}

    rows = len(sorted_s)
    cols = len(sorted_o)

    print(f"Creating dense matrix of shape {rows}x{cols} ({rows * cols} cells)...")
    M = np.zeros((rows, cols), dtype=np.int8)
    for s, o in triples:
        M[s_map[s], o_map[o]] = 1

    nz = np.count_nonzero(M)
    print(f"Matrix nonzeros: {nz} ({100.0 * nz / M.size:.4f}%)")
    M.astype(np.int32).tofile(out_file)
    print(f"Saved binary matrix to {out_file}. rows={rows} cols={cols}")

def build_sparse(dat_file, specs):
    """SPARSE builder: one pass over the triples file, emits a textual
    'row col' list (0-indexed, sorted by row then col) per predicate, ready for
    `matrepair --bool <name>.sparse <rows> <cols>`. Never materializes a dense
    matrix, so it scales to large Wikidata relations. `specs` is a list of
    (predicate_id, out_name)."""
    targets = {pid: name for pid, name in specs}
    edges = {pid: [] for pid in targets}
    subs = {pid: set() for pid in targets}
    objs = {pid: set() for pid in targets}

    print(f"Scanning {dat_file} for predicates {sorted(targets)} ...")
    with open(dat_file, 'r') as f:
        for line in f:
            parts = line.split()
            if len(parts) == 3 and parts[1] in targets:
                pid = parts[1]
                edges[pid].append((int(parts[0]), int(parts[2])))
                subs[pid].add(parts[0])
                objs[pid].add(parts[2])

    print(f"\n{'name':<22}{'pid':>6}{'rows':>12}{'cols':>12}{'nnz':>12}")
    for pid, name in specs:
        s_map = {s: i for i, s in enumerate(sorted(subs[pid], key=int))}
        o_map = {o: j for j, o in enumerate(sorted(objs[pid], key=int))}
        rows, cols = len(s_map), len(o_map)
        pairs = sorted((s_map[str(s)], o_map[str(o)]) for s, o in edges[pid])
        out = f"{name}.sparse"
        with open(out, 'w') as fo:
            fo.writelines(f"{r} {c}\n" for r, c in pairs)
        print(f"{name:<22}{pid:>6}{rows:>12}{cols:>12}{len(pairs):>12}")

def build_square(dat_file, pid, name):
    """SQUARE builder for an *endomorphic* relation (source and target the same
    entity type, e.g. subclassOf, cites). Remaps subjects and objects into ONE
    shared 0-indexed node space (their union), so the adjacency matrix is N x N
    and reachability/SSSP fixpoint iteration y = M (x) x is well defined. Emits
    a '<name>.sparse' edge list and prints N and nnz."""
    edges = []
    nodes = set()
    print(f"Scanning {dat_file} for predicate {pid} (square/union remap) ...")
    with open(dat_file, 'r') as f:
        for line in f:
            p = line.split()
            if len(p) == 3 and p[1] == pid:
                edges.append((p[0], p[2]))
                nodes.add(p[0]); nodes.add(p[2])
    nmap = {v: i for i, v in enumerate(sorted(nodes, key=int))}
    N = len(nmap)
    pairs = sorted((nmap[s], nmap[o]) for s, o in edges)
    with open(f"{name}.sparse", 'w') as fo:
        fo.writelines(f"{r} {c}\n" for r, c in pairs)
    print(f"{name} pid={pid} N={N} nnz={len(pairs)}  (square {N}x{N})")


def main():
    if len(sys.argv) < 3:
        print("Usage:")
        print("  Analyze:      process_wikidata.py <dat> <dat.P> analyze")
        print("  Build dense:  process_wikidata.py <dat> <dat.P> build <predicate_id> <out_matrix>")
        print("  Build sparse: process_wikidata.py <dat> <dat.P> sparse <id1>:<name1> [<id2>:<name2> ...]")
        print("  Build square: process_wikidata.py <dat> <dat.P> square <predicate_id> <out_name>")
        sys.exit(1)

    dat_file = sys.argv[1]
    p_file = sys.argv[2]
    mode = sys.argv[3]

    if mode == "analyze":
        analyze_yago(dat_file, p_file)
    elif mode == "build":
        target_p = int(sys.argv[4])
        out_file = sys.argv[5]
        build_predicate_matrix(dat_file, target_p, out_file)
    elif mode == "sparse":
        specs = [(a.split(":", 1)[0], a.split(":", 1)[1]) for a in sys.argv[4:]]
        build_sparse(dat_file, specs)
    elif mode == "square":
        build_square(dat_file, sys.argv[4], sys.argv[5])

if __name__ == "__main__":
    main()
