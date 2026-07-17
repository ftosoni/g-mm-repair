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

"""Convert a 1000G phased VCF into a genotype matrix for mm-repair.
Rows = individuals, columns = biallelic SNPs in genome order (so adjacent
columns are in linkage disequilibrium, and each column has values in {0,1,2}).
Output: dense int32, row-major (individual-major) -> feed matrepair with --i32.
This is the Chr20/21/22 slicer used to build the real genotype matrices of the
manuscript (REPRODUCIBILITY.md sec. A); [max_variants] caps the column count for
the 10^5-variant subsets vs. the full-width chromosomes.
Usage: vcf2mat.py <chr.vcf.gz> <out_matrix> [max_variants]
"""
import sys, gzip, numpy as np

vcf, out = sys.argv[1], sys.argv[2]
maxv = int(sys.argv[3]) if len(sys.argv) > 3 else 100000

cols = []
kept = 0
with gzip.open(vcf, 'rt') as f:
    for line in f:
        if line.startswith('#'):
            continue
        fields = line.rstrip('\n').split('\t')
        ref, alt = fields[3], fields[4]
        if len(ref) != 1 or len(alt) != 1:   # biallelic SNPs only
            continue
        gts = fields[9:]
        arr = np.fromiter(
            (((1 if g[0] == '1' else 0) + (1 if g[2] == '1' else 0)) for g in gts),
            dtype=np.int8, count=len(gts))
        cols.append(arr)
        kept += 1
        if kept % 20000 == 0:
            print(f"  ...{kept} variants", flush=True)
        if kept >= maxv:
            break

M = np.stack(cols, axis=1)        # (individuals, variants)
nz = np.count_nonzero(M)
print(f"matrix {M.shape}, nonzero {100.0*nz/M.size:.1f}% ({nz} nnz)")
M.astype(np.int32).tofile(out)
print(f"rows(individuals)={M.shape[0]} cols(variants)={M.shape[1]}")
