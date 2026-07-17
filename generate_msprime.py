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

"""
Coalescent genotype-matrix generator (citable, reproducible alternative to UK Biobank).

Simulates diploid genotypes under the coalescent-with-recombination using msprime
[Kelleher, Etheridge & McVean 2016; Baumdicker et al. 2022] and writes them in the
same dense int32 {0,1,2} row-major format (individuals x variants) consumed by
mm-repair / vcf2mat.py. This is the sole synthetic generator used in the manuscript
(a citable, seeded coalescent simulator).

Linkage disequilibrium (the repetitive structure RePair exploits) is governed by the
recombination-to-mutation ratio: a LOW recombination rate yields long shared haplotype
blocks (high LD, highly compressible); a HIGH recombination rate breaks them up
(low LD, poorly compressible). This replaces the ad-hoc founder-block model with a
standard population-genetic simulator that can be cited and exactly reproduced via --seed.

Usage:
  generate_msprime.py <rows> <cols> <out_matrix> [recomb_rate] [seed] [Ne] [mu]

  rows         number of diploid individuals (matrix rows)
  cols         number of variant sites (matrix columns)
  out_matrix   output path (raw int32, row-major, rows x cols)
  recomb_rate  per-base recombination rate (default 1e-8 ~ human).
               LD presets used in the paper:  high LD = 1e-9,  low LD = 1e-7.
  seed         RNG seed for full reproducibility (default 42)
  Ne           effective population size (default 10000)
  mu           per-base mutation rate (default 1e-8)
"""
import sys
import numpy as np

try:
    import msprime
except ImportError:
    sys.exit("ERROR: msprime is required. Install with:  pip install msprime")


def simulate_genotypes(rows, cols, recomb_rate, seed, Ne, mu):
    n_hap = 2 * rows  # haplotypes (diploid individuals)
    # Expected segregating sites S ~ theta * H_{n-1}, theta = 4*Ne*mu*L.
    # Invert for an initial sequence length giving ~2x the requested columns of margin.
    H = float(np.sum(1.0 / np.arange(1, n_hap)))
    seq_len = max(cols / (4.0 * Ne * mu * H) * 2.0, 1e4)

    for attempt in range(6):  # grow the genome until enough polymorphic sites appear
        ts = msprime.sim_ancestry(
            samples=rows, ploidy=2, population_size=Ne,
            sequence_length=int(seq_len), recombination_rate=recomb_rate,
            random_seed=seed,
        )
        mts = msprime.sim_mutations(ts, rate=mu, random_seed=seed)
        # genotype_matrix(): shape (num_sites, n_hap). Binarize to ancestral(0)/derived(1)
        # so the result is biallelic regardless of the mutation model.
        G = (mts.genotype_matrix() != 0).astype(np.int8)
        # Drop monomorphic sites (real SNP panels are polymorphic).
        seg = (G.sum(axis=1) > 0) & (G.sum(axis=1) < n_hap)
        G = G[seg]
        if G.shape[0] >= cols:
            break
        seq_len *= 2.0
        print(f"  [retry {attempt+1}] only {G.shape[0]} sites < {cols}; "
              f"doubling sequence length to {seq_len:.3g} bp")
    else:
        sys.exit(f"ERROR: could not reach {cols} sites; increase mu or Ne.")

    G = G[:cols]                       # exactly `cols` variants
    diploid = G[:, 0::2] + G[:, 1::2]  # (cols, rows) with values in {0,1,2}
    return diploid.T.astype(np.int32)  # (rows, cols)


def main():
    if len(sys.argv) < 4:
        print(__doc__)
        sys.exit(1)
    rows = int(sys.argv[1])
    cols = int(sys.argv[2])
    out_file = sys.argv[3]
    recomb_rate = float(sys.argv[4]) if len(sys.argv) > 4 else 1e-8
    seed = int(sys.argv[5]) if len(sys.argv) > 5 else 42
    Ne = int(sys.argv[6]) if len(sys.argv) > 6 else 10000
    mu = float(sys.argv[7]) if len(sys.argv) > 7 else 1e-8

    print(f"Simulating {rows} individuals x {cols} variants via msprime "
          f"(recomb={recomb_rate:g}, mu={mu:g}, Ne={Ne}, seed={seed})...")
    M = simulate_genotypes(rows, cols, recomb_rate, seed, Ne, mu)
    nz = int(np.count_nonzero(M))
    print(f"matrix {M.shape}, nonzero {100.0*nz/M.size:.2f}% ({nz} nnz)")
    M.tofile(out_file)
    print(f"Saved binary matrix to {out_file}. "
          f"rows(individuals)={M.shape[0]} cols(variants)={M.shape[1]}")


if __name__ == "__main__":
    main()
