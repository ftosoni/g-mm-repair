# Reproducibility Guide for "Streaming Right Multiplication over Grammar-Compressed Matrices"

> **This is an optional deep-dive.** To reproduce the paper you only need the
> [README](README.md): download the Zenodo package and run `./reproduce.sh` — that script,
> together with its canonical logs under `manuscript/logs/`, is the single source of truth.
> This guide is the *reference* behind it: the exact per-table commands `reproduce.sh`
> automates (§3), the log→artifact mapping (§4), the cross-implementation correctness protocol
> (§2.F), and — for the curious — how every dataset in the package was derived from scratch
> (§2.A–§2.E). None of §2.A–§2.E is required when you use the published Zenodo package.

This document gives step-by-step instructions to reproduce **all and only** the tables and figures presented in the manuscript. The manuscript reports two experiment families: **genotype matrices** (1000 Genomes + synthetic haplotypes, `PlusTimes` semiring) and **knowledge-graph relation matrices** (Wikidata relations, `Boolean`/`Tropical` semirings).

These experiments were originally profiled and measured on a remote prototyping node (**NVIDIA Grace-Blackwell GB10 node**, featuring unified coherent CPU-GPU memory of 119 GiB, CUDA 13.0, g++ 13.3, and Ubuntu 24.04 LTS). Per the manuscript (§ Limitations), time/energy figures are board-dependent; the structural figures ($|\mathcal{R}|$, $L$, $w^{*}$, $+\text{pt}$) are architecture-independent and reproducible on any host.

> ✅ **Status.** Every measured table and figure in this guide is backed by a canonical log under `manuscript/logs/`. The genotype/crossover family (§A–§C, Tables 4.1–4.3, A.1, B.1, Figs 4.3, A.1) was regenerated end-to-end with the `msprime` simulator on node `spark-a459` (2026-07-14); the Wikidata/SWH graph logs (§D–§E, Tables 5.1, 5.2, B.2) are unchanged.

The complete inventory of manuscript artifacts this guide reproduces:

| Artifact | Content | Datasets |
|---|---|---|
| Table 4.1 (`tab:geno_through`) | Structural figures, genotypes | 6 real + 5 synthetic genotype matrices |
| Table 4.2 (`tab:geno_time`) | Avg time/vector, genotypes | same 11 genotype matrices |
| Table 4.3 (`tab:geno`) | Space & energy vs. cuSPARSE | same 11 + `crossover_synth` |
| Table 5.1 (`tab:graph_scale`) | Scale/space at 10M–1.2G edges | 2 largest Wikidata relations + SWH software graph |
| Table 5.2 (`tab:graph`) | Boolean & Tropical graph product | 5 Wikidata relation matrices |
| Table A.1 (`tab:geno_spmm`) | Batched right product (SpMM) | 6 real + 5 synthetic genotype matrices |
| Table B.1 (`tab:build`) | Host-side construction cost + amortization | 11 genotype matrices |
| Table B.2 (`tab:graph_struct`) | Structural figures, Wikidata | 5 Wikidata relation matrices |
| Figure 4.3 (`fig:space`) | Memory footprint (incl. `crossover_synth` at billion-nnz scale) | genotype matrices + crossover |
| Figure A.1 (`fig:batched`) | Batched throughput vs. $B$ | `geno22full` |

(The remaining figures — `fig:matrix`, `fig:mmr_rs`, `fig:dag`, `fig:completion`, `fig:sweep`, `fig:trace` — are schematic TikZ/`includegraphics` illustrations of the running example, not measured results, and require no experiment to regenerate.)

> **Reproduction flow.** In order: **§1** compile the binaries → **§2** build/download the datasets (skip if already on the node) → **`./reproduce.sh all`**, which runs every core experiment into `manuscript/logs/`, then regenerates the tables (`manuscript/tables/*.tex`) and figure data (`manuscript/figures/data/*.dat`) via `extract_results.py`. `reproduce.sh` and its canonical logs are the source of truth; the per-table commands in §3 are the manual equivalents of what it automates, and §4 documents the log→artifact mapping. Three heavy families (`grammar`, `graphscale`, `crosscheck`) are **not** in `all` and are run on demand (§3/§4).

---

## 1. Compilation

Before running any experiments, compile the GPU engine, baseline tests, and CPU-side `mm-repair` binaries:

### Build the GPU Engine & cuSPARSE Baselines
```bash
cd gpu-engine
make clean
make                  # Compiles 'gpu_test' and 'host_scheduler'
make cusparse_test    # Compiles the cuSPARSE CSR SpMV/SpMM baseline
cd ..
```
> **GPU architecture.** `gpu-engine/Makefile` sets `-arch=sm_121` for the Grace-Blackwell GB10 node. On any other GPU, edit `NVCCFLAGS` to your compute capability (e.g. `sm_90` Hopper, `sm_89` Ada, `sm_80` Ampere) before building, or `gpu_test`/`cusparse_test` will fail to launch. The structural results are architecture-independent; time/energy scale with the board (§ Limitations).

### Build mm-repair CPU Baselines
`mm-repair` depends on [SDSL-lite](https://github.com/simongog/sdsl-lite) (packed `.iv` integer vectors); install it first — its own [`mm-repair/Readme.md`](mm-repair/Readme.md) lists the prerequisites and build steps. Then:
```bash
cd mm-repair
make clean
make all              # matrepair + csvmat2csrv + brepair/irepair0 + ANS/SDSL encoders
# For the tab:build timing column only, rebuild re32mm with detailed timing:
make re32mm CFLAGS="-Wall -std=c99 -g -O3 -DDETAILED_TIMING"
cd ..
```

---

## 2. Dataset Preparation

All datasets live in a single subfolder of this repo, `zenodo/`, laid out as
`zenodo/genotypes/`, `zenodo/wikidata/`, `zenodo/swh/` — this is exactly what `reproduce.sh`
reads (override the location with `ZENODO_DIR=/path/to/package`). There are two ways to
populate it:

- **Recommended — download the Zenodo data package** ([doi:10.5281/zenodo.22677746](https://doi.org/10.5281/zenodo.22677746),
  the concept DOI that always resolves to the latest version; the camera-ready used version
  [10.5281/zenodo.22677747](https://doi.org/10.5281/zenodo.22677747)). To stay under Zenodo's
  100-files-per-record limit the three data folders ship as **three uncompressed tar archives**
  (`genotypes.tar`, `wikidata.tar`, `swh.tar`; the payload is already compressed, so they are
  *not* gzipped). Extract all three in place to recreate the `zenodo/{genotypes,wikidata,swh}/`
  layout, then verify integrity against the shipped checksums:
  ```bash
  for t in genotypes.tar wikidata.tar swh.tar; do tar xf "$t"; done   # -> zenodo/{genotypes,wikidata,swh}/
  md5sum -c MANIFEST.md5                                              # expect: all files OK
  ```
  The package already contains every grammar (`.vc.C`, `.vc.R`, `.val`, …) and the
  Wikidata `.sparse` edge lists, so `./reproduce.sh crosscheck` / `struct` / `space` / `spmm` /
  `graph` run directly. The dense matrices ship zstd-compressed; decompress them only if you
  intend to *rebuild* a grammar (`./reproduce.sh grammar`) or run the from-scratch steps below:
  ```bash
  for f in zenodo/genotypes/*.zst zenodo/swh/*.zst; do zstd -d -k "$f"; done   # -> raw <base> alongside the grammar
  ```
  The three large `<base>.vc.zst` pre-RePair streams (`wd_country`, `wd_cites_work`, `swh_full`)
  need no manual handling: `./reproduce.sh graphscale` decompresses and back-dates them
  automatically as the lazy-rebuild gate for `matrepair -y` (§4, `tab:graph_scale`).
- **From scratch** — rebuild the same files into the same layout with the commands in §A–§E
  below. All matrices use the same dense-int32 format consumed by `mm-repair`.

If the datasets are already compiled and stored on the test server, you can skip this step.

**Python prerequisites:** `pip install numpy msprime==1.4.2` (msprime drives the synthetic genotype simulation in §B/§C; the VCF path in §A needs only `numpy` — `vcf2mat.py` reads the `.vcf.gz` with Python's standard-library `gzip`, no `pysam`/`cyvcf2`). The synthetic matrices are bit-for-bit reproducible only under the msprime version they were generated with — **`msprime==1.4.2`** (tskit 1.0.3, numpy 2.5.1); a different msprime release may change the coalescent RNG stream for the same seed.

> **§A–§E below are optional.** They document how each dataset in the package was derived
> from its public source, for transparency and independent regeneration. If you downloaded the
> Zenodo package you already have every file these steps produce — skip straight to §3/§4.

### A. Real Genotypes (1000 Genomes) — Chr20, Chr21, Chr22
The manuscript uses three human chromosomes, each at a $10^5$-variant subset and at full width. Download the phase-3 VCFs and slice them with `vcf2mat.py` into `zenodo/genotypes/`:
```bash
mkdir -p zenodo/genotypes
B=https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502
for chr in 20 21 22; do
  wget -O zenodo/genotypes/chr$chr.vcf.gz \
    $B/ALL.chr$chr.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz
done
# 100K-variant subset and full width (full width is capped high; the actual col count is printed):
python3 vcf2mat.py zenodo/genotypes/chr22.vcf.gz zenodo/genotypes/geno22     100000
python3 vcf2mat.py zenodo/genotypes/chr22.vcf.gz zenodo/genotypes/geno22full 1500000   # -> 1,055,454 cols
python3 vcf2mat.py zenodo/genotypes/chr21.vcf.gz zenodo/genotypes/geno21     100000
python3 vcf2mat.py zenodo/genotypes/chr21.vcf.gz zenodo/genotypes/geno21full 1500000   # -> 1,054,447 cols
python3 vcf2mat.py zenodo/genotypes/chr20.vcf.gz zenodo/genotypes/geno20     100000
python3 vcf2mat.py zenodo/genotypes/chr20.vcf.gz zenodo/genotypes/geno20full 2000000   # -> 1,739,315 cols
```
(`prepare_bio_datasets.py` is the authors' batch helper for the same steps, but its paths are hardcoded to the authors' `mm-grammar-gpu/` checkout — use the explicit commands above, or the Zenodo download.)

Resulting matrices and dimensions (rows = 2504 samples):

| base | rows | cols |
|---|---|---|
| `geno22`     | 2504 | 100000 |
| `geno22full` | 2504 | 1055454 |
| `geno21`     | 2504 | 100000 |
| `geno21full` | 2504 | 1054447 |
| `geno20`     | 2504 | 100000 |
| `geno20full` | 2504 | 1739315 |

### B. Synthetic Genotypes (Haplotypes)
The five synthetic configurations are simulated under the **coalescent with recombination** using [`msprime`](https://tskit.dev/msprime/) (`pip install msprime==1.4.2`), a standard, citable population-genetic simulator [Kelleher et al. 2016; Baumdicker et al. 2022]. Linkage disequilibrium is controlled by the recombination rate (low rate = long shared haplotype blocks = high LD = highly compressible), and every matrix is reproducible from a fixed `--seed`:
```bash
# args: <rows> <cols> <out_matrix> [recomb_rate] [seed] [Ne] [mu]   (seed=42 fixed for reproducibility)
python3 generate_msprime.py 2000  50000  zenodo/genotypes/geno_synth_small     1e-8 42   # synth_small
python3 generate_msprime.py 5000  200000 zenodo/genotypes/geno_synth_large     1e-8 42   # synth_large
python3 generate_msprime.py 5000  100000 zenodo/genotypes/geno_synth_ld_high   1e-9 42   # synth_ld_high (low recomb -> high LD)
python3 generate_msprime.py 5000  100000 zenodo/genotypes/geno_synth_ld_low    1e-7 42   # synth_ld_low  (high recomb -> low LD)
python3 generate_msprime.py 10000 50000  zenodo/genotypes/geno_synth_ind_large 1e-8 42   # synth_ind_large
```

### C. Large-scale / crossover matrix (Table 4.3, Figure 4.3)
The large crossover matrix (`crossover_synth`, $10{,}000 \times 700{,}000$, $\approx 1.00$ B nnz) is our billion-nonzero scale probe. Its CSR needs $\approx 8.0$ GB, which fits on the GB10's 119 GiB unified pool, so **cuSPARSE runs** (48.38 ms/vec, 1976 mJ/vec); the grammar engine stays resident at $1.01$ GB analytic / $0.98$ GB peak ($\approx 8.0\times$ smaller, $2.41\times$ faster, $3.15\times$ lower energy). Generated with the same `msprime` simulator as §B (this is a large run — the dense matrix is ~28 GB; produce it on the GB10 node's 119 GiB unified memory):
```bash
python3 generate_msprime.py 10000 700000 zenodo/genotypes/crossover_synth 1e-8 42
```

**`.val` value arrays (genotypes only).** `vcf2mat.py`/`generate_msprime.py` write the dense matrix but not the tiny `.val` file (the distinct-values array `{1.0, 2.0}` the engine reads for the `PlusTimes` semiring). Create one per genotype base (the Zenodo download already includes them):
```bash
for b in geno22 geno22full geno21 geno21full geno20 geno20full \
         geno_synth_small geno_synth_large geno_synth_ld_high geno_synth_ld_low \
         geno_synth_ind_large crossover_synth; do
  python3 -c "import struct,sys; open(sys.argv[1],'wb').write(struct.pack('dd',1.0,2.0))" zenodo/genotypes/$b.val
done
```

> ✅ **`crossover_synth` host-build cost (measured, node `spark-a459`).** The timed single-block grammar build
> ```bash
> ./mm-repair/matrepair -r --i32 zenodo/genotypes/crossover_synth 10000 700000
> ```
> takes **634.1 s total** (RePair 566.4 s + CSRV conversion 56.6 s) and yields a serialized grammar of **REANS 194.7 MB / RE32 250.7 MB** (the `>> REANS size` / `>> RE32 size` fields of the Compression Report; the lazy `-r -y` variant re-prints them without rebuilding). The space/energy/time figures above are read from the `## crossover_synth` block of `manuscript/logs/geno_space_energy.log` (produced by `reproduce.sh space`), and the run **self-verifies against the CPU reference within float precision** (`max_rel_diff` $\approx 9.4\mathrm{e}{-7}$; this is a $(+,\times)$ float run, so it is not bit-for-bit — that holds for the Boolean/Tropical semirings) — $10{,}000$ rows is below the driver's 30M-row CPU-ref threshold (`gpu_engine_test.cu:377`).

### D. Knowledge Graphs (Wikidata Relations)
The graph datasets are obtained from Zenodo: [10.5281/zenodo.7254968](https://zenodo.org/record/7254968) (Arroyuelo et al., *Datasets of Time- and Space-Efficient Regular Path Queries*).

**Wikidata (5 relations).** Download `wikidata.tar.gz` (4.95 GB), unpack to get `wikidata-enumerated.dat` (20 GB triples) + `.dat.P`. Wikidata relations have millions of subjects, so the **dense** builder is infeasible; use the **sparse** path (`sparse` mode → `matrepair --bool`), which never materializes the dense matrix. One pass emits all five `<name>.sparse` (`row col`) edge lists, then RePair compresses each:
```bash
mkdir -p zenodo/wikidata && tar -xzf wikidata.tar.gz

# 1. One-pass sparse extraction (predicate IDs are for this specific Zenodo dump):
python3 process_wikidata.py wikidata-enumerated.dat wikidata-enumerated.dat.P sparse \
  1107:wd_sports_team 708:wd_cast_member 205:wd_citizenship 206:wd_occupation 196:wd_subclass_of

# 2. Move each .sparse into the package and compress it into a grammar
#    (the base path == the .sparse; --bool = textual "row col" nonzeros):
for b in wd_sports_team wd_cast_member wd_citizenship wd_occupation wd_subclass_of; do
  mv -f $b.sparse zenodo/wikidata/$b.sparse; cp -f zenodo/wikidata/$b.sparse zenodo/wikidata/$b; done
./mm-repair/matrepair -r --bool zenodo/wikidata/wd_sports_team 332121  29854
./mm-repair/matrepair -r --bool zenodo/wikidata/wd_cast_member 173977  144095
./mm-repair/matrepair -r --bool zenodo/wikidata/wd_citizenship 2874250 2556
./mm-repair/matrepair -r --bool zenodo/wikidata/wd_occupation  3459933 10610
./mm-repair/matrepair -r --bool zenodo/wikidata/wd_subclass_of 1487709 73417
```
Wikidata relation dimensions (rows = subjects, cols = objects):

| base | Wikidata property | rows | cols | nnz |
|---|---|---:|---:|---:|
| `wd_sports_team` | P54 member of sports team | 332,121 | 29,854 | 1,136,249 |
| `wd_cast_member` | P161 cast member | 173,977 | 144,095 | 1,033,124 |
| `wd_citizenship` | P27 country of citizenship | 2,874,250 | 2,556 | 3,063,058 |
| `wd_occupation` | P106 occupation | 3,459,933 | 10,610 | 4,596,658 |
| `wd_subclass_of` | P279 subclass of | 1,487,709 | 73,417 | 2,024,347 |

**Large-scale relations (Table `tab:graph_scale`).** Two additional relations bracket the scale/compressibility spectrum (built the same way; the `sparse` extraction of `wd_cites_work` needs ~36 GB RAM to hold 166M edges). Their *dense* forms (2.2 GB / 10.8 TB) are unusable — the sparse path is what makes them addressable.
```bash
python3 process_wikidata.py wikidata-enumerated.dat wikidata-enumerated.dat.P sparse \
  159:wd_country 2831:wd_cites_work
for b in wd_country wd_cites_work; do
  mv -f $b.sparse zenodo/wikidata/$b.sparse; cp -f zenodo/wikidata/$b.sparse zenodo/wikidata/$b; done
./mm-repair/matrepair -r --bool zenodo/wikidata/wd_country    10058956 1747        # P17  country      (10.1M edges)
./mm-repair/matrepair -r --bool zenodo/wikidata/wd_cites_work 7072574  12245945    # P2860 cites-work  (166.7M edges)
```

| base | Wikidata property | rows | cols | nnz |
|---|---|---:|---:|---:|
| `wd_country` | P17 country | 10,058,956 | 1,747 | 10,089,284 |
| `wd_cites_work` | P2860 cites work | 7,072,574 | 12,245,945 | 166,682,725 |

### E. Software Heritage graph (SWH)

The software-provenance graph benchmark. Software Heritage's Merkle-DAG of software artifacts is **repetitive by construction** (deduplicated directories/contents shared across revisions) and ships in **WebGraph** format (the compressed-adjacency lineage already cited by the manuscript). We use the teaser export `2021-03-23-popular-3k-python` (3000 popular Python repositories) and build the forward artifact DAG as one **N×N boolean adjacency** ($N=45{,}691{,}499$ nodes, $1{,}218{,}488{,}928$ edges); its dense form (~261 TB) is impossible, so the sparse-to-grammar path is what makes it addressable. The forward DAG is extracted by `Dump.java` (below) against the WebGraph-big jars. Summary recipe (working dir `~/swh-work` on the node, kept out of the git tree; the node needed a portable JDK 21 + WebGraph-big jars):

```bash
cd ~/swh-work
B=https://softwareheritage.s3.amazonaws.com/graph/2021-03-23-popular-3k-python/compressed
# 1. Topology only (~430 MB): graph.graph/.offsets/.obl/.properties (+ the *.count.txt/stats)
for f in graph.graph graph.offsets graph.obl graph.properties; do curl -s -O "$B/$f"; done
# 2. Dump the forward DAG as a row-major "src dst" edge list, straight into the package
#    (WebGraph-big, memory-mapped). SW = <repo>/zenodo/swh :
REPO=/path/to/g-mm-repair    # absolute path to your checkout
SW=$REPO/zenodo/swh && mkdir -p "$SW" && cp -f graph.properties "$SW/"
export JAVA_HOME=~/swh-work/jdk && export PATH=$JAVA_HOME/bin:$PATH
javac -cp "jars/*" Dump.java && java -Xmx6g -cp "jars/*:." Dump graph "$SW/swh_full"   # 1.22 G edges, ~2 min
# 3. RePair-compress (no dense matrix; --bool = textual "row col" nonzeros):
$REPO/mm-repair/matrepair -r --bool "$SW/swh_full" 45691499 45691499                   # REANS 490 MB, 3.22 bpe
# 4. Structural check + Boolean device footprint:
export PATH=/usr/local/cuda/bin:$PATH
SEMIRING=boolean $REPO/gpu-engine/gpu_test "$SW/swh_full" 45691499 45691499 1 repair    # struct + device footprint (CPU ref auto-skipped >30M rows)
# The >30M skip is a speed guard, not a hard limit. To certify the *full* published graph
# bit-for-bit, force the CPU oracle on (slow, ~1.5 s/vector CPU seq, one-off):
FORCE_CPU_VERIFY=1 SEMIRING=boolean $REPO/gpu-engine/gpu_test "$SW/swh_full" 45691499 45691499 1 repair
# The crosscheck driver runs exactly this (Boolean, FORCE_CPU_VERIFY) plus cuSPARSE/CPU-CSR: ./reproduce.sh crosscheck
```
SWH dimensions and grammar/space figures (rows = cols = shared artifact node space):

| base | rows = cols | nnz | REANS | eng. device (analytic/peak) | verification |
|---|---|---|---|---|---|
| `swh_full`   | 45,691,499 | 1,218,488,928 | 490 MB (3.22 bpe) | 2.79 GB / 2.82 GB | **bit-for-bit exact** (Boolean; GPU vs sequential + OpenMP CPU sweeps, `max_abs_diff=0`) with `FORCE_CPU_VERIFY=1` (`./reproduce.sh crosscheck`) |

### F. Cross-implementation correctness (`./reproduce.sh crosscheck`)

Beyond each driver's own CPU-reference check, `crosscheck` certifies that **every** implementation
agrees with the engine on **one shared input vector** per dataset/semiring. The engine dumps the
vector `x` it used and its output `y` (`CROSSCHECK=<prefix>`); each independent tool then reloads that
same `x`, recomputes, and asserts agreement with the engine's `y`:

The columns of `tab:geno_time`, `tab:geno_spmm`, `tab:graph`, and `tab:graph_scale` are exactly these
implementations/configs:

| Semiring (tables) | Implementations checked against the engine (all must match) |
|---|---|
| `plustimes` (geno_time / geno_spmm) | CPU sequential sweep, CPU OpenMP-20 sweep, **`mm-repair` re32mm** CPU grammar mat-vec, cuSPARSE SpMV, CPU CSR SpMV, **batched engine + cuSPARSE SpMM** ($B{=}16$) |
| `boolean` (graph / graph_scale) | CPU seq, CPU OpenMP, cuSPARSE (count→bool), CPU CSR SpMV, SuiteSparse:GraphBLAS `lor_land`, batched engine ($B{=}16$) |
| `tropical` (graph) | CPU seq, CPU OpenMP, CPU CSR min-plus (same reconstructed matrix), SuiteSparse:GraphBLAS `min_plus`, batched engine ($B{=}16$) |

The cuSPARSE / CPU-CSR / re32mm paths reconstruct the matrix from the *same* grammar the engine reads (so
the matrix is identical by construction); GraphBLAS consumes the independent `.sparse` edge list. Run:

```bash
./reproduce.sh crosscheck        # -> manuscript/logs/crosscheck.log
```

Canonical result (`manuscript/logs/crosscheck.log`, node `spark-a459`): **25 blocks, 209 SUCCESS, 0 failures** —
all 12 genotype (`plustimes`), 5 Wikidata (`boolean`+`tropical`), and 2 largest-relation (`wd_country`,
`wd_cites_work`) rows pass, plus the full `swh_full` (1.22 G edges: engine bit-for-bit vs the sequential + OpenMP CPU sweeps and vs cuSPARSE / CPU CSR). Agreement is
within float precision for `plustimes` (`max_rel_diff ~1e-6`) and **bit-for-bit** (`mismatches=0`) for
`boolean`/`tropical`. The same CPU-sweep-vs-brute-force check runs GPU-free for all three semirings in CI
(`tests/run_integration_tests.py` → `cpu_selftest`).

**Choosing the input vector `x`.** By default the engine and each baseline draw `x` at random (the
engine then dumps its own `x` under `CROSSCHECK` so every implementation shares it — reproducibility does
*not* rely on a shared RNG seed). To pin a *specific* `x`, set `XVEC=<file>` on `gpu_test`: the file is a
raw binary array of exactly `cols` `float32` values — the **same format** the baselines consume in
`CROSSCHECK` mode (`<prefix>.<sr>.x`). `XVEC` composes with `CROSSCHECK`, so a hand-chosen `x` is dumped
and cross-validated across cuSPARSE / GraphBLAS / CPU automatically. Absent `XVEC`, behaviour is unchanged.

```bash
# little-endian float32, length = cols
python3 -c "import numpy; numpy.asarray([2,1,1,2,1,1,1],dtype='<f4').tofile('x.bin')"
XVEC=x.bin CROSSCHECK=out ./gpu-engine/gpu_test <base> <rows> <cols> 1 repair
```

---

## 3. Reproducing Tables and Figures

### Table 4.1 (Genotype structural figures) & Table B.2 (Wikidata structural figures)
The structural metrics (base rule count $|\mathcal{R}|$, depth $L$, maximum streaming width $w^{*}$, and pass-through nodes $+\text{pt}$) are computed by the host scheduler during graph loading and are architecture-independent. Run the test driver on each base path and read the layout from stdout:
```bash
# Table 4.1 — each genotype base (example shown for geno21):
./gpu-engine/gpu_test zenodo/genotypes/geno21 2504 100000 1

# Table B.2 — each Wikidata relation:
./gpu-engine/gpu_test zenodo/wikidata/wd_sports_team 332121 29854 1    # repeat for the 5 wd_* relations (dims above)
```
The driver prints one structural line, e.g.:
`Max depth (L): 6, Max width (w*): 49298, Total layered rules: 80896 (raw NTs: 70158, +pt: 10738)`
mapping to the table columns:
- `raw NTs`            -> $|\mathcal{R}|$
- `Max width (w*)`     -> $w^{*}$
- `Max depth (L)`      -> $L$
- `+pt`                -> $+\text{pt}$ (and `Total layered rules` = $|\mathcal{R}| + \text{pt}$)

### Table 4.2: Genotype Average Time (ms/vector)
Run the complete bioinformatics benchmarking suite (all 11 real + synthetic genotype datasets):
```bash
bash run_all_bio_baselines.sh
```
This benchmarks each dataset across:
1. GPU Level-Synchronous Sweep
2. CPU Level-Synchronous Sweep (OpenMP, 20 threads)
3. CPU Level-Synchronous Sweep (Sequential)
4. mm-repair (Sequential)
5. mm-repair (16 threads)

and the cuSPARSE CSR SpMV column. Results are saved to `manuscript/logs/bio_results.csv`.

### Table 4.3: Genotype Space & Energy vs. cuSPARSE (incl. crossover)
Execute the space/energy benchmark over the genotype matrices plus `crossover_synth`:
```bash
./reproduce.sh space
```
This runs `gpu_test` and `cusparse_test` per dataset (analytic device bytes, measured peak bytes, time, and GPU energy mJ/vector, grammar engine vs. cuSPARSE CSR) and writes the canonical log `manuscript/logs/geno_space_energy.log` — the file `extract_results.py` reads for this table. (The standalone `run_space_energy.py` driver is an alternate front-end that writes `manuscript/logs/space_energy_raw.txt`; it is *not* the canonical log and is not consumed by the extractor.) `crossover_synth` is the billion-nnz scale probe: cuSPARSE's CSR ($\approx 8.03$ GB) now fits and runs, while the engine reports its $1.01$ GB analytic / $0.98$ GB peak footprint ($\approx 8.0\times$ smaller).

### Table A.1: Batched Right Product (SpMM, $Y=MX$)
Run the GPU engine in batched mode (7th arg $B$) and cuSPARSE SpMM, sweeping batch sizes ($B = 16, 32, 64, 128, 256$, as in `reproduce.sh spmm`) and recording the best per-vector time:
```bash
# Grammar engine, batched (example: geno22full, B=32):
./gpu-engine/gpu_test zenodo/genotypes/geno22full 2504 1055454 100 repair 32

# cuSPARSE SpMM — sweeps ALG_DEFAULT / CSR_ALG2 / CSR_ALG3; CSR_ALG3 is the non-degrading baseline:
./gpu-engine/cusparse_test zenodo/genotypes/geno22full 2504 1055454 100 repair 32
```
Repeat for all 6 real + 5 synthetic genotype matrices (dimensions per §A/§B).

### Table 5.2: Graph Right Product (Boolean & Tropical Semirings)
Set the `SEMIRING` environment variable and run `gpu_test` (single-vector and batched $B{=}16$) on **each** of the five Wikidata relation base paths:
```bash
# --- Wikidata ---
for sr in boolean tropical; do
  SEMIRING=$sr ./gpu-engine/gpu_test zenodo/wikidata/wd_sports_team 332121  29854  50 repair 16
  SEMIRING=$sr ./gpu-engine/gpu_test zenodo/wikidata/wd_cast_member 173977  144095 50 repair 16
  SEMIRING=$sr ./gpu-engine/gpu_test zenodo/wikidata/wd_citizenship 2874250 2556   50 repair 16
  SEMIRING=$sr ./gpu-engine/gpu_test zenodo/wikidata/wd_occupation  3459933 10610  50 repair 16
  SEMIRING=$sr ./gpu-engine/gpu_test zenodo/wikidata/wd_subclass_of 1487709 73417  50 repair 16
done
```
The CPU reference (sequential and OpenMP) and the cuSPARSE Boolean CSR baseline reported in the table are produced by the same driver; results are verified bit-for-bit against the CPU reference (all five relations pass).

### Table 5.2 baselines: GraphBLAS (CPU) and cuGraph (GPU)
Because no vendor *dense* kernel exists for the Boolean/Tropical semirings, the graph baselines are the semiring-native **SuiteSparse:GraphBLAS** (CPU) and **cuGraph** (GPU BFS/SSSP). Both consume the same `<name>.sparse` edge lists (`row col`) produced by the `sparse` mode in §2.D.

```bash
# --- SuiteSparse:GraphBLAS (semiring-native mxv: lor_land / min_plus; 1 and 20 threads) ---
python3 -m venv gbvenv && ./gbvenv/bin/pip install python-graphblas
# Automated execution of all relations:
chmod +x gb_run.sh && ./gb_run.sh > manuscript/logs/graphblas_bench.log 2>&1
# Or run manually:
# args: <sparse_file> <rows> <cols> <bool|tropical> [iters]
./gbvenv/bin/python graphblas_bench.py zenodo/wikidata/wd_sports_team.sparse 332121 29854 bool 50
./gbvenv/bin/python graphblas_bench.py zenodo/wikidata/wd_sports_team.sparse 332121 29854 tropical 50
#   (repeat for the 5 Wikidata relations; dims per §2.D)

# --- cuGraph (GPU end-to-end BFS/SSSP; RAPIDS cu12 wheels) ---
python3 -m venv cgvenv
./cgvenv/bin/pip install --extra-index-url=https://pypi.nvidia.com cugraph-cu12 cudf-cu12
# Automated execution of all relations:
chmod +x cg_run.sh && ./cg_run.sh > manuscript/logs/cugraph_bench.log 2>&1
# Or run manually (RAPIDS pip wheels need their bundled libs on the loader path):
export LD_LIBRARY_PATH="$(find cgvenv/lib/python3.12/site-packages -type d \( -name lib -o -name lib64 \) | tr '\n' ':'):/usr/local/cuda/lib64"
./cgvenv/bin/python cugraph_bench.py zenodo/wikidata/wd_sports_team.sparse 332121 29854 bool 20      # BFS
./cgvenv/bin/python cugraph_bench.py zenodo/wikidata/wd_sports_team.sparse 332121 29854 tropical 20  # SSSP
```
Note: GraphBLAS `mxv` is one semiring mat-vec (directly comparable to the engine's single-vector sweep), whereas cuGraph BFS/SSSP run the *full* traversal to convergence (an end-to-end reference, not per mat-vec). On the node, RAPIDS is the CUDA-12 build running on CUDA 13 / `sm_121` via PTX-JIT.

### Table 5.1: Scale at 10M–1.2G edges (`tab:graph_scale`)
The two largest Wikidata relations (`wd_country`, `wd_cites_work`; built in §2.D) **plus the billion-edge Software Heritage software graph** `swh_full` (§2.E). Per relation: the Boolean single-vector engine and cuSPARSE times, the engine/CSR analytic device footprints, and the serialized REANS grammar size. The `matrepair -r -y` call is *lazy* — `-y` skips recompression and just prints the size report (`>> REANS size: N bytes`) from the existing grammar, so it is fast:
```bash
# All three live in the Zenodo package layout (zenodo/wikidata, zenodo/swh):
for e in "zenodo/wikidata/wd_country 10058956 1747" "zenodo/wikidata/wd_cites_work 7072574 12245945" \
         "zenodo/swh/swh_full 45691499 45691499"; do set -- $e
  SEMIRING=boolean ./gpu-engine/gpu_test  "$1" "$2" "$3" 50 repair   # struct, eng dev MB, Bool eng ms
  ./gpu-engine/cusparse_test              "$1" "$2" "$3" 50          # CSR dev MB, Bool cuS ms
  ./mm-repair/matrepair -r -y --bool      "$1" "$2" "$3"             # serialized REANS size
done
```
`reproduce.sh graphscale` runs exactly this into `manuscript/logs/graph_scale.log`. It is **not** part of `./reproduce.sh all` because `wd_cites_work` (166.7M edges) and `swh_full` (1.22G edges) are heavy runs. **Note (SWH row):** in `tab:graph_scale` the space columns (REANS, CSR, eng.\ dev.) are analytic/deterministic, while the two time columns (Bool eng/cuS ms) carry the usual profiling variance of the shared unified-memory node. The full graph exceeds the driver's 30M-row CPU-reference threshold (`gpu_engine_test.cu`) so it does not self-verify by default; `./reproduce.sh crosscheck` runs it with `FORCE_CPU_VERIFY=1`, certifying the full graph bit-for-bit against the sequential + OpenMP CPU sweeps (and against cuSPARSE / CPU-CSR SpMV) — see §2.F.

---

## 4. From logs to tables and figures (automated extraction)

The manuscript tables and the two measured figures are **generated from the logs**, not transcribed by hand. The chain is:

```
run experiment  ->  manuscript/logs/<experiment>.log   (raw, provenance-headed, "## <key>" blocks)
                ->  extract_results.py       ->  manuscript/tables/tab_*.tex   (\input by main.tex)
                                                 manuscript/figures/data/*.dat (read by TikZ)
                ->  pdflatex fig_*.tex        ->  manuscript/figures/fig_*.pdf
```

**One command (on the GB10 node)** runs the core experiment families into their canonical logs, then extracts and plots:
```bash
./reproduce.sh all      # = struct + time + space + spmm + graph, then extract + plot
```
Or run a single stage: `./reproduce.sh {struct|time|space|spmm|graph}` re-runs one experiment family; `./reproduce.sh extract` re-derives all tables/figure-data from the **existing** logs without recomputing; `./reproduce.sh plot` recompiles the TikZ figures.

Three heavy families are **not** part of `all` and are run on demand: `./reproduce.sh grammar` (→ `grammar_build.log`, Table B.1 grammar column), `./reproduce.sh graphscale` (→ `graph_scale.log`, Table 5.1), and `./reproduce.sh crosscheck` (→ `crosscheck.log`, cross-implementation correctness). The GraphBLAS and cuGraph baseline logs are produced separately by `gb_run.sh` / `cg_run.sh` (§ *Table 5.2 baselines*).

### Canonical logs (kept under `manuscript/logs/`, gitignored)
Every log begins with a provenance header (date, host, git commit) and separates datasets with `## <key>` markers so a single parser can slice it. Log → artifact map:

| Canonical log | Feeds | Produced by |
|---|---|---|
| `manuscript/logs/geno_struct.log` | Table 4.1 (`tab:geno_through`) | `reproduce.sh struct` |
| `manuscript/logs/bio_results.csv` + `manuscript/logs/geno_space_energy.log` | Table 4.2 (`tab:geno_time`) | `reproduce.sh time`/`space` |
| `manuscript/logs/geno_space_energy.log` | Table 4.3 (`tab:geno`), Table B.1 (`tab:build`, $+$pt col), Fig. 4.3 (`fig:space`) | `reproduce.sh space` |
| `manuscript/logs/grammar_build.log` | Table B.1 (`tab:build`, grammar col) | `reproduce.sh grammar` |
| `manuscript/logs/geno_spmm.log` + `manuscript/logs/geno_cusparse_alg.log` | Table A.1 (`tab:geno_spmm`), Fig. A.1 (`fig:batched`) | `reproduce.sh spmm` |
| `manuscript/logs/graph_struct.log` | Table B.2 (`tab:graph_struct`) | `reproduce.sh struct` |
| `manuscript/logs/graph_semiring.log` | Table 5.2 (`tab:graph`) | `reproduce.sh graph` |
| `manuscript/logs/graph_scale.log` | Table 5.1 (`tab:graph_scale`) | `reproduce.sh graphscale` |

The following logs are **not** parsed by `extract_results.py` (they back correctness claims and baselines, not table cells), but are kept alongside the canonical logs for completeness:

| Log | Role | Produced by |
|---|---|---|
| `manuscript/logs/crosscheck.log` | Cross-implementation correctness (§2.F); certifies the full SWH graph bit-for-bit | `reproduce.sh crosscheck` |
| `manuscript/logs/graphblas_bench.log` | Standalone SuiteSparse:GraphBLAS baseline (Table 5.2's GB numbers are extracted from `graph_semiring.log`, which embeds the same run) | `gb_run.sh` |
| `manuscript/logs/cugraph_bench.log` | cuGraph BFS/SSSP end-to-end reference, discussed in the manuscript (§ Limitations) | `cg_run.sh` |

### Host-side construction cost (Table B.1 / `tab:build`)
`tab:build` reports **two distinct, additive** host construction costs per genotype matrix:

1. **grammar RePair (s)** — the one-time offline `mm-repair` compressor build, *shared with the CPU baseline*. `matrepair -r` prints a "Compression Report" (`>> total time: X`); `reproduce.sh grammar` runs it per dataset into `manuscript/logs/grammar_build.log`. This recompresses the grammars, so it is **expensive and kept out of `./reproduce.sh all`** — run it once to populate the column.

   **Decision — this build is sequential, single-block, by design.** In `mm-repair`, the `-p` flag parallelizes `irepair` *only across `-b` row-blocks* and is hard-capped at `-b` (`matrepair` lines 75–77: `if args.p > args.b: args.p = args.b`). Passing `-p 16` to a single-block build (`-b 1`) is silently reduced to `-p 1`; the underlying `irepair0` is a single sequential process, not internally threaded. The genotype grammars are **single-block** — that is the exact grammar the engine consumes (the driver reads one `.vc.C`) and the one whose `|R|`, `L`, `w*`, `+pt` are reported in every other table. We therefore time the **sequential single-block** construction: `matrepair -r --i32 <path> <rows> <cols>` (`-b 1`, `-p 1`). Using `-b > 1` to unlock parallelism would build a *structurally different* multi-block grammar, no longer the object benchmarked elsewhere — so we do not.
2. **+pt build (ms)** — the engine's *marginal* completion step (pass-through completion + level-bucketing + terminal compaction, i.e. `build_schedule`), timed inside `gpu_test`, which prints a canonical line every run:
   ```
   BUILD schedule ms: <t> peak host MB: <m>
   ```
   `extract_results.py` reads it from `manuscript/logs/geno_space_energy.log` (engine half of each `## <key>` block). No extra experiment is needed — it comes free from `reproduce.sh space`; just rebuild `gpu_test` first so the `BUILD` line is present.

The **≈ mat-vecs** column is `+pt build / eval` (single right-product time, `tab:geno`): the number of products after which the *marginal* cost is amortized — a handful for every matrix. The grammar cost is the one-time price of the compressed format, not part of that threshold. To populate the full table: `./reproduce.sh grammar && ./reproduce.sh space && python3 extract_results.py`.

### Extraction only
```bash
python3 extract_results.py          # manuscript/logs/ -> manuscript/tables/*.tex + figures/data/*.dat
```
Missing logs are skipped with a warning (never a hard error), so the pipeline can be run incrementally as experiments land. A regenerated table normalizes number formatting but carries the same measured values; the existing `manuscript/tables/tab_*.tex` and `figures/data/*.dat` are the last-known-good fallbacks until a fresh node run overwrites them.

### Figures (TikZ / pgfplots, not matplotlib)
- **Figure 4.3** (`fig:space`) and **Figure A.1** (`fig:batched`) are standalone pgfplots sources — `manuscript/figures/fig_space.tex` and `fig_batched.tex` — that read `manuscript/figures/data/fig_space.dat` / `fig_batched.dat`. They compile to `manuscript/figures/fig_space.pdf` / `fig_batched.pdf`, which `main.tex` includes via `\includegraphics` (same convention as the schematic figures).
```bash
./reproduce.sh plot                 # or: cd manuscript/figures && pdflatex fig_space.tex
```
