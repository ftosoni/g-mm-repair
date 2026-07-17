# Reproducibility Guide for "Streaming Right Multiplication over Grammar-Compressed Matrices"

This document provides step-by-step instructions to reproduce **all and only** the tables and figures presented in the manuscript. The manuscript reports two experiment families: **genotype matrices** (1000 Genomes + synthetic haplotypes, `PlusTimes` semiring) and **knowledge-graph relation matrices** (Wikidata relations, `Boolean`/`Tropical` semirings). No machine-learning matrices (e.g. `higgs`, `covtype`) are used in the manuscript and they are intentionally excluded here.

These experiments were originally profiled and measured on a remote prototyping node (**NVIDIA Grace-Blackwell GB10 node**, featuring unified coherent CPU-GPU memory of 119 GiB, CUDA 13.0, g++ 13.3, and Ubuntu 24.04 LTS). Per the manuscript (§ Limitations), time/energy figures are board-dependent; the structural figures ($|\mathcal{R}|$, $L$, $w^{*}$, $+\text{pt}$) are architecture-independent and reproducible on any host.

> ✅ **Status.** Every measured table and figure in this guide is backed by a canonical log under `manuscript/logs/`. The genotype/crossover family (§A–§C, Tables 1–4, B, Figs 3–4) was regenerated end-to-end with the `msprime` simulator on node `spark-a459` (2026-07-14); the Wikidata/SWH graph logs (§D–§E, Tables 5–7) are unchanged.

The complete inventory of manuscript artifacts this guide reproduces:

| Artifact | Content | Datasets |
|---|---|---|
| Table 1 (`tab:geno_through`) | Structural figures, genotypes | 6 real + 5 synthetic genotype matrices |
| Table 2 (`tab:geno_time`) | Avg time/vector, genotypes | same 11 genotype matrices |
| Table 3 (`tab:geno`) | Space & energy vs. cuSPARSE | same 11 + `crossover_synth` |
| Table 4 (`tab:geno_spmm`) | Batched right product (SpMM) | 6 real + 5 synthetic genotype matrices |
| Table 5 (`tab:graph_struct`) | Structural figures, Wikidata | 5 Wikidata relation matrices |
| Table 6 (`tab:graph`) | Boolean & Tropical graph product | 5 Wikidata relation matrices |
| Table 7 (`tab:graph_scale`) | Scale/space at 10M–1.2G edges | 2 largest Wikidata relations + SWH software graph |
| Table B (`tab:build`) | Host-side construction cost + amortization | 11 genotype matrices |
| Figure 3 (`fig:space`) | Memory footprint (incl. `crossover_synth` at billion-nnz scale) | genotype matrices + crossover |
| Figure 4 (`fig:batched`) | Batched throughput vs. $B$ | `geno22full` |

(The remaining figures — `fig:matrix`, `fig:mmr_rs`, `fig:dag`, `fig:completion`, `fig:sweep`, `fig:trace` — are schematic TikZ/`includegraphics` illustrations of the running example, not measured results, and require no experiment to regenerate.)

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

### Build mm-repair CPU Baselines
```bash
cd mm-repair
make clean
# Build re32mm with detailed timing enabled
make re32mm CFLAGS="-Wall -std=c99 -g -O3 -DDETAILED_TIMING"
cd ..
```

---

## 2. Dataset Preparation

If datasets are already compiled and stored on the test server, you can skip this step. Otherwise, follow these instructions. All matrices use the same dense-int32 format consumed by `mm-repair`.

**Python prerequisites:** `pip install numpy msprime` (msprime drives the synthetic genotype simulation in §B/§C; the VCF path in §A additionally needs `pysam`/`cyvcf2` per `vcf2mat.py`).

### A. Real Genotypes (1000 Genomes) — Chr20, Chr21, Chr22
The manuscript uses three human chromosomes, each at a $10^5$-variant subset and at full width. `prepare_bio_datasets.py` automates Chr20/Chr21 (and the synthetic sets in §B); reproduce **Chr22** the same way. Download the phase-3 VCFs and slice them with `vcf2mat.py`:
```bash
# Chr22 (running-example chromosome; geno22 / geno22full)
wget -O geno/chr22.vcf.gz \
  https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/ALL.chr22.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz
python3 geno/vcf2mat.py geno/chr22.vcf.gz mm-repair/data/geno22     100000
python3 geno/vcf2mat.py geno/chr22.vcf.gz mm-repair/data/geno22full 1500000   # yields 1,055,454 cols

# Chr20 + Chr21 + the synthetic genotypes of §B (downloads VCFs, runs vcf2mat, writes .val files):
python3 prepare_bio_datasets.py
```
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
The five synthetic configurations are simulated under the **coalescent with recombination** using [`msprime`](https://tskit.dev/msprime/) (`pip install msprime`), a standard, citable population-genetic simulator [Kelleher et al. 2016; Baumdicker et al. 2022]. Linkage disequilibrium is controlled by the recombination rate (low rate = long shared haplotype blocks = high LD = highly compressible), and every matrix is reproducible from a fixed `--seed`. `prepare_bio_datasets.py` runs all five; the explicit commands are:
```bash
# args: <rows> <cols> <out_matrix> [recomb_rate] [seed] [Ne] [mu]   (seed=42 fixed for reproducibility)
python3 generate_msprime.py 2000  50000  mm-repair/data/geno_synth_small     1e-8 42   # synth_small
python3 generate_msprime.py 5000  200000 mm-repair/data/geno_synth_large     1e-8 42   # synth_large
python3 generate_msprime.py 5000  100000 mm-repair/data/geno_synth_ld_high   1e-9 42   # synth_ld_high (low recomb -> high LD)
python3 generate_msprime.py 5000  100000 mm-repair/data/geno_synth_ld_low    1e-7 42   # synth_ld_low  (high recomb -> low LD)
python3 generate_msprime.py 10000 50000  mm-repair/data/geno_synth_ind_large 1e-8 42   # synth_ind_large
```

### C. Large-scale / crossover matrix (Table 3, Figure 3)
The large crossover matrix (`crossover_synth`, $10{,}000 \times 700{,}000$, $\approx 1.00$ B nnz) is our billion-nonzero scale probe. Its CSR needs $\approx 8.0$ GB, which fits on the GB10's 119 GiB unified pool, so **cuSPARSE runs** (48.38 ms/vec, 1976 mJ/vec); the grammar engine stays resident at $1.01$ GB analytic / $0.98$ GB peak ($\approx 8.0\times$ smaller, $2.41\times$ faster, $3.15\times$ lower energy). Generated with the same `msprime` simulator as §B (this is a large run — the dense matrix is ~28 GB; produce it on the GB10 node's 119 GiB unified memory):
```bash
python3 generate_msprime.py 10000 700000 mm-repair/data/crossover_synth 1e-8 42
```

> ✅ **`crossover_synth` host-build cost (measured, node `spark-a459`).** The timed single-block grammar build
> ```bash
> ./mm-repair/matrepair -r --i32 mm-repair/data/crossover_synth 10000 700000
> ```
> takes **634.1 s total** (RePair 566.4 s + CSRV conversion 56.6 s) and yields a serialized grammar of **REANS 194.7 MB / RE32 250.7 MB** (the `>> REANS size` / `>> RE32 size` fields of the Compression Report; the lazy `-r -y` variant re-prints them without rebuilding). The space/energy/time figures above are read from the `## crossover_synth` block of `manuscript/logs/geno_space_energy.log` (produced by `reproduce.sh space`), and the run **self-verifies bit-for-bit** against the CPU reference — $10{,}000$ rows is below the driver's 30M-row CPU-ref threshold (`gpu_engine_test.cu:377`).

### D. Knowledge Graphs (Wikidata Relations)
The graph datasets are obtained from Zenodo: [10.5281/zenodo.7254968](https://zenodo.org/record/7254968) (Arroyuelo et al., *Datasets of Time- and Space-Efficient Regular Path Queries*).

**Wikidata (5 relations).** Download `wikidata.tar.gz` (4.95 GB), unpack to get `wikidata-enumerated.dat` (20 GB triples) + `.dat.P`. Wikidata relations have millions of subjects, so the **dense** builder is infeasible; use the **sparse** path (`sparse` mode → `matrepair --bool`), which never materializes the dense matrix. One pass emits all five `<name>.sparse` (`row col`) edge lists, then RePair compresses each:
```bash
tar -xzf wikidata.tar.gz

# 1. One-pass sparse extraction (predicate IDs are for this specific Zenodo dump):
python3 process_wikidata.py wikidata-enumerated.dat wikidata-enumerated.dat.P sparse \
  1107:wd_sports_team 708:wd_cast_member 205:wd_citizenship 206:wd_occupation 196:wd_subclass_of

# 2. Compress each into a grammar (rename to a clean base, then --bool = textual "row col" nonzeros):
for b in wd_sports_team wd_cast_member wd_citizenship wd_occupation wd_subclass_of; do cp -f $b.sparse $b; done
./mm-repair/matrepair -r --bool wd_sports_team 332121  29854
./mm-repair/matrepair -r --bool wd_cast_member 173977  144095
./mm-repair/matrepair -r --bool wd_citizenship 2874250 2556
./mm-repair/matrepair -r --bool wd_occupation  3459933 10610
./mm-repair/matrepair -r --bool wd_subclass_of 1487709 73417
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
for b in wd_country wd_cites_work; do cp -f $b.sparse $b; done
./mm-repair/matrepair -r --bool wd_country    10058956 1747        # P17  country      (10.1M edges)
./mm-repair/matrepair -r --bool wd_cites_work 7072574  12245945    # P2860 cites-work  (166.7M edges)
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
# 2. Dump the forward DAG as a row-major "src dst" edge list (WebGraph-big, memory-mapped):
export JAVA_HOME=~/swh-work/jdk && export PATH=$JAVA_HOME/bin:$PATH
javac -cp "jars/*" Dump.java && java -Xmx6g -cp "jars/*:." Dump graph swh_full     # 1.22 G edges, ~2 min
# 3. RePair-compress (no dense matrix; --bool = textual "row col" nonzeros):
REPO=/mnt/nfs/home/tosoni/mm-grammar-gpu
$REPO/mm-repair/matrepair -r --bool swh_full 45691499 45691499                      # REANS 490 MB, 3.22 bpe
# 4. Structural check + Boolean device footprint:
export PATH=/usr/local/cuda/bin:$PATH
SEMIRING=boolean $REPO/gpu-engine/gpu_test swh_full 45691499 45691499 1 repair       # struct + device footprint (CPU ref auto-skipped >30M rows)
# The >30M skip is a speed guard, not a hard limit. To certify the *full* published graph
# bit-for-bit, force the CPU oracle on (slow, ~1.5 s/vector CPU seq, one-off):
FORCE_CPU_VERIFY=1 SEMIRING=boolean $REPO/gpu-engine/gpu_test swh_full 45691499 45691499 1 repair
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

### Table 1 (Genotype structural figures) & Table 5 (Wikidata structural figures)
The structural metrics (base rule count $|\mathcal{R}|$, depth $L$, maximum streaming width $w^{*}$, and pass-through nodes $+\text{pt}$) are computed by the host scheduler during graph loading and are architecture-independent. Run the test driver on each base path and read the layout from stdout:
```bash
# Table 1 — each genotype base (example shown for geno21):
./gpu-engine/gpu_test mm-repair/data/geno21 2504 100000 1

# Table 5 — each Wikidata relation:
./gpu-engine/gpu_test wd_sports_team 332121 29854 1        # repeat for the 5 wd_* relations (dims above)
```
The driver prints one structural line, e.g.:
`Max depth (L): 6, Max width (w*): 49298, Total layered rules: 80896 (raw NTs: 70158, +pt: 10738)`
mapping to the table columns:
- `raw NTs`            -> $|\mathcal{R}|$
- `Max width (w*)`     -> $w^{*}$
- `Max depth (L)`      -> $L$
- `+pt`                -> $+\text{pt}$ (and `Total layered rules` = $|\mathcal{R}| + \text{pt}$)

### Table 2: Genotype Average Time (ms/vector)
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

### Table 3: Genotype Space & Energy vs. cuSPARSE (incl. crossover)
Execute the space/energy benchmark over the genotype matrices plus `crossover_synth`:
```bash
python3 run_space_energy.py
```
This compiles the comparative table (analytic device bytes, measured peak bytes, time, and GPU energy mJ/vector) for the grammar engine vs. cuSPARSE CSR, and stores raw results in `manuscript/logs/space_energy_results.txt`. **Note:** ensure the `datasets` list inside `run_space_energy.py` is set to the 11 genotype matrices plus `crossover_synth` (rows/cols per the tables above) — the script must not be pointed at any non-manuscript matrices. `crossover_synth` is the billion-nnz scale probe: cuSPARSE's CSR ($\approx 10.3$ GB) now fits and runs, while the engine reports its $1.30$ GB analytic / $1.59$ GB peak footprint ($\approx 7.9\times$ smaller).

### Table 4: Batched Right Product (SpMM, $Y=MX$)
Run the GPU engine in batched mode (7th arg $B$) and cuSPARSE SpMM, sweeping batch sizes ($B = 32, 64, 128, 256$) and recording the best per-vector time:
```bash
# Grammar engine, batched (example: geno22full, B=32):
./gpu-engine/gpu_test mm-repair/data/geno22full 2504 1055454 100 repair 32

# cuSPARSE SpMM — sweeps ALG_DEFAULT / CSR_ALG2 / CSR_ALG3; CSR_ALG3 is the non-degrading baseline:
./gpu-engine/cusparse_test mm-repair/data/geno22full 2504 1055454 100 repair 32
```
Repeat for all 6 real + 5 synthetic genotype matrices (dimensions per §A/§B).

### Table 6: Graph Right Product (Boolean & Tropical Semirings)
Set the `SEMIRING` environment variable and run `gpu_test` (single-vector and batched $B{=}16$) on **each** of the five Wikidata relation base paths:
```bash
# --- Wikidata ---
for sr in boolean tropical; do
  SEMIRING=$sr ./gpu-engine/gpu_test wd_sports_team 332121  29854  50 repair 16
  SEMIRING=$sr ./gpu-engine/gpu_test wd_cast_member 173977  144095 50 repair 16
  SEMIRING=$sr ./gpu-engine/gpu_test wd_citizenship 2874250 2556   50 repair 16
  SEMIRING=$sr ./gpu-engine/gpu_test wd_occupation  3459933 10610  50 repair 16
  SEMIRING=$sr ./gpu-engine/gpu_test wd_subclass_of 1487709 73417  50 repair 16
done
```
The CPU reference (sequential and OpenMP) and the cuSPARSE Boolean CSR baseline reported in the table are produced by the same driver; results are verified bit-for-bit against the CPU reference (all five relations pass).

### Table 6 baselines: GraphBLAS (CPU) and cuGraph (GPU)
Because no vendor *dense* kernel exists for the Boolean/Tropical semirings, the graph baselines are the semiring-native **SuiteSparse:GraphBLAS** (CPU) and **cuGraph** (GPU BFS/SSSP). Both consume the same `<name>.sparse` edge lists (`row col`) produced by the `sparse` mode in §2.D.

```bash
# --- SuiteSparse:GraphBLAS (semiring-native mxv: lor_land / min_plus; 1 and 20 threads) ---
python3 -m venv gbvenv && ./gbvenv/bin/pip install python-graphblas
# Automated execution of all relations:
chmod +x gb_run.sh && ./gb_run.sh > manuscript/logs/graphblas_bench.log 2>&1
# Or run manually:
# args: <sparse_file> <rows> <cols> <bool|tropical> [iters]
./gbvenv/bin/python graphblas_bench.py wd_sports_team.sparse 332121 29854 bool 50
./gbvenv/bin/python graphblas_bench.py wd_sports_team.sparse 332121 29854 tropical 50
#   (repeat for the 5 Wikidata relations; dims per §2.D)

# --- cuGraph (GPU end-to-end BFS/SSSP; RAPIDS cu12 wheels) ---
python3 -m venv cgvenv
./cgvenv/bin/pip install --extra-index-url=https://pypi.nvidia.com cugraph-cu12 cudf-cu12
# Automated execution of all relations:
chmod +x cg_run.sh && ./cg_run.sh > manuscript/logs/cugraph_bench.log 2>&1
# Or run manually (RAPIDS pip wheels need their bundled libs on the loader path):
export LD_LIBRARY_PATH="$(find cgvenv/lib/python3.12/site-packages -type d \( -name lib -o -name lib64 \) | tr '\n' ':'):/usr/local/cuda/lib64"
./cgvenv/bin/python cugraph_bench.py wd_sports_team.sparse 332121 29854 bool 20      # BFS
./cgvenv/bin/python cugraph_bench.py wd_sports_team.sparse 332121 29854 tropical 20  # SSSP
```
Note: GraphBLAS `mxv` is one semiring mat-vec (directly comparable to the engine's single-vector sweep), whereas cuGraph BFS/SSSP run the *full* traversal to convergence (an end-to-end reference, not per mat-vec). On the node, RAPIDS is the CUDA-12 build running on CUDA 13 / `sm_121` via PTX-JIT.

### Table 7: Scale at 10M–1.2G edges (`tab:graph_scale`)
The two largest Wikidata relations (`wd_country`, `wd_cites_work`; built in §2.D) **plus the billion-edge Software Heritage software graph** `swh_full` (§2.E). Per relation: the Boolean single-vector engine and cuSPARSE times, the engine/CSR analytic device footprints, and the serialized REANS grammar size. The `matrepair -r -y` call is *lazy* — `-y` skips recompression and just prints the size report (`>> REANS size: N bytes`) from the existing grammar, so it is fast:
```bash
# Wikidata relations live in the repo dir; swh_full lives under ~/swh-work (absolute path):
for e in "wd_country 10058956 1747" "wd_cites_work 7072574 12245945" \
         "$HOME/swh-work/swh_full 45691499 45691499"; do set -- $e
  SEMIRING=boolean ./gpu-engine/gpu_test  "$1" "$2" "$3" 50 repair   # struct, eng dev MB, Bool eng ms
  ./gpu-engine/cusparse_test              "$1" "$2" "$3" 50          # CSR dev MB, Bool cuS ms
  ./mm-repair/matrepair -r -y --bool      "$1" "$2" "$3"             # serialized REANS size
done
```
`reproduce.sh graphscale` runs exactly this into `manuscript/logs/graph_scale.log`. It is **not** part of `./reproduce.sh all` because `wd_cites_work` (166.7M edges) and `swh_full` (1.22G edges) are heavy runs. **Note (SWH row):** the current `tab:graph_scale` SWH figures were measured under multi-tenant node load, so the two time columns (Bool eng/cuS ms) are **preliminary**; the space columns (REANS, CSR, eng.\ dev.) are analytic/deterministic. The full graph exceeds the driver's 30M-row CPU-reference threshold (`gpu_engine_test.cu`) so it does not self-verify by default; `./reproduce.sh crosscheck` runs it with `FORCE_CPU_VERIFY=1`, certifying the full graph bit-for-bit against the sequential + OpenMP CPU sweeps (and against cuSPARSE / CPU-CSR SpMV) — see §2.F.

---

## 4. From logs to tables and figures (automated extraction)

The manuscript tables and the two measured figures are **generated from the logs**, not transcribed by hand. The chain is:

```
run experiment  ->  manuscript/logs/<experiment>.log   (raw, provenance-headed, "## <key>" blocks)
                ->  extract_results.py       ->  manuscript/tables/tab_*.tex   (\input by main.tex)
                                                 manuscript/figures/data/*.dat (read by TikZ)
                ->  pdflatex fig_*.tex        ->  manuscript/figures/fig_*.pdf
```

**One command (on the GB10 node)** runs every experiment into its canonical log, then extracts and plots:
```bash
./reproduce.sh all
```
Or run a single stage: `./reproduce.sh {struct|time|space|spmm|graph}` re-runs one experiment family; `./reproduce.sh extract` re-derives all tables/figure-data from the **existing** logs without recomputing; `./reproduce.sh plot` recompiles the TikZ figures.

### Canonical logs (kept under `manuscript/logs/`, gitignored)
Every log begins with a provenance header (date, host, git commit) and separates datasets with `## <key>` markers so a single parser can slice it. Log → artifact map:

| Canonical log | Feeds | Produced by |
|---|---|---|
| `manuscript/logs/geno_struct.log` | Table 1 (`tab:geno_through`) | `reproduce.sh struct` |
| `manuscript/logs/bio_results.csv` + `manuscript/logs/geno_space_energy.log` | Table 2 (`tab:geno_time`) | `reproduce.sh time`/`space` |
| `manuscript/logs/geno_space_energy.log` | Table 3 (`tab:geno`), Table B (`tab:build`, $+$pt col), Fig. 3 (`fig:space`) | `reproduce.sh space` |
| `manuscript/logs/grammar_build.log` | Table B (`tab:build`, grammar col) | `reproduce.sh grammar` |
| `manuscript/logs/geno_spmm.log` + `manuscript/logs/geno_cusparse_alg.log` | Table 4 (`tab:geno_spmm`), Fig. 4 (`fig:batched`) | `reproduce.sh spmm` |
| `manuscript/logs/graph_struct.log` | Table 5 (`tab:graph_struct`) | `reproduce.sh struct` |
| `manuscript/logs/graph_semiring.log` | Table 6 (`tab:graph`) | `reproduce.sh graph` |
| `manuscript/logs/graph_scale.log` | Table 7 (`tab:graph_scale`) | `reproduce.sh graphscale` |

### Host-side construction cost (Table B / `tab:build`)
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
- **Figure 3** (`fig:space`) and **Figure 4** (`fig:batched`) are standalone pgfplots sources — `manuscript/figures/fig_space.tex` and `fig_batched.tex` — that read `manuscript/figures/data/fig_space.dat` / `fig_batched.dat`. They compile to `manuscript/figures/fig_space.pdf` / `fig_batched.pdf`, which `main.tex` includes via `\includegraphics` (same convention as the schematic figures).
```bash
./reproduce.sh plot                 # or: cd manuscript/figures && pdflatex fig_space.tex
```
