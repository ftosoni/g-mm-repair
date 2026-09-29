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

"""Reproduce Table 4.3 (tab:geno): run the grammar engine and the cuSPARSE CSR
baseline across the genotype matrices plus the crossover/OOM matrix, collecting
time, analytic device bytes, measured peak bytes, and energy/vector. Produces the
combined space/energy table for the paper (see REPRODUCIBILITY.md section on
Table 4.3). crossover_synth is where cuSPARSE OOMs while the engine stays resident."""
import subprocess, os, re

# Genotype matrices + the crossover/OOM matrix that make up Table 4.3 (tab:geno) of the
# manuscript. (path, rows, cols); paths mirror run_all_bio_baselines.sh.
datasets = [
    ("geno/geno22",                        2504,  100000),
    ("geno/geno22full",                    2504,  1055454),
    ("mm-repair/data/geno21",              2504,  100000),
    ("mm-repair/data/geno21full",          2504,  1054447),
    ("mm-repair/data/geno20",              2504,  100000),
    ("mm-repair/data/geno20full",          2504,  1739315),
    ("mm-repair/data/geno_synth_small",    2000,  50000),
    ("mm-repair/data/geno_synth_large",    5000,  200000),
    ("mm-repair/data/geno_synth_ld_high",  5000,  100000),
    ("mm-repair/data/geno_synth_ld_low",   5000,  100000),
    ("mm-repair/data/geno_synth_ind_large",10000, 50000),
    ("mm-repair/data/crossover_synth",     10000, 900000),
]

env = os.environ.copy()
env["PATH"] = f"/usr/local/cuda/bin:{env.get('PATH','')}"

def parse(out, time_re):
    g = lambda rx: (float(m.group(1)) if (m := re.search(rx, out)) else None)
    return {
        "ms":       g(time_re),
        "omp":      g(r"CPU reference \(OpenMP\): ([\d.]+) ms/vector"),
        "seq":      g(r"CPU reference \(sequential\): ([\d.]+) ms/vector"),
        "analytic": g(r"MEM analytic bytes: (\d+)"),
        "peak":     g(r"MEM peak bytes: (\d+)"),
        "mj":       g(r"ENERGY mJ/vector: ([\d.eE+-]+)"),
        "afull":    g(r"alpha (\d+) ->"),
        "anew":     g(r"-> (\d+) \("),
        "ok":       "SUCCESS" in out,
    }

def run(cmd):
    try:
        r = subprocess.run(cmd, capture_output=True, text=True, env=env, timeout=1800)
        return r.stdout + "\n" + r.stderr
    except Exception as e:
        return f"ERROR {e}"

rows = []
for path, nr, nc in datasets:
    name = os.path.basename(path)
    print(f"== {name} ==", flush=True)
    eng = parse(run(["./gpu-engine/gpu_test", path, str(nr), str(nc), "100"]),
                r"GPU kernel: .* = ([\d.]+) ms/vector")
    cus = parse(run(["./gpu-engine/cusparse_test", path, str(nr), str(nc), "100"]),
                r"Average time: ([\d.]+) ms/vector")
    rows.append((name, eng, cus))
    print(f"   engine: {eng}\n   cusparse: {cus}", flush=True)

# Results/logs go under manuscript/logs/ (gitignored), co-located with the tables/figures.
os.makedirs("manuscript/logs", exist_ok=True)

# raw dump for transcription into the paper
with open("manuscript/logs/space_energy_raw.txt", "w") as fh:
    for name, e, c in rows:
        fh.write(f"{name}\n  engine={e}\n  cusparse={c}\n")

def mb(b):  return f"{b/1e6:.1f}" if b else "N/A"
def f(x, p=2): return f"{x:.{p}f}" if x is not None else "N/A"

hdr = (f"{'matrix':<11}| "
       f"{'eng ms':>8} {'eng MB':>8} {'eng pkMB':>9} {'eng mJ':>8} | "
       f"{'cus ms':>8} {'cus MB':>8} {'cus pkMB':>9} {'cus mJ':>8}")
lines = ["SPACE / ENERGY: grammar engine vs cuSPARSE, unified memory, GB10",
         "=" * len(hdr), hdr, "-" * len(hdr)]
for name, e, c in rows:
    lines.append(f"{name:<11}| "
                 f"{f(e['ms']):>8} {mb(e['analytic']):>8} {mb(e['peak']):>9} {f(e['mj']):>8} | "
                 f"{f(c['ms']):>8} {mb(c['analytic']):>8} {mb(c['peak']):>9} {f(c['mj']):>8}")
lines.append("=" * len(hdr))
out = "\n".join(lines)
print(out)
with open("manuscript/logs/space_energy_results.txt", "w") as fh:
    fh.write(out + "\n")
