#!/bin/bash

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

cd "$(dirname "$0")"

export PATH=/usr/local/cuda/bin:$PATH

# Results go under manuscript/logs/ (gitignored), co-located with the tables/figures.
mkdir -p manuscript/logs
RESULTS=manuscript/logs/bio_results.csv
echo "Dataset,GPU,OMP,CPU_Seq,mm_seq,mm_16" > "$RESULTS"

run_ds() {
    name=$1
    rows=$2
    cols=$3
    base=$4
    
    echo "Running benchmark for $name..."
    
    # Run single-vector benchmarks
    out_single=$(python3 run_benchmark.py $base $rows $cols 100)
    
    # Parse times
    gpu=$(echo "$out_single" | grep -E "Level-Synchronous GPU Sweep" | grep -o -E "[0-9.]+ ms" | cut -d' ' -f1)
    omp=$(echo "$out_single" | grep -E "Level-Synchronous CPU Sweep \(OpenMP" | grep -o -E "[0-9.]+ ms" | cut -d' ' -f1)
    cpu_seq=$(echo "$out_single" | grep -E "Level-Synchronous CPU Sweep \(Sequential" | grep -o -E "[0-9.]+ ms" | cut -d' ' -f1)
    mm_seq=$(echo "$out_single" | grep -E "mm-repair \(re32mm" | grep -o -E "[0-9.]+ ms" | cut -d' ' -f1)
    
    # Run 16-thread benchmark
    out_16=$(python3 run_benchmark_16.py $base $rows $cols 100 --i32)
    mm_16=$(echo "$out_16" | grep -E "mm-repair \(re32mm" | grep -o -E "[0-9.]+ ms" | cut -d' ' -f1)
    
    echo "$name,$gpu,$omp,$cpu_seq,$mm_seq,$mm_16" >> "$RESULTS"
}

run_ds "geno22" 2504 100000 "geno/geno22"
run_ds "geno22full" 2504 1055454 "geno/geno22full"
run_ds "geno21" 2504 100000 "mm-repair/data/geno21"
run_ds "geno21full" 2504 1054447 "mm-repair/data/geno21full"
run_ds "geno20" 2504 100000 "mm-repair/data/geno20"
run_ds "geno20full" 2504 1739315 "mm-repair/data/geno20full"
run_ds "geno_synth_small" 2000 50000 "mm-repair/data/geno_synth_small"
run_ds "geno_synth_large" 5000 200000 "mm-repair/data/geno_synth_large"
run_ds "geno_synth_ld_high" 5000 100000 "mm-repair/data/geno_synth_ld_high"
run_ds "geno_synth_ld_low" 5000 100000 "mm-repair/data/geno_synth_ld_low"
run_ds "geno_synth_ind_large" 10000 50000 "mm-repair/data/geno_synth_ind_large"

cat "$RESULTS"
