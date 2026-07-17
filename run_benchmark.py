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

"""Single-vector benchmark for one matrix: compiles the GPU engine + mm-repair's
re32mm, runs y = M x, and prints the GPU / OpenMP / CPU-sequential / mm-repair
(re32mm) times side by side. This is the per-dataset driver that
run_all_bio_baselines.sh loops over to fill Table 2 (tab:geno_time); the summary
labels (1)-(4) map to that table's columns. The GPU/OMP/seq times are parsed from
gpu_test's stdout, the mm-repair time from re32mm's DETAILED_TIMING stderr."""
import os
import sys
import subprocess
import struct
import random
import re

def print_help():
    print("Usage: python3 run_benchmark.py <matrix_base_path> <rows> <cols> [iters]")
    print("Example: python3 run_benchmark.py mm-repair/data/geno22 2504 100000 100")

def main():
    if len(sys.argv) < 4:
        print_help()
        sys.exit(1)

    matrix_base = sys.argv[1]
    rows = int(sys.argv[2])
    cols = int(sys.argv[3])
    iters = int(sys.argv[4]) if len(sys.argv) >= 5 else 100

    print("=" * 70)
    print(f"Grammar-Compressed Matrix Multiplication Benchmark")
    print(f"Matrix Base: {matrix_base}")
    print(f"Dimensions:  {rows} x {cols}")
    print(f"Iterations:  {iters}")
    print("=" * 70)

    # 1. Compile targets
    print("\n[1/3] Compiling executables...")
    
    # Compile GPU engine
    print("-> Compiling gpu-engine...")
    try:
        subprocess.run(["make", "clean"], cwd="gpu-engine", check=True, stdout=subprocess.DEVNULL)
        subprocess.run(["make"], cwd="gpu-engine", check=True, stdout=subprocess.DEVNULL)
    except subprocess.CalledProcessError as e:
        print(f"Error compiling gpu-engine: {e}")
        sys.exit(1)

    # Compile mm-repair re32mm with DETAILED_TIMING
    print("-> Compiling mm-repair (re32mm with DETAILED_TIMING)...")
    try:
        subprocess.run(["make", "clean"], cwd="mm-repair", check=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
        # Pass CFLAGS as a make argument: a command-line assignment overrides the
        # makefile's own `CFLAGS=` (an environment CFLAGS would be ignored).
        subprocess.run(["make", "re32mm", "CFLAGS=-Wall -std=c99 -g -O3 -DDETAILED_TIMING"],
                       cwd="mm-repair", check=True, stdout=subprocess.DEVNULL)
    except subprocess.CalledProcessError as e:
        print(f"Error compiling mm-repair: {e}")
        sys.exit(1)

    # 2. Generate random input vector
    print("\n[2/3] Generating random input vector...")
    vector_path = "temp_x.bin"
    # re32mm expects double (float64) values in binary format
    with open(vector_path, "wb") as f:
        for _ in range(cols):
            val = random.random()
            f.write(struct.pack("d", val))

    # 3. Run benchmarks
    print("\n[3/3] Running benchmarks...")
    
    # GPU / CPU Sweep benchmark
    gpu_cmd = ["./gpu-engine/gpu_test", matrix_base, str(rows), str(cols), str(iters)]
    print(f"-> Running: {' '.join(gpu_cmd)}")
    gpu_output = ""
    try:
        res = subprocess.run(gpu_cmd, capture_output=True, text=True, check=True)
        gpu_output = res.stdout
    except subprocess.CalledProcessError as e:
        print(f"Error running GPU benchmark: {e}")
        print(e.stdout)
        print(e.stderr)
        if os.path.exists(vector_path):
            os.remove(vector_path)
        sys.exit(1)

    # mm-repair benchmark
    mm_cmd = ["./mm-repair/re32mm", "-n", str(iters), "-b", "1", matrix_base, str(rows), str(cols), vector_path]
    print(f"-> Running: {' '.join(mm_cmd)}")
    mm_output = ""
    try:
        res = subprocess.run(mm_cmd, capture_output=True, text=True, check=True)
        mm_output = res.stderr # DETAILED_TIMING output goes to stderr
    except subprocess.CalledProcessError as e:
        print(f"Error running mm-repair benchmark: {e}")
        print(e.stdout)
        print(e.stderr)
        if os.path.exists(vector_path):
            os.remove(vector_path)
        sys.exit(1)

    # Clean up temporary vector
    if os.path.exists(vector_path):
        os.remove(vector_path)

    # 4. Parse results
    # GPU Time
    # e.g.: GPU kernel: 202.16 ms total over 100 vectors = 2.0216 ms/vector.
    gpu_match = re.search(r"GPU kernel: .* = ([\d.]+) ms/vector", gpu_output)
    gpu_ms = float(gpu_match.group(1)) if gpu_match else None

    # CPU OMP Time
    # e.g.: CPU reference (OpenMP): 35.12 ms/vector (64 threads).
    omp_match = re.search(r"CPU reference \(OpenMP\): ([\d.]+) ms/vector", gpu_output)
    omp_ms = float(omp_match.group(1)) if omp_match else None

    # CPU Sequential Time (our implementation)
    # e.g.: CPU reference (sequential): 120.45 ms/vector.
    seq_match = re.search(r"CPU reference \(sequential\): ([\d.]+) ms/vector", gpu_output)
    seq_ms = float(seq_match.group(1)) if seq_match else None

    # mm-repair Sequential Time
    # e.g.: Average mult time (secs) Ax: 0.125432  xA: 0.118932
    mm_match = re.search(r"Average mult time \(secs\) Ax:\s+([\d.]+)", mm_output)
    mm_ms = float(mm_match.group(1)) * 1000.0 if mm_match else None

    # Print Summary Table
    print("\n" + "=" * 70)
    print("                    BENCHMARK RESULTS SUMMARY                    ")
    print("=" * 70)
    print(f"{'Method / Implementation':<45} | {'Avg Time / Vector':<20}")
    print("-" * 70)
    
    if gpu_ms is not None:
        print(f"{'(1) Level-Synchronous GPU Sweep (Emit-on-the-spot)':<45} | {gpu_ms:>13.4f} ms")
    else:
        print(f"{'(1) Level-Synchronous GPU Sweep (Emit-on-the-spot)':<45} | {'N/A':>16}")

    if omp_ms is not None:
        print(f"{'(3) Level-Synchronous CPU Sweep (OpenMP parallel)':<45} | {omp_ms:>13.4f} ms")
    else:
        print(f"{'(3) Level-Synchronous CPU Sweep (OpenMP parallel)':<45} | {'N/A':>16}")

    if seq_ms is not None:
        print(f"{'(2) Level-Synchronous CPU Sweep (Sequential)':<45} | {seq_ms:>13.4f} ms")
    else:
        print(f"{'(2) Level-Synchronous CPU Sweep (Sequential)':<45} | {'N/A':>16}")

    if mm_ms is not None:
        print(f"{'(4) mm-repair (re32mm, CPU Sequential)':<45} | {mm_ms:>13.4f} ms")
    else:
        print(f"{'(4) mm-repair (re32mm, CPU Sequential)':<45} | {'N/A':>16}")

    print("=" * 70)

    # Calculate Speedups
    if mm_ms and gpu_ms and gpu_ms > 0:
        print(f"Speedup GPU vs mm-repair (re32mm): {mm_ms / gpu_ms:.2f}x")
    if seq_ms and gpu_ms and gpu_ms > 0:
        print(f"Speedup GPU vs CPU Sequential:     {seq_ms / gpu_ms:.2f}x")
    if seq_ms and omp_ms and omp_ms > 0:
        print(f"Speedup OpenMP vs CPU Sequential:  {seq_ms / omp_ms:.2f}x")
    print("=" * 70)

if __name__ == "__main__":
    main()