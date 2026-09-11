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

"""16-thread mm-repair benchmark for one matrix: rebuilds the RePair grammar with
16 row-blocks (-b 16 -p 16) and times re32mm at 16 threads. Produces the
'mmr (16th)' column of Table 4.2 (tab:geno_time); run_all_bio_baselines.sh calls it
alongside run_benchmark.py, which supplies the other columns. mm-repair parallelizes
coarsely (one independent grammar per row-block), which is the baseline the paper
contrasts with its fine-grained, inside-one-grammar GPU scheme."""
import os
import sys
import subprocess
import struct
import random
import re

def print_help():
    print("Usage: python3 run_benchmark_16.py <matrix_base_path> <rows> <cols> [iters] [format_option]")
    print("Format options: --i32, --f32, --f64 (default: --f64)")
    print("Example: python3 run_benchmark_16.py mm-repair/data/geno22 2504 100000 100 --i32")

def main():
    if len(sys.argv) < 4:
        print_help()
        sys.exit(1)

    matrix_base = sys.argv[1]
    rows = int(sys.argv[2])
    cols = int(sys.argv[3])
    iters = int(sys.argv[4]) if len(sys.argv) >= 5 and not sys.argv[4].startswith("-") else 100
    
    # Parse format option
    format_opt = "--f64"
    for arg in sys.argv[4:]:
        if arg in ["--i32", "--f32", "--f64"]:
            format_opt = arg
            break

    print("=" * 70)
    print(f"Grammar-Compressed Matrix Multithreaded Benchmark (16 Threads)")
    print(f"Matrix Base: {matrix_base}")
    print(f"Dimensions:  {rows} x {cols}")
    print(f"Iterations:  {iters}")
    print(f"Format:      {format_opt}")
    print("=" * 70)

    # 1. Compile targets
    print("\n[1/4] Compiling executables...")
    
    # Compile mm-repair targets with DETAILED_TIMING
    print("-> Compiling mm-repair targets with DETAILED_TIMING)...")
    try:
        subprocess.run(["make", "clean"], cwd="mm-repair", check=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
        subprocess.run(["make", "CFLAGS=-Wall -std=c99 -g -O3 -DDETAILED_TIMING"],
                       cwd="mm-repair", check=True, stdout=subprocess.DEVNULL)
    except subprocess.CalledProcessError as e:
        print(f"Error compiling mm-repair: {e}")
        sys.exit(1)

    # 2. Build the grammar with 16 blocks using matrepair
    # The 16-block grammar (mm_16 column) is rebuilt from the bare int32 dense, which
    # matrepair stats as its source even in lazy (-y) mode. The Zenodo package ships the
    # dense zstd-compressed (<base>.zst), so inflate it on demand if the bare file is
    # absent -- keeps run_all_bio_baselines.sh working out-of-the-box on the package.
    if not os.path.exists(matrix_base) and os.path.exists(matrix_base + ".zst"):
        print(f"-> Bare dense '{matrix_base}' missing; decompressing {matrix_base}.zst ...")
        try:
            subprocess.run(["zstd", "-d", "-k", "-f", matrix_base + ".zst"], check=True)
        except (subprocess.CalledProcessError, FileNotFoundError) as e:
            print(f"Error decompressing {matrix_base}.zst: {e}")
            sys.exit(1)

    print("\n[2/4] Constructing RePair grammar with 16 row-blocks (16 threads)...")
    matrepair_cmd = [
        "python3", "matrepair", "-r", "-b", "16", "-p", "16",
        format_opt, os.path.abspath(matrix_base), str(rows), str(cols)
    ]
    print(f"-> Running: {' '.join(matrepair_cmd)}")
    try:
        # Run inside the mm-repair directory where matrepair resides
        subprocess.run(matrepair_cmd, cwd="mm-repair", check=True)
    except subprocess.CalledProcessError as e:
        print(f"Error running matrepair: {e}")
        sys.exit(1)

    # 3. Generate random input vector
    print("\n[3/4] Generating random input vector...")
    vector_path = "temp_x.bin"
    # re32mm expects double (float64) values in binary format
    with open(vector_path, "wb") as f:
        for _ in range(cols):
            val = random.random()
            f.write(struct.pack("d", val))

    # 4. Run benchmarks
    print("\n[4/4] Running benchmark...")
    
    # mm-repair benchmark at 16 threads
    mm_cmd = ["./mm-repair/re32mm", "-n", str(iters), "-b", "16", matrix_base, str(rows), str(cols), vector_path]
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

    # 5. Parse results
    # mm-repair Sequential/Parallel Time
    # e.g.: Average mult time (secs) Ax: 0.125432  xA: 0.118932
    mm_match = re.search(r"Average mult time \(secs\) Ax:\s+([\d.]+)", mm_output)
    mm_ms = float(mm_match.group(1)) * 1000.0 if mm_match else None

    # Print Summary Table
    print("\n" + "=" * 70)
    print("                    BENCHMARK RESULTS SUMMARY                    ")
    print("=" * 70)
    print(f"{'Method / Implementation':<45} | {'Avg Time / Vector':<20}")
    print("-" * 70)
    
    if mm_ms is not None:
        print(f"{'mm-repair (re32mm, CPU 16 Threads)':<45} | {mm_ms:>13.4f} ms")
    else:
        print(f"{'mm-repair (re32mm, CPU 16 Threads)':<45} | {'N/A':>16}")

    print("=" * 70)

if __name__ == "__main__":
    main()