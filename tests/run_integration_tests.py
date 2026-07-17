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

import os
import sys
import subprocess
import random
import struct
import shutil

def run_cmd(cmd, cwd=None):
    print(f"Running: {' '.join(cmd)}")
    res = subprocess.run(cmd, cwd=cwd, capture_output=True, text=True)
    if res.returncode != 0:
        print(f"Command failed with exit code {res.returncode}")
        print("Stdout:")
        print(res.stdout)
        print("Stderr:")
        print(res.stderr)
        sys.exit(res.returncode)
    return res

def main():
    base_dir = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    gpu_engine_dir = os.path.join(base_dir, "gpu-engine")
    mm_repair_dir = os.path.join(base_dir, "mm-repair")
    
    print("=== [1/4] Compiling CPU Test Suite and Tools ===")
    
    # 1. Compile cpu_selftest
    print("Compiling cpu_selftest...")
    cxx = os.environ.get("CXX", "g++")
    cpu_selftest_cmd = [
        cxx, "-std=c++17", "-O3",
        "-o", os.path.join(gpu_engine_dir, "cpu_selftest"),
        os.path.join(gpu_engine_dir, "cpu_selftest.cpp"),
        os.path.join(gpu_engine_dir, "build_schedule.cpp"),
        "-fopenmp"
    ]
    # On Windows, we append .exe suffix to target output
    if sys.platform == "win32":
        cpu_selftest_cmd[4] += ".exe"
    run_cmd(cpu_selftest_cmd)

    # 2. Compile irepair0
    print("Compiling irepair0...")
    cc = os.environ.get("CC", "gcc")
    irepair_src = ["irepair0.c", "array.c", "hash.c", "heap.c", "records.c", "basics.c"]
    irepair_cmd = [
        cc, "-O3", "-o", os.path.join(mm_repair_dir, "brepair", "irepair0"),
    ] + [os.path.join(mm_repair_dir, "brepair", src) for src in irepair_src]
    if sys.platform == "win32":
        irepair_cmd[3] += ".exe"
    run_cmd(irepair_cmd)

    # 3. Compile bin2csrv
    print("Compiling bin2csrv...")
    bin2csrv_cmd = [
        cxx, "-std=c++17", "-O3",
        "-o", os.path.join(mm_repair_dir, "bin2csrv"),
        os.path.join(mm_repair_dir, "bin2csrv.cpp")
    ]
    if sys.platform == "win32":
        bin2csrv_cmd[4] += ".exe"
    run_cmd(bin2csrv_cmd)

    # 3.5. Compile re32mm CPU reference
    print("Compiling re32mm...")
    re32mm_cmd = [
        cc, "-Wall", "-std=c99", "-O3",
        "-o", os.path.join(mm_repair_dir, "re32mm"),
        os.path.join(mm_repair_dir, "remm.c"),
        "-pthread"
    ]
    if sys.platform == "win32":
        re32mm_cmd[5] += ".exe"
    run_cmd(re32mm_cmd)

    # 4. Compile dummy encoders to satisfy matrepair requirements without sdsl/ans libraries
    print("Compiling dummy SDSL and ANS encoders...")
    dummy_c_path = os.path.join(base_dir, "tests", "dummy.c")
    os.makedirs(os.path.dirname(dummy_c_path), exist_ok=True)
    with open(dummy_c_path, "w") as f:
        f.write("int main() { return 0; }\n")
        
    os.makedirs(os.path.join(mm_repair_dir, "sdsl"), exist_ok=True)
    os.makedirs(os.path.join(mm_repair_dir, "ans"), exist_ok=True)
    
    sdsl_dummy_cmd = [cc, "-O3", "-o", os.path.join(mm_repair_dir, "sdsl", "encode.x"), dummy_c_path]
    ans_dummy_cmd = [cc, "-O3", "-o", os.path.join(mm_repair_dir, "ans", "encode.x"), dummy_c_path]
    if sys.platform == "win32":
        sdsl_dummy_cmd[3] += ".exe"
        ans_dummy_cmd[3] += ".exe"
    run_cmd(sdsl_dummy_cmd)
    run_cmd(ans_dummy_cmd)

    print("\n=== [2/4] Generating Small Synthetic Dataset ===")
    rows = 50
    cols = 50
    
    # Generate repeated patterns to make it easily compressible
    patterns = []
    random.seed(42)
    for _ in range(5):
        patterns.append([float(random.choice([0, 1, 2])) for _ in range(cols)])
        
    matrix = []
    for _ in range(rows):
        matrix.extend(random.choice(patterns))
        
    matrix_file = os.path.join(base_dir, "tests", "small_matrix.bin")
    with open(matrix_file, "wb") as f:
        for val in matrix:
            f.write(struct.pack("d", val))
    print(f"Generated synthetic matrix and saved to {matrix_file}")

    print("\n=== [3/4] Compressing with RePair ===")
    matrepair_py = os.path.join(mm_repair_dir, "matrepair")
    # run matrepair
    matrepair_cmd = [
        sys.executable, matrepair_py,
        matrix_file, str(rows), str(cols),
        "--f64"
    ]
    run_cmd(matrepair_cmd)
    print("RePair compression completed successfully.")

    print("\n=== [3.5/5] Generating Input Vector and Running re32mm CPU Reference ===")
    # Generate random input vector x as double
    x_vector = [random.uniform(-1.0, 1.0) for _ in range(cols)]
    x_file = os.path.join(base_dir, "tests", "xvector.bin")
    with open(x_file, "wb") as f:
        for val in x_vector:
            f.write(struct.pack("d", val))
            
    # Run re32mm
    re32mm_exe = os.path.join(mm_repair_dir, "re32mm")
    if sys.platform == "win32":
        re32mm_exe += ".exe"
    y_ref_file = os.path.join(base_dir, "tests", "yvector_re32mm.bin")
    
    run_cmd([re32mm_exe, "-y", y_ref_file, matrix_file, str(rows), str(cols), x_file])
    print("re32mm CPU reference multiplication completed successfully.")

    print("\n=== [4/4] Running Unified CPU Self-Test (RePair construction) ===")
    cpu_selftest_exe = os.path.join(gpu_engine_dir, "cpu_selftest")
    if sys.platform == "win32":
        cpu_selftest_exe += ".exe"
        
    # First, run against the brute-force baseline
    run_cmd([cpu_selftest_exe, matrix_file, str(rows), str(cols)])
    # Second, run against the re32mm CPU reference baseline
    run_cmd([cpu_selftest_exe, matrix_file, str(rows), str(cols), x_file, y_ref_file])
    # Third, run programmatic unit test covering all 4 kernel cases
    run_cmd([cpu_selftest_exe, "--unittest"])

    print("\n=== Cleanup ===")
    # Remove generated matrix and compression outputs
    extensions = [
        ".bin", ".bin.val", ".bin.vc", ".bin.vc.R", ".bin.vc.C",
        ".bin.vc.C.iv", ".bin.vc.R.iv", ".bin.vc.C.ansf.1", ".bin.log"
    ]
    for ext in extensions:
        path = os.path.join(base_dir, "tests", "small_matrix" + ext)
        if os.path.exists(path):
            os.remove(path)
            
    if os.path.exists(dummy_c_path):
        os.remove(dummy_c_path)
        
    if os.path.exists(x_file):
        os.remove(x_file)
        
    if os.path.exists(y_ref_file):
        os.remove(y_ref_file)
        
    print("\nALL INTEGRATION TESTS PASSED SUCCESSFULLY!")

if __name__ == "__main__":
    main()
