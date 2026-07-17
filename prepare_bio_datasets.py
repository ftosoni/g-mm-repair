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

"""Prepare the genotype datasets of the manuscript (REPRODUCIBILITY.md sec. A/B):
download the 1000 Genomes Chr21/Chr20 VCFs and slice them with vcf2mat.py (100K and
full-width subsets), simulate the five synthetic haplotype matrices with
generate_msprime.py (coalescent-with-recombination, LD tuned by recombination rate,
seed fixed at 42), and write the shared {1.0, 2.0} .val value arrays. Chr22 is the
running-example chromosome and is prepared the same way in the REPRODUCIBILITY.md
commands. NOTE: paths are hardcoded relative to a PARENT directory containing the
`mm-grammar-gpu/` checkout, so run this from that parent, not from inside the repo."""
import urllib.request
import os
import subprocess
import struct

def download_file(url, dest):
    if os.path.exists(dest):
        print(f"File {dest} already exists. Skipping download.")
        return
    print(f"Downloading {url} to {dest}...")
    urllib.request.urlretrieve(url, dest)
    print("Download completed.")

def main():
    os.makedirs("mm-grammar-gpu/geno", exist_ok=True)
    os.makedirs("mm-grammar-gpu/mm-repair/data", exist_ok=True)
    
    # 1. Download Chr21 and Chr20
    vcf21_url = "https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/ALL.chr21.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz"
    vcf20_url = "https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/ALL.chr20.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz"
    
    vcf21_dest = "mm-grammar-gpu/geno/chr21.vcf.gz"
    vcf20_dest = "mm-grammar-gpu/geno/chr20.vcf.gz"
    
    download_file(vcf21_url, vcf21_dest)
    download_file(vcf20_url, vcf20_dest)
    
    # 2. Run vcf2mat.py for Chr21 and Chr20
    print("Processing VCFs to matrices...")
    vcf2mat_script = "mm-grammar-gpu/geno/vcf2mat.py"
    
    # Chr21
    subprocess.run(["python3", vcf2mat_script, vcf21_dest, "mm-grammar-gpu/mm-repair/data/geno21", "100000"], check=True)
    subprocess.run(["python3", vcf2mat_script, vcf21_dest, "mm-grammar-gpu/mm-repair/data/geno21full", "1500000"], check=True)
    
    # Chr20
    subprocess.run(["python3", vcf2mat_script, vcf20_dest, "mm-grammar-gpu/mm-repair/data/geno20", "100000"], check=True)
    subprocess.run(["python3", vcf2mat_script, vcf20_dest, "mm-grammar-gpu/mm-repair/data/geno20full", "2000000"], check=True)
    
    # 3. Generate synthetic genotype matrices under the coalescent with recombination
    #    via msprime (citable, seeded). LD is tuned by the recombination rate:
    #    low recomb (1e-9) => long shared haplotype blocks (high LD, compressible);
    #    high recomb (1e-7) => broken-up blocks (low LD). Seed fixed for reproducibility.
    synth_script = "mm-grammar-gpu/generate_msprime.py"
    print("Generating synthetic datasets via msprime...")
    SEED = "42"
    #                                    rows     cols      out                                                  recomb  seed
    subprocess.run(["python3", synth_script, "2000",  "50000",  "mm-grammar-gpu/mm-repair/data/geno_synth_small",     "1e-8", SEED], check=True)
    subprocess.run(["python3", synth_script, "5000",  "200000", "mm-grammar-gpu/mm-repair/data/geno_synth_large",     "1e-8", SEED], check=True)
    subprocess.run(["python3", synth_script, "5000",  "100000", "mm-grammar-gpu/mm-repair/data/geno_synth_ld_high",   "1e-9", SEED], check=True)
    subprocess.run(["python3", synth_script, "5000",  "100000", "mm-grammar-gpu/mm-repair/data/geno_synth_ld_low",    "1e-7", SEED], check=True)
    subprocess.run(["python3", synth_script, "10000", "50000",  "mm-grammar-gpu/mm-repair/data/geno_synth_ind_large", "1e-8", SEED], check=True)

    # 4. Generate .val files
    print("Writing .val files...")
    geno_val_datasets = ["geno21", "geno21full", "geno20", "geno20full",
                         "geno_synth_small", "geno_synth_large",
                         "geno_synth_ld_high", "geno_synth_ld_low", "geno_synth_ind_large"]
    for ds in geno_val_datasets:
        val_path = f"mm-grammar-gpu/mm-repair/data/{ds}.val"
        with open(val_path, "wb") as f:
            f.write(struct.pack("dd", 1.0, 2.0))
            
    print("Dataset preparation complete.")

if __name__ == "__main__":
    main()
