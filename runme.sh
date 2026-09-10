#!/bin/bash
# This file is part of g-mm-repair <https://github.com/ftosoni/g-mm-repair>.
# Copyright (c) 2026 Francesco Tosoni. Apache-2.0.
#
# High-level one-shot artifact driver (SIAM/ALENEX Artifact Evaluation).
# It only orchestrates steps that already exist and are documented in README.md /
# REPRODUCIBILITY.md -- nothing new:
#
#   1. initialise the mm-repair submodule
#   2. build the GPU engine + cuSPARSE baseline and the mm-repair toolchain
#   3. fetch + verify the Zenodo data package into ./zenodo/ (skipped if present)
#   4. run ./reproduce.sh all  -> experiments -> tables (manuscript/tables/*.tex)
#      + figure data (manuscript/figures/data/*.dat) + compiled figures
#
# Before running on a non-GB10 GPU, set the arch in gpu-engine/Makefile (NVCCFLAGS
# -arch=sm_XX) to your compute capability, and make sure SDSL-lite is installed.
#
# Usage:   ./runme.sh
# Env:     SKIP_DATA=1 ./runme.sh     # do not touch ./zenodo/ (data already in place)
#          ZENODO_DIR=/path ./runme.sh
set -euo pipefail
cd "$(dirname "$0")"
export PATH=/usr/local/cuda/bin:$PATH
ZENODO_DIR=${ZENODO_DIR:-zenodo}
say() { printf '\n\033[1m==== %s ====\033[0m\n' "$*"; }

say "1/4  Initialise submodule (mm-repair)"
git submodule update --init --recursive

say "2/4  Build (GPU engine + cuSPARSE baseline, then mm-repair)"
echo "    NOTE: gpu-engine/Makefile targets -arch=sm_121 (GB10). Edit NVCCFLAGS for your GPU"
echo "          (e.g. sm_90 Hopper, sm_89 Ada, sm_80 Ampere) or the binaries will not launch."
echo "    NOTE: mm-repair needs SDSL-lite installed (see mm-repair/Readme.md)."
make -C gpu-engine all
make -C mm-repair all

say "3/4  Data package (Zenodo doi:10.5281/zenodo.22677746)"
if [ "${SKIP_DATA:-0}" = "1" ]; then
  echo "    SKIP_DATA=1 -> assuming datasets already in $ZENODO_DIR/"
elif [ -d "$ZENODO_DIR/genotypes" ] && [ -d "$ZENODO_DIR/wikidata" ] && [ -d "$ZENODO_DIR/swh" ]; then
  echo "    $ZENODO_DIR/ already populated -> skipping download."
else
  echo "    Downloading + extracting + verifying into $ZENODO_DIR/ ..."
  mkdir -p "$ZENODO_DIR"; ( cd "$ZENODO_DIR"
    REC=https://zenodo.org/api/records/22677747/files
    for f in genotypes.tar wikidata.tar swh.tar README.txt MANIFEST.md5; do
      curl -L --retry 5 --retry-all-errors -o "$f" "$REC/$f/content"
    done
    for t in genotypes.tar wikidata.tar swh.tar; do tar xf "$t"; done
    md5sum -c MANIFEST.md5 )
fi

say "4/4  Reproduce (experiments -> extract tables/figure-data -> plot figures)"
ZENODO_DIR="$ZENODO_DIR" ./reproduce.sh all

say "DONE"
echo "  Tables : manuscript/tables/tab_*.tex"
echo "  Figures: manuscript/figures/data/*.dat + manuscript/figures/fig_*.pdf"
echo "  Heavy, on-demand (not in 'all'): ./reproduce.sh grammar | graphscale | crosscheck"
