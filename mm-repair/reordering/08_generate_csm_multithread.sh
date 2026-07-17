#!/bin/bash -l

# generate the locally-pruned CMS (multithreaded version)

# sample invocation:
#   08_generate_csm_multithread.sh /data/mm/covtype 581012 54 16 4 1024 16
# data path must be absolute
# must be executed from the mm-repair/reordering directory 

#args
DATASET_PATH=$1
NR=$2
NC=$3
K_PAR=$4
PARDEGREE=$5
MEMLIM=$6
NB=$7

#params
let NR_FIRST=($NR+$NB-1)/$NB
let NR_LAST=($NB-1)*$NR_FIRST
let NR_LAST=$NR-$NR_LAST
let NB_LAST=$NB-1


set -eu
pids=()

#first blocks
for (( i=0; i+1<$NB; ++i ))
do
    ./build/tsp_generator_pruned_local_mt $DATASET_PATH.$NB.$i $NR_FIRST $NC $K_PAR $PARDEGREE $MEMLIM  &
    pids+=($!)
done
#last block
./build/tsp_generator_pruned_local_mt $DATASET_PATH.$NB.$NB_LAST $NR_LAST $NC $K_PAR $PARDEGREE $MEMLIM  &
pids+=($!)

for pid in "${pids[@]}"; do
  wait "$pid"
done