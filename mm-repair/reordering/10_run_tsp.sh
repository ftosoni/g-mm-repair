#!/bin/bash -l

# reorder each row block via LKH

# sample invocation:
#   10_run_tsp.sh /data/mm/covtype 8 16
# data path must be absolute
# must be executed from the mm-repair/reordering directory 

#args
DATASET_PATH=$1
NB=$2
K_PAR=$3
LKH_VERSION=$4

set -eu
pids=()

for (( i=0; i<$NB; ++i ))
do
    $LKH_VERSION/LKH $DATASET_PATH.$NB.$i.pruned_local_$K_PAR.par > $DATASET_PATH.$NB.$i.pruned_local_$K_PAR.log &
    pids+=($!)
done

for pid in "${pids[@]}"; do
  wait "$pid"
done