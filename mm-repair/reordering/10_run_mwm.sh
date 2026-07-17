#!/bin/bash -l

# reorder each row block via MWM

# sample invocation:
#   10_run_mwm.sh /data/mm/covtype 8 16
# data path must be absolute
# must be executed from the mm-repair/reordering directory 

#args
DATASET_PATH=$1
NB=$2
K_PAR=$3

set -eu
pids=()

for (( i=0; i<$NB; ++i ))
do
    ./build/mwm $DATASET_PATH.$NB.$i pruned_local_$K_PAR &&
    python3 ./mwm_sol_from_pairs.py $DATASET_PATH.$NB.$i pruned_local_$K_PAR & 
    pids+=($!)
done

for pid in "${pids[@]}"; do
  wait "$pid"
done