#!/bin/bash -l

# reorder each row block via PathCover+

# sample invocation:
#   10_run_cover2.sh /data/mm/covtype 8 16
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
    python3 ./cover2.py $DATASET_PATH.$NB.$i.pruned_local_$K_PAR.tsp &
    pids+=($!)
done

for pid in "${pids[@]}"; do
  wait "$pid"
done