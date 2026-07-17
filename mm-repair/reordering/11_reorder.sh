#!/bin/bash -l

# reorder each row block

# sample invocation:
#   11_reorder.sh /data/mm/covtype 581012 54 8 16 tsp
# data path must be absolute
# must be executed from the mm-repair/reordering directory 

#args
DATASET_PATH=$1
NR=$2
NC=$3
NB=$4
K_PAR=$5
ALGO=$6

#params
let NR_FIRST=($NR+$NB-1)/$NB
let NR_LAST=($NB-1)*$NR_FIRST
let NR_LAST=$NR-$NR_LAST
let NB_LAST=$NB-1

set -eu
pids=()

# .vco are the original unpermuted .vc files 
for (( i=0; i+1<$NB; ++i ))
do
    ./vc_reorder.x $DATASET_PATH.$NB.$i $NR_FIRST $NC $DATASET_PATH.$NB.$i.pruned_local_$K_PAR.$ALGO.solution && 
    echo " --- block $i done." &
    pids+=($!)
done
./vc_reorder.x $DATASET_PATH.$NB.$NB_LAST $NR_LAST $NC $DATASET_PATH.$NB.$NB_LAST.pruned_local_$K_PAR.$ALGO.solution &&
echo " --- block $NB_LAST done." &
pids+=($!)

for pid in "${pids[@]}"; do
  wait "$pid"
done
