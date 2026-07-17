#!/bin/bash -l

# transpose each row block

# sample invocation:
#   07_transpose.sh /data/mm/covtype 8
# data path must be absolute
# must be executed from the mm-repair/reordering directory 

#args
DATASET_PATH=$1
NB=$2

for (( i=0; i+1<=$NB; ++i ))
do
    python3 ./column_major.py $DATASET_PATH.$NB.$i &
    pids+=($!)
done

for pid in "${pids[@]}"; do
  wait "$pid"
done