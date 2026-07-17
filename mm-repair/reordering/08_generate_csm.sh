
#!/bin/bash -l

# generate the CSM for each row block

# sample invocation:
#   08_generate_csm.sh /data/mm/covtype 581012 54 8 16
# data path must be absolute
# must be executed from the mm-repair/reordering directory 

#args
DATASET_PATH=$1
NR=$2
NC=$3
NB=$4
K_PAR=$5

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
    ./build/tsp_generator_pruned_local $DATASET_PATH.$NB.$i $NR_FIRST $NC $K_PAR &
    pids+=($!)
done
#last block
./build/tsp_generator_pruned_local $DATASET_PATH.$NB.$NB_LAST $NR_LAST $NC $K_PAR &
pids+=($!)

for pid in "${pids[@]}"; do
  wait "$pid"
done