#!/bin/bash
cd /mnt/nfs/home/tosoni/mm-grammar-gpu
SP=cgvenv/lib/python3.12/site-packages
export LD_LIBRARY_PATH="$(find $SP -type d \( -name lib -o -name lib64 \) | tr "\n" :):/usr/local/cuda/lib64"
export PATH=/usr/local/cuda/bin:$PATH
declare -A DIM=( [yago_plays_for]="91397 13765" [yago_is_citizen_of]="45823 309" [yago_acted_in]="26665 39559" [wd_sports_team]="332121 29854" [wd_cast_member]="173977 144095" [wd_citizenship]="2874250 2556" [wd_occupation]="3459933 10610" [wd_subclass_of]="1487709 73417" )
for b in yago_plays_for yago_is_citizen_of yago_acted_in wd_sports_team wd_cast_member wd_citizenship wd_occupation wd_subclass_of; do
  echo "##### $b #####"
  timeout 300 ./cgvenv/bin/python cugraph_bench.py ${b}.sparse ${DIM[$b]} bool 20 2>&1 | grep -E "cuGraph|Error|error" | tail -3
  timeout 300 ./cgvenv/bin/python cugraph_bench.py ${b}.sparse ${DIM[$b]} tropical 20 2>&1 | grep -E "cuGraph|Error|error" | tail -3
done
echo CG_DONE > cg.done
