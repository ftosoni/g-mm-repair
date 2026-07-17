#!/bin/bash
cd "$(dirname "$0")"
echo "=== generate YAGO sparse edge lists ==="
python3 process_wikidata.py yago-2s.dat yago-2s.dat.P sparse 8:yago_plays_for 16:yago_is_citizen_of 15:yago_acted_in
declare -A DIM=( [yago_plays_for]="91397 13765" [yago_is_citizen_of]="45823 309" [yago_acted_in]="26665 39559" [wd_sports_team]="332121 29854" [wd_cast_member]="173977 144095" [wd_citizenship]="2874250 2556" [wd_occupation]="3459933 10610" [wd_subclass_of]="1487709 73417" )
for b in yago_plays_for yago_is_citizen_of yago_acted_in wd_sports_team wd_cast_member wd_citizenship wd_occupation wd_subclass_of; do
  echo "##### $b #####"
  ./gbvenv/bin/python graphblas_bench.py ${b}.sparse ${DIM[$b]} bool 50 2>&1
  ./gbvenv/bin/python graphblas_bench.py ${b}.sparse ${DIM[$b]} tropical 50 2>&1
done
echo GB_DONE > gb.done
