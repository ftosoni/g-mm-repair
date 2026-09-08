#!/bin/bash
cd "$(dirname "$0")"
Z=${ZENODO_DIR:-zenodo}/wikidata   # Wikidata .sparse edge lists from the Zenodo package (sec 2.D)
declare -A DIM=( [wd_sports_team]="332121 29854" [wd_cast_member]="173977 144095" [wd_citizenship]="2874250 2556" [wd_occupation]="3459933 10610" [wd_subclass_of]="1487709 73417" )
for b in wd_sports_team wd_cast_member wd_citizenship wd_occupation wd_subclass_of; do
  echo "##### $b #####"
  ./gbvenv/bin/python graphblas_bench.py $Z/${b}.sparse ${DIM[$b]} bool 50 2>&1
  ./gbvenv/bin/python graphblas_bench.py $Z/${b}.sparse ${DIM[$b]} tropical 50 2>&1
done
echo GB_DONE > gb.done
