#!/usr/bin/env python3

import sys, time, argparse, subprocess, os.path, struct

Description = """
Given a matrix (not necessarily of numbers) in csv format,
replace the missing or NA entries with a given string 
"""

shasum_exe = "sha256sum"

def main():
  parser = argparse.ArgumentParser(description=Description, formatter_class=argparse.RawTextHelpFormatter)
  parser.add_argument('input', help='input file name', type=str)
  parser.add_argument('cols', help='number of columns', type=int)
  parser.add_argument('value', help='entry to be used for missing values', type=str)
  args = parser.parse_args()
    
  # take start time 
  start0 = start = time.time()
  with open(args.input,"rt") as f:
    outname = args.input + ".missing.%s" %(args.value) 
    outmode = "wt"
    with open(outname,outmode) as g:
      nonz = 0
      wr = 0 # number of written rows
      for line in f:
        wr += 1
        # get items in csv file
        a = [p.strip() for p in line.split(",")]
        if len(a)!=args.cols:
          print("Row", wr, "has", len(a), "columns")
          sys.exit(1) 
        b = [args.value if (len(p)==0 or p=='NA') else p for p in a]
        print(",".join(b),file=g)
  print("Written", wr, "rows")
  print("==== Done")


if __name__ == '__main__':
    main()
