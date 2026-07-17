#!/usr/bin/env python3

import sys, argparse

Description = """
Strip the initial or final columns and/or the inital rows
from a csv matrix producing a smaller csv matrix.

The stripping is done at the textual level: non-stripped values 
are not modified (eg integers remain integers, they are not 
converted to floats with a 0 fractional part). 
"""

def main():
  parser = argparse.ArgumentParser(description=Description, formatter_class=argparse.RawTextHelpFormatter)
  parser.add_argument('input', help='input file name', type=str)
  parser.add_argument('rows', help='number of expected rows in outfile', type=int)
  parser.add_argument('cols', help='number of expected columns in outfile', type=int)
  parser.add_argument('-c', help='initial cols to skip, if <0 skip final cols (def. 0)',type=int,default=0 )
  parser.add_argument('-r', help='initial rows to skip (def. 0)',type=int,default=0 )
  parser.add_argument('-o', help='output file name (def. input.RxC)',type=str,default="" )
  args = parser.parse_args()
  # minimal check on the input values
  if args.c==0 and args.r==0:
    print("Nothing to strip. Exiting...");
    return 
    
  with open(args.input,"rt") as f:
    # set output file name and output file mode
    if args.o!="": 
      outname = args.o
    else:
      outname = args.input + ".%dx%d" %(args.rows,args.cols)
    with open(outname,"w") as g:
      r = 0 # number of read rows
      wr = 0 # number of written rows
      for line in f:
        r += 1
        if r>args.r: # skip first args.r rows
          a = line.rstrip().split(",")
          if args.c>0:
            a = a[args.c:] # remove initial args.c columns
          elif args.c<0:
            a = a[:args.c] # remove final args.c columns
          if len(a)!=args.cols:
            print("row", r,"has", len(a), "elements")
            sys.exit(1);
          # create and write to g stripped row
          print(",".join(a),file=g)
          wr += 1
  if wr!=args.rows:
    print("Warning! Written", wr, "rows instead of", args.rows)
  print("==== Done")


if __name__ == '__main__':
    main()
