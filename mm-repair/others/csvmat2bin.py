#!/usr/bin/env python3

import sys, argparse, subprocess, struct

Description = """
Tool to convert a matrix written in text csv format (one line per row)
into binary form, ie rows x cols float64 (default) or float32 (option -f) 
or int32 (option -i).
Optionally, a set of leading rows and/or columns can be removed before 
the conversion. If the number of columns to be removed is negative, 
minus that number of trailing columns are removed. 

All matrix entries are represented so the outfile has size 
    rows*cols*sizeof(entry). 
Note that when using the int32 or float32 output format some information 
will be lost if the input values are not of the right type.
"""

shasum_exe = "sha1sum"

def main():
  parser = argparse.ArgumentParser(description=Description, formatter_class=argparse.RawTextHelpFormatter)
  parser.add_argument('input', help='input file name', type=str)
  parser.add_argument('rows', help='number of expected rows in outfile', type=int)
  parser.add_argument('cols', help='number of expected columns in outfile', type=int)
  parser.add_argument('-c', help='initial cols to skip, if <0 skip final cols (def. 0)',type=int,default=0 )
  parser.add_argument('-r', help='initial rows to skip (def. 0)',type=int,default=0 )
  parser.add_argument('-o', help='output file name (def. input.[if]bin)',type=str,default="" )
  parser.add_argument('-f', help='output values as single precision floats',action='store_true')
  parser.add_argument('-i', help='output values as int32 (trucate input values)',action='store_true')
  parser.add_argument('--sum', help='show input/output file shasum',action='store_true')
  #parser.add_argument('-v',  help='verbose',action='store_true')
  args = parser.parse_args()
  # minimal check on the input values
  if args.f and args.i:
    print("Error: Options -f and -i are mutually exclusive");
    return 
  
  # report input file name and input file hash digest
  if args.sum:
    print("SHA1  Infile:", file_digest(args.input), args.input);
  with open(args.input,"rt") as f:
    # set output file name and output file mode
    if args.o!="": 
      outname = args.o
    elif args.f: outname = args.input + ".fbin"
    elif args.i: outname = args.input + ".ibin" 
    else: outname = args.input + ".bin"
    with open(outname,"wb") as g:
      nonz = 0
      r = 0 # number of read rows
      wr = 0 # number of written rows
      values  = set()
      for line in f:
        r += 1
        if r>args.r: # skip first args.r rows
          a = line.rstrip().split(",")
          if args.c>0:
            a = a[args.c:] # remove initial args.c columns
          elif args.c<0:
            a = a[:args.c] # remove final args.c columns
          if len(a)!=args.cols:
            for i in range(len(a)):
              print(i,a[i])
            print("row", r,"has", len(a), "elements")
            sys.exit(1);
          ## convert all string values to float
          if args.i: 
            b = [int(x) for x in a]
          else:
            b = [float(x) for x in a]
          for x in b:
            if x!=0:
              nonz +=1
              values.add(x)
          if args.f: g.write(struct.pack("%df" % len(b), *b)) # single precision 
          elif args.i: g.write(struct.pack("%di" % len(b), *b)) # int32
          else:      g.write(struct.pack("%dd" % len(b), *b)) # double precision
          wr += 1
  if args.sum:        
    print("SHA1 Outfile:", file_digest(outname), outname);
  if wr!=args.rows:
    print("Warning! Written", wr, "rows instead of", args.rows)
  print("Number of nonzeros:",nonz," Nonzero ratio: %.4f" % (nonz/(wr*args.cols)))  
  print(len(values), "Distinct nonzeros values")  
  print("==== Done")


# compute hash digest for a file 
def file_digest(name):
    try:
      hash_command = "{exe} {infile}".format(exe=shasum_exe, infile=name)
      hashsum = subprocess.check_output(hash_command.split())
      hashsum = hashsum.decode("utf-8").split()[0]
    except:
      hashsum = "Error!" 
    return hashsum  

if __name__ == '__main__':
    main()
