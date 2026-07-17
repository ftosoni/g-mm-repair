#!/usr/bin/env python3

import sys, time, argparse, math, os.path, struct

Description = """
Tool to create a vector and write it to a file in binary format

The default format is 64 bit floats; use -i/-f to write 32 bit integers/floats 
"""


def main():
  parser = argparse.ArgumentParser(description=Description, formatter_class=argparse.RawTextHelpFormatter)
  parser.add_argument('outfile', help='output file name', type=str)
  parser.add_argument('size',    help='number of elements', type=int)
  parser.add_argument('val',      help='values to be repeated cyclically', nargs='+')
  parser.add_argument('-f', help='vector entries are 32 bit floats',action='store_true')
  parser.add_argument('-i', help='vector entries are 32 bit integers',action='store_true')
  #parser.add_argument('-v',  help='verbose',action='store_true')
  args = parser.parse_args()
  if args.i and args.f:
    print("Error: Options -f and -i are mutually exclusive");
    return 
  
  with open(args.outfile,"wb") as f:
    for i in range(args.size):
      x = args.val[ i%len(args.val)]
      if args.i:   f.write(struct.pack("<i", int(x)))
      elif args.f: f.write(struct.pack("<f", float(x)))
      else:        f.write(struct.pack("<d", float(x)))
  print("==== Done")



if __name__ == '__main__':
    main()
