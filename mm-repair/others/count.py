#!/usr/bin/env python3

import sys, time, argparse, math, os.path, struct

Description = """
Count the number of occurrences of a set of values in a binary file 
"""


def main():
  parser = argparse.ArgumentParser(description=Description, formatter_class=argparse.RawTextHelpFormatter)
  parser.add_argument('infile', help='input file name', type=str)
  parser.add_argument('val',    help='values to be counted', type=int, nargs='+')
  parser.add_argument('-i',   help='input integer size (def. 4)', type=int,default=4)
  parser.add_argument('-b',   help="input is big endian (def false)",  action='store_true')
  parser.add_argument('-v',  help='verbose',action='store_true')
  args = parser.parse_args()
  # -- check
  isize = args.i
  if not 0<isize<=8:
    sys.sdterr.write("Invalid input integer size, must be in [1,8]\n");
    return
  if args.b!=0 and not(0<args.b<=8):
    sys.sdterr.write("Invalid output integer size, must be in [1,8]\n");
    return
  # -- init counter
  count = {}
  for x in args.val: count[x] = []  

  # read data from file 
  sfill =  bytearray([0]*(8-isize))
  tot = 0
  with open(args.infile,"rb") as f:
    while True:
      s = f.read(isize)
      if len(s)==0:
        break
      if len(s)!=isize:
        raise "Input file len not multiple of " + str(isize)
      # decode and update count   
      if(args.b):
        # big endian (untested)
        v = struct.unpack('>Q',sfill+s)[0]
      else:
        # little endian
        v = struct.unpack('<Q',s+sfill)[0]
      if v in count:
        count[v].append(tot)
      tot += 1
  for i in count:
    if args.v: print( "%d:" % i,  count[i])
    else: print("%d: %d" % (i, len(count[i])))
  print(tot, "integers read");
  print("==== Done")



if __name__ == '__main__':
    main()
