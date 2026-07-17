#!/usr/bin/env python3
import sys, time, argparse, math, os.path, struct, array

# input file extension indicating that the input matrix is in binary 
dbl_extension = ".dbl"

Description = """
Tool to convert a matrix written in csv or double format into
the CRSV representation, generating a .vc file and a .val file

If the input file ends with {ext} then it is assumed to be in binary form 
with each entry represented as an 8 byte double; otherwise the matrix is
assumed to be in csv textual format with one line per row
""".format(ext=dbl_extension)



def main():
  show_command_line(sys.stderr)  
  parser = argparse.ArgumentParser(description=Description, formatter_class=argparse.RawTextHelpFormatter)
  parser.add_argument('input', help='input file name', type=str)
  parser.add_argument('rows', help='number of rows', type=int)
  parser.add_argument('cols', help='number of columns', type=int)
  parser.add_argument('-b', help='number of row blocks (def. 1)', type=int, default=1)
  parser.add_argument('-c', help='initial columns to skip (def. 0)',type=int,default=0 )
  parser.add_argument('-r', help='initial rows to skip (def. 0)',type=int,default=0 )
  parser.add_argument('-f', help='store matrix entries as 32-bit floats',action='store_true')
  parser.add_argument('-i', help='store matrix entries as 32-bit integers',action='store_true')
  #parser.add_argument('--sum', help='compute output file shasum',action='store_true')
  #parser.add_argument('-v',  help='verbose',action='store_true')
  args = parser.parse_args()
  if args.i and args.f:
    print("Error: Options -f and -i are mutually exclusive",file=sys.stderr);
    return 
  # establish if the input is csv or binary   
  if args.input.endswith(dbl_extension):
    read_mode = "rb"
    print("Input file format is binary",file=sys.stderr)
    csv = False
  else:
    read_mode = "rt"
    print("Input file format is csv",file=sys.stderr)
    csv = True
  assert args.b >0, "Number of row blocks must be greater than 0"
  assert args.b<= args.rows, "Too many row blocks!"
  # compute size of each block 
  block_size = (args.rows+args.b-1)//args.b
  print("Block size", block_size);
  
  
  start = time.time()
  with open(args.input,read_mode) as f:
    r = 0   # number of read rows
    wr = 0  # number of written rows
    nonz = 0        # total number of nonzeros  
    maxcode = 0     # largest code in a .vc file
    values  = {}    # dictionary of distinct nonzero 
    outname_val = args.input +".val"
    with open(outname_val,"wb") as g_val:
      for bn in range(args.b):
        outname = args.input + filext_multipart(args.b,bn)+".vc"
        with open(outname,"wb") as g:
          # read one line at a time from f
          while True:
            # read a text or binary matrix row 
            if csv:
              line = f.readline()
              if not line:
                break
            else:
              b = array.array('d')   # create empty array
              try:
                b.fromfile(f,args.cols)
              except EOFError as e:
                break
            # check row 
            r += 1
            if wr>=args.rows:
              print("Warning: more matrix rows than expected",file=sys.stderr)
              break 
            if r<=args.r: continue   # skip first args.r rows
            # convert text numerical values
            if csv:
              a = line.rstrip().split(",")
              a = a[args.c:] # remove initial arg.c columns
              if len(a)!=args.cols:
                # row with wrong number of elements: print error msg and exit
                print(a,file=sys.stderr); print("row", r,"has", len(a), "elements",file=sys.stderr)
                sys.exit(1)
              ## convert to double or integer (the latter if option -i was given)
              if args.i: b = [  int(s) for s in a]
              else:      b = [float(s) for s in a]
            assert len(b)==args.cols
            for i in range(len(b)):
              # get (or create) integer id associated to value x=b[i]
              x = b[i]
              if x!=0:
                nonz +=1
                if x not in values:
                  newid = len(values)
                  values[x] = newid
                  if   args.i: g_val.write(struct.pack("<i", x)) # 32-bit int
                  elif args.f: g_val.write(struct.pack("<f", x)) # 32-bit float
                  else:        g_val.write(struct.pack("<d", x)) # 64-bit double 
                assert x in values
                # create code combining value id and column number
                code = values[x]*args.cols + i 
                code += 1              # shift by 1 to allow code for endrow code 0 
                if code>= 2**30:
                  print("Code", code, "larger than 2**30. We are in trouble",file=sys.stderr)
                  sys.exit(1)
                if code>maxcode: maxcode=code
                g.write(struct.pack("<I", code))
            # row wr completed, write unique end of row code, and update wr
            g.write(struct.pack("<I", 0)) # not-so-unique EOR code (was wr now is 0 for all rows)
            wr += 1
            # if a block is full stop and start the next 
            if (wr%block_size) == 0: break
  if wr!=args.rows:
    print("Warning! Written", wr, "rows instead of", args.rows,file=sys.stderr)
  print("Elapsed time: {0:.4f} secs".format(time.time()-start),file=sys.stderr);
  print("Number of nonzeros:",nonz," Nonzero ratio: %.4f" % (nonz/(wr*args.cols)),file=sys.stderr)  
  print(len(values), "distinct nonzeros values",file=sys.stderr)
  print("Largest codeword:", maxcode, " bits:", math.ceil(math.log(1+maxcode,2)),file=sys.stderr)
  print("==== Done",file=sys.stderr)


# return the extension for multipart file
def filext_multipart(n,i):
  assert i<n, "Illegal parameters"
  if n==1:
    return ""
  return ".{tot}.{part}".format(tot=n,part=i)


def show_command_line(f):
  f.write("==== Command line:\n  ") 
  for x in sys.argv:
     f.write(x+" ")
  f.write("\n")   

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
