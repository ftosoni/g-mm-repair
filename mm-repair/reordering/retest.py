#!/usr/bin/env python3
import subprocess, os.path, sys, argparse, time, struct
from pathlib import Path

Description = """
Split into row blocks and reorder a collection of files

Currently supported algorithms:
   pc:   path cover
  mwm:   maximum weighted matching
  lkh:   Lin-Kernighan
  pc+:   PathCover+
  
"""

Files = ['susy','higgs','airline78','covtype', 'census', 'optical', 'mnist2m']

Sizes = {'covtype':(581012, 54), 'census':(2458285, 68), 'optical':(325834, 174),
         'susy':(5000000, 18), 'higgs': (11000000,  28), 'mnist2m':(2000000,784),  
         'airline78':(14462943, 29)}

Data_dir = '../data/'
Logfile_name = "retest.log"
Time_exe = "/usr/bin/time"
Timelimit = 180000

# return the extension for multipart file
def filext_multipart(n,i):
  assert i<n, "Illegal parameters"
  if n==1:
    return ""
  return ".{tot}.{part}".format(tot=n,part=i)


# execute command: return True if everything OK, False otherwise
def execute_command(command,logfile):
  try:
    ris = subprocess.run(command.split(),stdout=logfile,
                         stderr=logfile,timeout=Timelimit,check=True)
  except subprocess.TimeoutExpired:
    # handle time out
      print("## ERROR: no result after %d seconds. Command:" % Timelimit)
      print("\t"+ command)
      return False
  except subprocess.CalledProcessError as ex:
    print("## Test failed: non-zero exit code ", ex.returncode,file=logfile)
    # write stdout/stderr to a separate logfile
    print("## Error executing:", command,file=logfile);
    return False
  except Exception as ex:
    print("## Error executing:", command,file=logfile);
    return False
  return True


def makerow_abs(f, a):
  s = "{name:10.9}& {col:<5}".format(name=f,col=Sizes[f][1])
  for p in a:
    s += "&{:>12} &{:12} &{:>12} &{:>12} ".format(p[0],p[1],p[2],p[3])
  s += "\\\\\n"
  return s

def makerow_per(f, a):
  s = "{name:10.9}& {col:<5}".format(name=f,col=Sizes[f][1])
  d = 8*Sizes[f][0]*Sizes[f][1]/100
  for p in a:
    s += "&{:8.4f} &{:8.4f} &{:8.4f} &{:8.4f} ".format(p[0]/d,p[1]/d,p[2]/d,p[3]/d)
  s += "\\\\\n"
  return s


def getsize_multipart(base,num,ext):
  tot = 0
  for i in range(num):
    name = base + filext_multipart(num,i)+ext
    tot += os.path.getsize(name)
  return tot

def show_command_line(f):
  f.write("=== Command line: ") 
  for x in sys.argv:
     f.write(x+" ")
  f.write("\n")   

def main():
  show_command_line(sys.stderr)
  parser = argparse.ArgumentParser(description=Description, formatter_class=argparse.RawTextHelpFormatter)
  parser.add_argument('algo',  help='algorithm to test: pc|mwm|lkh|pc+', type=str)  
  parser.add_argument('-d',    help='data directory (def. %s)' % Data_dir, type=str, default=Data_dir)
  parser.add_argument('-b',    help='number of row blocks (def 4)', default=4, type=int)
  parser.add_argument('-k',    help='pruning parameter (default: 16)', type=int, default=16)
  parser.add_argument('--files',help="colon separated list of file names eg matA:matB",type=str,default="")
  parser.add_argument('--sizes',help="colon separated list of matrix sizes eg. 1000x30:5000x20\n"
                                     "SIZES must contain the sizes of the matrices in FILES",type=str,default="")

  args = parser.parse_args()

  #params
  algolist = ['pc', 'mwm', 'lkh', 'pc+']
  if args.algo not in algolist :
    algos_s = '|'.join(algolist)
    print('Unknown algorithm: must be', algos_s)
    sys.exit(1)

  #check --sizes was not given without --files
  if len(args.sizes)>0 and len(args.files)==0:
    print("Error: option --sizes requires --files")
    sys.exit(1)
     
  absdir = os.path.abspath(args.d)
  if os.path.isdir(absdir):
    print("Data directory:", absdir)
  else:
    print("ERROR: Invalid data directory")
    sys.exit(2)

  # get file list if available
  global Files, Sizes
  if len(args.files)>0:
    Files = args.files.split(':')
    # check if argument --size was given as well
    if len(args.sizes)>0:
      sizes = args.sizes.split(':')
      if len(Files)!=len(sizes):
        print("Error: --files and --sizes must have the same number of elements")
        sys.exit(1)
      # add sizes to global dictionary Sizes
      for i in range(len(Files)):
        Sizes[Files[i]] = tuple([int(x) for x in sizes[i].split('x')])
        Sizes[Files[i]] = tuple([int(x) for x in sizes[i].split('x')])
    
  if args.b<1:
    print("ERROR: Invalid number of blocks (must be >0)") 
    sys.exit(3)    

  table_abs = ["### csrv and repair size; %d row-blocks\n" % args.b, 
           " file     & rows &        crsv &        re32 &        reiv &       reans \\\\\n"]   # latex table containing the results 
  table_per = ["### csrv and repair size (percentage); %d row-blocks\n" % args.b, 
           " file     & rows &    crsv &    re32 &    reiv &   reans \\\\\n"]   # latex table containing the results 
  with open(Logfile_name,"w",1) as logfile:
    print("Sending logging messages to file:", logfile.name)
    for f in Files:
      print("Processing file:", f)
      name  = os.path.join(args.d, f)
      rows,cols = Sizes[f]
      tablerow = []  # row of the results table
      
      # --- conversion to vc/val (CSRV) format      
      command = "../{exe} --norepair -b {blocks} {name} {r} {c}".format(
                  exe = "matrepair", blocks = args.b, name=name, r=rows, c=cols)
      if(not execute_command(command,logfile)):
        print("matrepair (1) failed")
        print("Check log file: " + logfile.name)
        sys.exit(4)
        
      # --- apply permutation
      command = "./{exe} -b {blocks} -k {prun} {alg} {name} {r} {c}".format(
                  exe = "reorder.py", alg=args.algo, blocks = args.b, 
                  prun = args.k, name=name, r=rows, c=cols)
      if(not execute_command(command,logfile)):
        print("reorder failed")
        print("Check log file: " + logfile.name)
        sys.exit(5)

      # --- complete compression     
      command = "../{exe} --noconv -b {blocks} {name} {r} {c}".format(
                  exe = "matrepair", blocks = args.b, name=name, r=rows, c=cols)
      if(not execute_command(command,logfile)):
        print("matrepair (2) failed")
        print("Check log file: " + logfile.name)
        sys.exit(6)
        
      # --- get size of compressed files   
      v = os.path.getsize(name+".val")
      vcsize = getsize_multipart(name,args.b,".vc") 
      csize = getsize_multipart(name,args.b,".vc.C") 
      rsize = getsize_multipart(name,args.b,".vc.R") 
      csizeiv = getsize_multipart(name,args.b,".vc.C.iv") 
      rsizeiv = getsize_multipart(name,args.b,".vc.R.iv") 
      ans_csize = getsize_multipart(name,args.b,".vc.C.ansf.1")
      tablerow.append((v+vcsize,v+csize+rsize,v+csizeiv+rsizeiv,
                      v+ans_csize+rsizeiv))
      # tests for current file completed
      table_abs.append(makerow_abs(f, tablerow))
      table_per.append(makerow_per(f, tablerow))
    # output tables  
    for s in table_abs:
      print(s,end="")
    print()
    for s in table_per:
      print(s,end="")
  return

       
if __name__ == '__main__':
  print("Setting the working directory to:", Path(__file__).resolve().parent) 
  os.chdir(Path(__file__).resolve().parent)
  main()
