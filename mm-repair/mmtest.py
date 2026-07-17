#!/usr/bin/env python3
import subprocess, os.path, sys, argparse, time, struct

Description = """
Tool to create a Latex table containing the results of a set of experiments

Currently supported tests:
   mg: compression with gzip/xz
   mz: compression to CSRV format followed by RePair grammar compression
   md: compression to DRV  format followed by RePair grammar compression
   mm: test matrix-vector multiplication algorithms (only for f64 matrices)

By default the input matrix is assumed to contain doubles in csv format,
the options --i32, --f32, and --f64 signal that the input is in binary
format with entries respectively of type int32, float32 and float64.
"""

Files = ['susy','higgs','airline78','covtype', 'census', 'optical', 'mnist2m', 'imagenet']

Data_dir = 'data/'
Logfile_name = "errors.log"
Time_exe = "/usr/bin/time"

Sizes = {'covtype':(581012, 54), 'census':(2458285, 68), 'optical':(325834, 174),
         'susy':(5000000, 18), 'higgs': (11000000,  28), 'mnist2m':(2000000,784),  
         'airline78':(14462943, 29),'imagenet':(1262102, 900), 
         'glove.q64.csv':(400000,300), 'glove.q256.csv':(400000,300),'glove.q1k.csv':(400000,300)}

# matrix mutiplication algorithms to be tested
Algos = ['csrvmm', 're32mm','reivmm','remm']

# name of files containing the input/output vectors
# these files are created by this script at each execution
Xvname = "x1.dbl"
Yvname = "y.dbl"
Zvname = "z.dbl"
Evname = "ein.dbl"
Timelimit = 36000
TmpFilename = "tmp_mmtest"
# name (and options and extension) for zip-like compression
zip1 = "gzip"
zip1ext = ".gz"
# zip2 = "xz --lzma2=dict=1GiB" works but requires more than 10GiB RAM ...
# zip2 = "xz --lzma2=dict=128MiB" # require approx 1.2GiB RAM
zip2 = "zstd"
zip2ext = ".zst"

# check that the test files exist and sizes are defined
def check_testfiles(args,sufxs):
  ok = True
  for f in Files:
    if f not in Sizes:
      print("  Size for test file", f, "missing",file=sys.stderr)
      ok = False
    
    for sx in sufxs:
      path = os.path.join(args.d,f + sx)
      if not os.path.exists(path):
        print("  Test file", path, "missing",file=sys.stderr)
        ok = False
  if not ok:
    print("Check -d option and the constants Files and Sizes",file=sys.stderr)
    sys.exit(1)


# write to file a vector cointaining cols copies of value
# in big-endian float64 format
def createx(cols,value=1):
  with open(Xvname,"wb") as f:
    for i in range(cols):
      f.write(struct.pack("<d", float(value)))


# convert a csv file to binary and compress it with gzip and xz
def test_gzip(args,logfile):
  global TmpFilename
  table = [f"### {zip1} and {zip2} size vs dense uncompressed size (absolute and percentage)\n", 
           f" file     & rows &   dense size % &&     {zip1} % &&   {zip2} % &\\\\\n"]
  for f in Files:
    name= os.path.join(args.d,f)
    exe_name = os.path.join(args.main_dir,"others/csvmat2bin.py")

    rows,cols = Sizes[f]
    tablerow = []  # row of the results table
    try:
      # convert csv matrix to binary if necessary
      if not (args.i32 or args.f32 or args.f64):
        command = f"{exe_name} {name} -o {TmpFilename} {rows} {cols}"
        ris = subprocess.run(command.split(),stdout=subprocess.PIPE,
                           stderr=subprocess.PIPE,timeout=Timelimit,check=True)
      else:
        TmpFilename = name
      # check number of rows and columns
      if os.path.getsize(TmpFilename)!=rows*cols*args.entry_size:
        print(f"# of rows/columns and entry size don't match with {TmpFilename} file size")
        sys.exit(1)
      # run gzip and xz
      command = f"{zip1} -kf {args.extra} {TmpFilename}"            
      ris = subprocess.run(command.split(),stdout=subprocess.PIPE,
                           stderr=subprocess.PIPE,timeout=Timelimit,check=True)
      command = f"{zip2} -kf {args.extra} {TmpFilename}"            
      ris = subprocess.run(command.split(),stdout=subprocess.PIPE,
                           stderr=subprocess.PIPE,timeout=Timelimit,check=True)
    except subprocess.TimeoutExpired:
      # handle time out
      print("Command:", command)      
      print(" Test failed: no result after %d seconds" % Timelimit)
      continue
    except subprocess.CalledProcessError as ex:
      print("Command:", command)
      print(" Test failed: non-zero exit code ", ex.returncode)
      # write stdout/stderr to a separate logfile
      print("## Error executing:", command,file=logfile);
      print("## stdout:\n", ex.stdout.decode() ,file=logfile);
      print("## stderr:\n", ex.stderr.decode() ,file=logfile);
      sys.exit(2)
    except Exception as ex:
      print("Command:", command)      
      print(" Test failed:", str(ex))
      sys.exit(2)
    tablerow.append((os.path.getsize(TmpFilename),os.path.getsize(TmpFilename+zip1ext),
                     os.path.getsize(TmpFilename+zip2ext)))
    # elapsed = int(ris.stdout.split()[-2])
    # peakmem = int(ris.stderr.split()[3])
    # tablerow.append((a, elapsed/n,peakmem, e))
    # tests for current file completed
    table.append(makerow_mgzip(args,f, tablerow))
  # all files processed
  return table


# compress with matrepair obtaining CSRV and grammar representation
# if drv is True use the drv representation
def test_compress(args, logfile, drv=False):
  # set different behavior for no-column-id option 
  if drv:
    args.name = "drv "
    args.mext = ".dv"
    args.val_ext += "d"
    args.extra += " --drv"
  else:
    args.name = "csrv"           
    args.mext = ".vc"
  # init latex table containing the results
  table = [f"### {args.name} + repair +iv/ans size; {args.b} row-blocks\n", 
           f" file     & rows &        {args.name} &        re32 &        reiv &       reans & reans%\\\\\n"]  
  for f in Files:
    name  = os.path.join(args.d,f)
    exe_name = os.path.join(args.main_dir,"matrepair")
    rows,cols = Sizes[f]
    tablerow = []  # row of the results table
    command = f"{exe_name} -r -y -b {args.b} -p {args.p} {args.extra} {args.bin} {name} {rows} {cols}"
    try:
      ris = subprocess.run(command.split(),stdout=logfile,
                           stderr=logfile,timeout=Timelimit,check=True)
    except subprocess.TimeoutExpired:
      # caso time out
      print("Command:", command)
      print(" Test failed: no result after %d seconds" % Timelimit)
      continue
    except subprocess.CalledProcessError as ex:
      print("Command:", command)
      print(" Test failed: non-zero exit code ", ex.returncode)
      # write stdout/stderr to a separate logfile
      print("## Error executing:", command,file=logfile);
      print("## stdout:\n", ex.stdout.decode() ,file=logfile);
      print("## stderr:\n", ex.stderr.decode() ,file=logfile);
      sys.exit(2)
    except Exception as ex:
      print("Command:", command)
      print(" Test failed:", str(ex))
      sys.exit(2)
    v = os.path.getsize(name+args.val_ext)
    vcsize = getsize_multipart(name,args.b,args.mext) 
    csize = getsize_multipart(name,args.b,args.mext+".C") 
    rsize = getsize_multipart(name,args.b,args.mext+".R") 
    csizeiv = getsize_multipart(name,args.b,args.mext+".C.iv") 
    rsizeiv = getsize_multipart(name,args.b,args.mext+".R.iv") 
    ans_csize = getsize_multipart(name,args.b,args.mext+".C.ansf.1")
    tablerow.append((v+vcsize,v+csize+rsize,v+csizeiv+rsizeiv,
                    v+ans_csize+rsizeiv))
    # tests for current file completed
    table.append(makerow_mz(args,f, tablerow))
  # all files processed
  return table

# return the extension for multipart file
def filext_multipart(n,i):
  assert i<n, "Illegal parameters"
  if n==1:
    return ""
  return ".{tot}.{part}".format(tot=n,part=i)

def getsize_multipart(base,num,ext):
  tot = 0
  for i in range(num):
    name = base + filext_multipart(num,i)+ext
    tot += os.path.getsize(name)
  return tot



# test running times and space for matrix multiplication 
def test_time(args,logfile):
  # build latex table containing the reuslts  
  table = ["### time x iteration and peak memory usage in kb for matrix multiplication; %d row-blocks\n" % args.b, 
           " file     & cols "]
  for a in Algos:
    table[1] += "& {name:12.9}&".format(name=a)         
  table[1] += "\\\\\n" 
  for f in Files:
    name  = os.path.join(args.d, f)
    rows,cols = Sizes[f]
    createx(cols)  # create file containing x vector
    tablerow = []  # row of the results table
    for a in Algos:
      exe_name = os.path.join(args.main_dir,a)
      # only save the eigenvalue
      command = "{exe} -n {num} -b {blocks} -e {ename} -z {z} {name} {r} {c} {x}".format(ename=Evname,
                exe = exe_name, num=args.n, blocks=args.b, name=name, r=rows, c=cols, x=Xvname, z=f+Zvname)
      command = Time_exe + ' -f%e:%M ' + command          
      try:
        ris = subprocess.run(command.split(),timeout=Timelimit,check=True,
              stdout=subprocess.PIPE,stderr=subprocess.PIPE)
      except subprocess.TimeoutExpired:
        # caso time out
        print("Command:", command)
        print(" Test failed: no result after %d seconds" % Timelimit)
        continue
      except subprocess.CalledProcessError as ex:
        print("Command:", command)
        print(" Test failed: non-zero exit code ", ex.returncode)
        # write stdout/stderr to a separate logfile
        print("## Error executing:", command,file=logfile);
        print("## stdout:\n", ex.stdout.decode() ,file=logfile);
        print("## stderr:\n", ex.stderr.decode() ,file=logfile);
        sys.exit(2)
      except Exception as ex:
        print("Command:", command)
        print(" Test failed:", str(ex))
        sys.exit(2)
      # convert /bin/time output to elapsed, peak memory   
      timespace= str(ris.stderr,'utf-8').split()[-1].split(":")
      elapsed = float(str(timespace[0]))
      peakmem = int(str(timespace[1]))
      with open(Evname,"rb") as evf:
        e = struct.unpack("d",evf.read(8))[0]
      # a = algo, elapsed e=eigenvalue
      tablerow.append((a, elapsed/args.n,peakmem, e))
    # tests for current file completed
    table.append(makerow_mm(f, tablerow))
  # all files processed
  return table

def makerow_mm(f, a):
  s = "{name:10.9}& {col:<5}".format(name=f,col=Sizes[f][1])
  for p in a:
    # no eigenvalue
    # s += "&{:6.2f} &{:4.0f}  ".format(p[1],p[2]/1000000)
    # with eigenvalue
    s += "&{:6.2f} &{:4.0f}&{:.5g} ".format(p[1],p[2],p[3])
  s += "\\\\\n"
  return s

def makerow_mgzip(args,f, a):
  s = "{name:10.9}& {col:<5}".format(name=f,col=Sizes[f][0])
  d = args.entry_size*Sizes[f][0]*Sizes[f][1]/100
  for p in a:
    s += "&{:11.0f} &{:6.2f} &{:11.0f} &{:6.2f} &{:11.0f} &{:6.2f}".format(
          p[0],p[0]/d,p[1],p[1]/d,p[2],p[2]/d)
  s += "\\\\\n"
  return s


def makerow_mz(args,f, a):
  s = "{name:10.9}& {col:<5}".format(name=f,col=Sizes[f][0])
  d = args.entry_size*Sizes[f][0]*Sizes[f][1]/100
  for p in a:
    s += "&{:>12} &{:12} &{:>12} &{:>12} &{:6.2f}".format(p[0],p[1],p[2],p[3],p[3]/d)
  s += "\\\\\n"
  return s



def show_command_line(f):
  f.write("=== Command line: ") 
  for x in sys.argv:
     f.write(x+" ")
  f.write("\n")   

def main():
  show_command_line(sys.stderr)
  parser = argparse.ArgumentParser(description=Description, formatter_class=argparse.RawTextHelpFormatter)
  parser.add_argument('op', help='operation to test: mg|mz|md|mm', type=str)
  parser.add_argument('-d', help='data directory (def. %s)' % Data_dir, type=str, default=Data_dir)
  parser.add_argument('-n', help='number of iterations (def 3)', default=3, type=int)  
  parser.add_argument('-b', help='number of row blocks (def 1)', default=1, type=int)    
  parser.add_argument('-p', help='number of parallel procs (def. 1)',type=int, default=1)
  parser.add_argument('--i32', help='input contains int32 entries in binary format',action='store_true')
  parser.add_argument('--f32', help='input contains float32 entries in binary format',action='store_true')
  parser.add_argument('--f64', help='input contains float64 entries in binary format',action='store_true')
  parser.add_argument('--files',help="colon separated list of file names eg matA:matB",type=str,default="")
  parser.add_argument('--sizes',help="colon separated list of matrix sizes eg. 1000x30:5000x20\n"
                                     "SIZES must contain the sizes of the matrices in FILES",type=str,default="")
  parser.add_argument('--extra',help="extra options to pass to the tested tool",type=str,default="")
  args = parser.parse_args()

  if (args.i32 and args.f32) or (args.i32 and args.f64) or (args.f64 and args.f32):
    print("Options --i32 --f32 and --f64 are incompatible")
    sys.exit(1)
  args.entry_size = 8
  args.val_ext = ".val"
  args.bin = ""
  if args.i32:
    args.val_ext = ".ival"
    args.bin = "--i32" 
    args.entry_size = 4   # size in byte of a int32
  elif args.f32:
    args.val_ext = ".fval"
    args.bin = "--f32"
    args.entry_size = 4   # size in byte of a float32
  elif args.f64:
    args.bin = "--f64" # matrix of doubles in binary form
  
  #check --sizes was not given without --files
  if len(args.sizes)>0 and len(args.files)==0:
    print("Error: option --sizes requires --files")
    sys.exit(1)

  # check data directory   
  if not  os.path.isdir(args.d):
    print("Invalid data directory: "+ args.d,file=sys.stderr) 
    sys.exit(1)
  args.d = os.path.abspath(args.d)
  
  # get directory containing this file   
  args.main_dir = os.path.dirname(os.path.abspath(__file__))
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

  # run the task   
  s1 = time.time()
  print("Sending logging messages to file:", Logfile_name)  
  with open(Logfile_name,"a", 1) as logfile: # 1 for line buffering: keeps the logs ordered
    if args.op=='mm':    # matrix multiplication (supported only for f64)
      check_testfiles(args,[".val"])
      table = test_time(args, logfile)
    elif args.op=='mg':   # matrix zip-like compression
      table = test_gzip(args, logfile)
    elif args.op=='mz' or args.op=='md':   # matrix compression
      table = test_compress(args,logfile, args.op=='md')
    else: 
      print("Unknown operation! Must be mg, mz, md, or mm",file=sys.stderr)
      sys.exit(1)
  e1 = time.time()
  print("Elapsed time: %.3f\n" % (e1-s1),file=sys.stderr)
  if args.op=='mm' or args.op=='mg' or args.op=='mz' or args.op=='md':  
    for s in table:
      print(s,end="")


  
if __name__ == '__main__':
  main()
