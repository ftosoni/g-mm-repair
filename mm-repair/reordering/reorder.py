#!/usr/bin/env python3
import subprocess, os.path, sys, argparse, time, struct

Description = """
Reorder the individual blocks of a matrix in CSRV format 

Currently supported reordering algorithms:
   pc:   PathCover
  mwm:   Maximum weighted matching
  lkh:   Lin-Kernighan
  pc+:   PathCover+
  
Usage examples: 
  reorder.py pc  /data/mm/covtype 581012 54 -b 8 -k 16
  reorder.py pcm /data/mm/covtype 581012 54 -b 8 -k 16
"""

Timelimit = 180000

# return the extension for multipart file
def filext_multipart(n,i):
  assert i<n, "Illegal parameters"
  if n==1:
    return ""
  return ".{tot}.{part}".format(tot=n,part=i)

# execute command: return True if everything OK, False otherwise
def execute_command(command):
  try:
    subprocess.check_call(command.split(),timeout=Timelimit)
  except subprocess.TimeoutExpired:
      # caso time out
      print("ERROR: no result after %d seconds. Command:" % Timelimit)
      print("\t"+ command)
      return False
  except subprocess.CalledProcessError:
    print("Error executing command:")
    print("\t"+ command)
    return False
  print("OK    ", command)
  return True

def execute_command_verbose(cmd, exit_code, msg='Something went wrong: please contact the maintainers') :
  print(":::executing command:", cmd)
  if execute_command(cmd):
    print('All done.')
  else:
    print(msg)
    sys.exit(exit_code)

# check that the input file and the corresponding .vco exist
# note .vco files are obtained from .vc if necessary and possible 
def check_input_files(fullname,n):
  if not os.path.isfile(fullname):
    print("ERROR File", fullname, "missing")
    return  False
  else: 
    print("OK     File",fullname)
  # test/convert vco files   
  for i in range(n):
    fex = fullname + filext_multipart(n,i)
    # if vco is there nothing to do
    if os.path.isfile(fex + ".vco"):
      print("OK     File",fex+".vco")
    # if vc and not vco copy it 
    elif os.path.isfile(fex + ".vc"):
      command = "cp -p {name}.vc {name}.vco".format(name=fex)
      if not execute_command(command):
        return False
    # if both are missing error
    else:
      print("ERROR: Files", fullname, ".vc/.vco missing")
      return  False
  return True
  

def show_command_line(f):
  f.write("=== Command line: ") 
  for x in sys.argv:
     f.write(x+" ")
  f.write("\n")   

def main():
  #args
  show_command_line(sys.stderr)
  parser = argparse.ArgumentParser(description=Description, formatter_class=argparse.RawTextHelpFormatter)
  parser.add_argument('algo',  help='algorithm to test: pc|mwm|lkh|pc+', type=str)  
  parser.add_argument('input', help='matrix file name (must be in csv format))', type=str)
  parser.add_argument('rows',  help='number of rows', type=int)
  parser.add_argument('cols',  help='number of columns', type=int)  
  parser.add_argument('-b', help='number of row blocks (default: 1)', type=int, default=1)
  parser.add_argument('-k', help='pruning parameter (default: 16)', type=int, default=16)
  parser.add_argument('-m', help='memory limit', type=int, default=0)
  parser.add_argument('-p', help='parallelism degree (default: 1)', type=int, default=1)
  args = parser.parse_args()

  #params
  algomap = {'pc':'cover', 'pc+':'cover2', 'mwm':'mwm', 'lkh':'tsp'}
  if args.algo not in algomap :
    algos_s = '|'.join(algomap.keys())
    print('Unknown algorithm: must be', algos_s)
    sys.exit(1)

  fullname = os.path.abspath(args.input)
         
  if args.b<1:
    print("ERROR: Invalid number of blocks (must be >0)") 
    sys.exit(3) 
     
  if args.rows<1 or args.cols<1:
    print("ERROR: Number of rows and columns must be positive") 
    sys.exit(4) 
     
  if not check_input_files(fullname,args.b):
    sys.exit(5)

  # arguments to be passed to reordering algorithms 
  cr_args = {'a_descr':args.algo, 'a':algomap[args.algo],'f':args.input, 'r':args.rows, 'c':args.cols, 'b':args.b, 
      'k':args.k, 'b_lst':args.b-1, 'pardegree':args.p, 'memlim':args.m, 
      'lkh_version':'LKH-3.0.7', 'lkh_link':'http://webhotel4.ruc.dk/~keld/research/LKH-3/LKH-3.0.7.tgz' }

  if args.b == 1 :
    reorder_matrix(cr_args)
  else :
    reorder_matrix_blocks(cr_args)
  
  return 0


# reorder a matrix taken as a single block
def reorder_matrix(cr_args) :

  is_csm_multithreaded = (cr_args['memlim']!=0 or cr_args['pardegree']!=1)

  if not is_csm_multithreaded :
    print('Transposing the matrix...')
    cmd = 'python3 ./column_major.py {f}'.format(**cr_args)
    execute_command_verbose(cmd, 7)
##
  print('Generating the CSM...')
  if is_csm_multithreaded :
    cmd = './build/tsp_generator_pruned_local_mt {f} {r} {c} {k} {pardegree} {memlim} '.format(**cr_args)
    execute_command_verbose(cmd, 8)
  else :
    cmd = './build/tsp_generator_pruned_local {f} {r} {c} {k} '.format(**cr_args)
    execute_command_verbose(cmd, 8)

##
  if cr_args['a_descr'] == 'lkh' :
    print('Checking LKH in local dir...')
    cmd = 'bash 09_get_lkh.sh {lkh_version} {lkh_link} '.format(**cr_args)
    execute_command_verbose(cmd, 9)
##
  print('Running', cr_args['a_descr'], '...')
  if False :
    pass
  elif cr_args['a_descr'] == 'pc' :
    cmd = 'python3 ./cover.py {f}.pruned_local_{k}.tsp '.format(**cr_args)
    execute_command_verbose(cmd, 10)
  elif cr_args['a_descr']=='pc+' :
    cmd = 'python3 ./cover2.py {f}.pruned_local_{k}.tsp '.format(**cr_args)
    execute_command_verbose(cmd, 10)
  elif cr_args['a_descr']=='mwm' :
    cmd = './build/mwm {f} pruned_local_{k} '.format(**cr_args)
    execute_command_verbose(cmd, 10)
    cmd = 'python3 ./mwm_sol_from_pairs.py {f} pruned_local_{k} '.format(**cr_args)
    execute_command_verbose(cmd, 10)
  elif cr_args['a_descr']=='lkh' :
    cmd = '{lkh_version}/LKH {f}.pruned_local_{k}.par'.format(**cr_args)
    execute_command_verbose(cmd, 10)
##
  print('Reordering the .vc file...')
  cmd = './vc_reorder.x {f} {r} {c} {f}.pruned_local_{k}.{a}.solution '.format(**cr_args)
  print(cmd)
  execute_command_verbose(cmd, 12)
##
  print('Cleaning...')
  cmd = 'rm -r '
  if not is_csm_multithreaded :
    cmd += '{f}_cols '.format(**cr_args)
  cmd += '{f}.pruned_local_{k}.par '.format(**cr_args)
  cmd += '{f}.pruned_local_{k}.tsp '.format(**cr_args)
  cmd += '{f}.pruned_local_{k}.{a}.solution '.format(**cr_args)
  #if cr_args['a_descr'] == 'lkh' :
  #  cmd += '{f}.pruned_local_{k}.log '.format(**cr_args)
  
  execute_command_verbose(cmd, 14)
  return   

# reoder a matrix splitted into b blocks
def reorder_matrix_blocks(cr_args) :

  is_csm_multithreaded = (cr_args['memlim']!=0 or cr_args['pardegree']!=1)

  print('Dividing into row chunks...')
  cmd ='python3 csv_splitter.py {f} {r} {b} '.format(**cr_args)
  execute_command_verbose(cmd, 6)
##
  if not is_csm_multithreaded :
    print('Transposing each row block...')
    cmd = 'bash 07_transpose.sh {f} {b} '.format(**cr_args)
    execute_command_verbose(cmd, 7)
##
  print('Generating the CSM for each row block...')
  if is_csm_multithreaded :
    cmd = 'bash 08_generate_csm_multithread.sh {f} {r} {c} {k} {pardegree} {memlim} {b} '.format(**cr_args)
    execute_command_verbose(cmd, 10)
  else :
    cmd = 'bash 08_generate_csm.sh {f} {r} {c} {b} {k}'.format(**cr_args)
    execute_command_verbose(cmd, 8)
##
  if cr_args['a_descr'] == 'lkh' :
    print('Checking LKH in local dir...')
    cmd = 'bash 09_get_lkh.sh {lkh_version} {lkh_link} '.format(**cr_args)
    execute_command_verbose(cmd, 9)
##
  print('Running', cr_args['a_descr'], 'upon each row block...')
  cmd = 'bash 10_run_{a}.sh {f} {b} {k} {lkh_version} '.format(**cr_args)
  execute_command_verbose(cmd, 10)
##
  print('Reordering each .vc file...')
  cmd = 'bash 11_reorder.sh {f} {r} {c} {b} {k} {a}'.format(**cr_args)
  execute_command_verbose(cmd, 11) 
##
  print('Cleaning...')
  for i in range(cr_args['b']) :
    cmd = 'rm -r '
    if not is_csm_multithreaded :
      cmd += '{f}.{b}.{i}_cols '.format(**cr_args, i=i)
    cmd += '{f}.{b}.{i}.pruned_local_{k}.par '.format(**cr_args, i=i)
    cmd += '{f}.{b}.{i}.pruned_local_{k}.tsp '.format(**cr_args, i=i)
    cmd += '{f}.{b}.{i}.pruned_local_{k}.{a}.solution '.format(**cr_args, i=i)
    if cr_args['a_descr'] == 'lkh' :
      cmd += '{f}.{b}.{i}.pruned_local_{k}.log '.format(**cr_args, i=i)
    execute_command_verbose(cmd, 12)
  return     

       
if __name__ == '__main__':
  main()
