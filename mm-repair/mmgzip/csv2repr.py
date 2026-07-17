import struct, sys
import argparse
import subprocess

Description = """
Convert a CSV into a block-wise dense representation.
  
Example usage: 
  reorder.py pc /data/mm/covtype 581012 16 unco
"""

# execute command: return True if everything OK, False otherwise
def execute_command(command):
  Timelimit = 180000
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
  if execute_command(cmd):
    print('All done.')
  else:
    print(msg)
    sys.exit(exit_code)

def get_nums(infilepath, sep=',') :
    infile = open(infilepath, 'r')
    line = infile.readline()
    count = 0
    while line :
        assert(line[-1] == '\n')
        tokens = line[:-1].split(sep)
        for token in tokens :
            yield count, float(token)
        #update
        line = infile.readline()
        count += 1
    infile.close()

def to_byte(d) :
    return struct.pack('d', d)

def show_command_line(f):
  f.write("=== Command line: ") 
  for x in sys.argv:
     f.write(x+" ")
  f.write("\n")   

if __name__ == '__main__' :
    #args
    show_command_line(sys.stderr)
    parser = argparse.ArgumentParser(description=Description, formatter_class=argparse.RawTextHelpFormatter)
    parser.add_argument('input', help='matrix file name (must be in csv format)', type=str)
    parser.add_argument('rows',  help='number of rows', type=int)
    parser.add_argument('blocks', help='number of row blocks', type=int)
    parser.add_argument('repr', help='representation (unco|gzip)', type=str)
    args = parser.parse_args()
    if args.repr not in ['unco','gzip'] :
        print("The repr field must be either 'unco' or 'uncompressed'")
        exit(-1)
    
    infilepath = args.input 
    rows = args.rows 
    blocks = args.blocks 

    min_rows = (rows + blocks - 1 ) // blocks
    end_rows = 0
    curr_block = 0

    outfile = None
    for l,d in get_nums(infilepath) :
        if l == end_rows :
            if outfile is not None :
                outfile.close()
            outfilepath = f'{infilepath}_b.{blocks}.{curr_block}'
            outfile = open(outfilepath, 'wb')
            end_rows += min_rows
            curr_block += 1
        assert(l < end_rows)
        b = to_byte(d)
        outfile.write(b)
    outfile.close()
    assert(curr_block == blocks)

    #gzip
    if args.repr == 'gzip' :
        for outfilepath in [f'{infilepath}_b.{blocks}.{bid}' for bid in range(blocks)] :
            cmd = f'gzip -v {outfilepath}'
            execute_command_verbose(cmd, 1)
