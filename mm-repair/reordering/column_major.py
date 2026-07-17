import os, sys
import hashlib
import re 

Description = """
We read a matrix stored in row-major order, and we store it
using a sparse representation using a column-major order.
"""

#check whether a string represents an integer or not
def check_int(s):
    if s[0] in ('-', '+'):
        return s[1:].isdigit()
    return s.isdigit()

#hash a string into a 32-bit integer
def hash(s, separator='\.') :
    if check_int(s) :
        return s

    #check whether s represent 0
    regex = '^[+-]?([0-9]+)({separator}([0-9]+))?'.format(separator=separator)
    res = re.search(regex, s)
    assert(res), f"Called with [{s}], res={res}"
    
    is_zero = True
    for ch in res.group(1) :
        if not is_zero :
            break
        if ch != '0' :
            is_zero = False
    if res.group(3) is not None :
        for ch in  res.group(3) :
            if not is_zero :
                break
            if ch != '0' :
                is_zero = False

    if is_zero :
        #a zero remains a zero
        return 0 
    else :
        #hash
        b = bytes(s, 'utf-8')
        h = int.from_bytes(hashlib.sha256(b).digest()[:4], 'little') # 32-bit int
        return h


def count_cols(infilepath) :
    infile = open(infilepath, 'r')
    line = infile.readline()
    tokens = line.split(',')
    ncols = len(tokens)
    infile.close()
    return ncols

def _readlines(infile) :
    line = infile.readline()
    n = 0
    while line :
        yield n,line
        n +=1
        line = infile.readline()

if __name__ == '__main__' :
    if len(sys.argv) != 1+1 :
        print('Usage is:',sys.argv[0],'<infilepath>')
        exit(-1)
    
    infilename = sys.argv[1]

    infilepath = infilename
    dirpath = infilename+'_cols'
    os.system('if [ -d {d} ] ; then rm -rv {d} ; fi ; '.format(d=dirpath))
    os.mkdir(dirpath)
    ncols = count_cols(infilepath)
    outfiles = [open(path,'w') for path in ['{d}/{x}.txt'.format(d=dirpath,x=str(n)) for n in range(ncols)]]

    infile = open(infilepath, 'r')
    for i,line in _readlines(infile) :
        values = [hash(s.lstrip()) for s in line[:-1].split(',')]
        assert(len(values) == ncols)
        for v,f in zip(values,outfiles) :
            if v :
                s = '{}@{}\n'.format(str(v),str(i))
                f.write(s)
    for f in outfiles :
        f.flush()
        f.close()
    


