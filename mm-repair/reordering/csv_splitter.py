import sys

def iter_lines(infilename) :
    infile = open(infilename, 'r')
    line = infile.readline()
    while line :
        yield line
        line = infile.readline()
    infile.close()

def split_csv(infilename, m, b) :
    infile = open(infilename, 'r')
    maxr = (m + b - 1) // b

    #first b-1 blocks
    for i in range(b-1) :
        outfilename = '{ifn}.{b}.{i}'.format(ifn=infilename, b=str(b), i=str(i))
        outfile = open(outfilename, 'w')
        for _ in range(maxr) :
            line = infile.readline()
            outfile.write(line)
        outfile.close()
    
    #last block
    outfilename = '{ifn}.{b}.{i}'.format(ifn=infilename, b=str(b), i=str(b-1))
    outfile = open(outfilename, 'w')
    line = infile.readline()
    while line :
        outfile.write(line)
        line = infile.readline()
    outfile.close()

if __name__ == '__main__' :
    if len(sys.argv) != 3 + 1 :
        print('Usage is:', sys.argv[0], '<data set> <rows> <num blocks>')
        exit(-1)
    
    infilename = sys.argv[1]
    m = int(sys.argv[2])
    b = int(sys.argv[3])

    split_csv(infilename,m,b)
