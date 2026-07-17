import sys

END = -1

Description = """
Herein we read teh list of pairs computed by the maximum weighed matching (MWM)
algorithm. The pairs are then concatenated, as to form a single column tour.
"""

#reading the pair list from file
def read_pairs(infilepath) :
    infile = open(infilepath, 'r')
    line = infile.readline()
    d = dict()
    while line :
        tokens = line[:-1].split(" ")
        assert(len(tokens) == 2)
        fst = int(tokens[0])
        snd = int(tokens[1])
        if snd>=0 :
            d[fst] = snd
        else :
            d[fst] = END
        line = infile.readline()

    infile.close()

    return d

#concatenate pairs ans determine the column permutation
def create_sequence(d) :
    ss = []
    for i in range(len(d)) :
        if d[i] is not None :
            last = i
            seq = [i]
            while d[last] != END :
                prev = last
                last = d[last]
                d[prev] = None
                seq.append(last)
            d[last] = None
            ss.append(seq)
    # print(ss)
    return ss

if __name__ == '__main__' :
    if len(sys.argv) != 2+1 :
        print('Usage is:',sys.argv[0],'<dataset> <generator ext.>')
        exit(-1)
    dataset = sys.argv[1]
    gen = sys.argv[2]

    name = '{}.{}'.format(dataset,gen)

    d = read_pairs(name+'.pairs')
    dimension = len(d)
    ss = create_sequence(d)

    #write the solution file
    solfile = open(name + '.mwm.solution', 'w')
    solfile.write('NAME : {name}.tour\n'.format(name=name))
    solfile.write('TYPE : TOUR\n')
    solfile.write('DIMENSION : {}\n'.format(str(dimension)))
    solfile.write('TOUR_SECTION\n')
    for s in ss :
        for e in s :
            solfile.write(str(e+1)) #1-based (instead of 0-based integers)
            solfile.write('\n')
    solfile.write('-1\nEOF\n')
    solfile.close()



