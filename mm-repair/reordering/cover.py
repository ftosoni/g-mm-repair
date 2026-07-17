import sys

Description = """
The PathCover approach to column reordering. 
"""

class Graph:
    def __init__(self, name):
        tspfile = open(name+'.tsp','r')
        line = tspfile.readline()
        while not line.startswith('DIMENSION : '):
            line = tspfile.readline()
        numv = line[len('DIMENSION : '):-1]

        self.graph = []
        self.name = name
        self.V = int(numv)

        #reading the similarity matrix
        line = tspfile.readline()
        while line != 'EDGE_WEIGHT_SECTION\n':
            line = tspfile.readline()
        for c1 in range(self.V-1) :
            line = tspfile.readline()[:-1]
            ws = [int(x) for x in line.split(' ') if x!='']
            for i,w in enumerate(ws) :
                c2 = c1+1+i
                self.add_edge(c1,c2,w)
        line = tspfile.readline()
        assert(line == 'EOF\n')
        tspfile.close()
 
    def add_edge(self, u, v, w):
        self.graph.append([u, v, w])

    def search(self, parent, i):
        if parent[i] == i:
            return i
        return self.search(parent, parent[i])
 
    def apply_union(self, parent, rank, x, y):
        xroot = self.search(parent, x)
        yroot = self.search(parent, y)
        if rank[xroot] < rank[yroot]:
            parent[xroot] = yroot
        elif rank[xroot] > rank[yroot]:
            parent[yroot] = xroot
        else:
            parent[yroot] = xroot
            rank[xroot] += 1
    
    """
    We implemented PathCover modifying the Kruskal's MST algorithm.
    Kruskal's algorithm computes the minimum spanning tree (MST) over a graph; cf. https://en.wikipedia.org/wiki/Kruskal%27s_algorithm
    """
    def kruskal(self) :
        total_cost = 0
        res = []
        parent = []
        rank = []
        for node in range(self.V):
            parent.append(node)
            rank.append(0)
        i,e = 0,0
        #sorting edges by decreasing weight
        self.graph = sorted(self.graph, key=lambda item: item[2])
        #every node in a (column) path has degree two at most
        neigh = [2] * self.V

        while e < self.V-1:
            #extract one edge
            u, v, w = self.graph[i]
            i = i + 1
            
            is_selected = False

            x = self.search(parent, u)
            y = self.search(parent, v)
            
            #we want our algorithm to produce an acyclic path.
            #for this reason we make sure that:
            # - we don't create loops
            # - every node in the spanning tree has degree at most 2.
            if x != y:
                if neigh[u] and neigh[v] :
                    neigh[u] -= 1
                    neigh[v] -= 1
                    is_selected = True

            if is_selected :
                e += 1
                self.apply_union(parent, rank, x, y)
                res.append((u,v))
                total_cost += w
        
        #at this point we have determined the the acyclic path representing the column reordering.
        #the path is represented by means of edges (i.e. node pairs) in 'res'.
        #we look for the beginning of the path, i.e., any of the two nodes with degree one, instead of two.
        head = None
        for i,v in enumerate(neigh) :
            if v==1 :
                head = i
                break
        assert(head is not None)

        #in the following we scan the set of selected edges and construct the column path.
        fwd,bwd = dict(),dict()
        for e1,e2 in res :
            if e1 not in fwd :
                fwd[e1] = set([e2])
            else :
                fwd[e1].add(e2)
            if e2 not in bwd :
                bwd[e2] = set([e1])
            else :
                bwd[e2].add(e1)
        
        seq = [head]
        empty_set = set()
        while  True:
            if head in fwd :
                new = fwd[head].pop()
                bwd[new].remove(head)
                if fwd[head] == empty_set :
                    del fwd[head]
                if bwd[new] == empty_set :
                    del bwd[new]
            elif head in bwd:
                new = bwd[head].pop()
                fwd[new].remove(head)
                if bwd[head] == empty_set :
                    del bwd[head]
                if fwd[new] == empty_set :
                    del fwd[new]
            else :
                assert(len(seq) == self.V)
                break
            seq.append(new)
            head = new
        
        #do we have concatenated all nodes (i.e. columns)?
        assert(len(seq) == self.V)

        #writing solution file
        solfile = open(self.name + '.cover.solution', 'w')
        tc = str(total_cost)
        solfile.write('NAME : {name}.{tc}.tour\n'.format(name=self.name,tc=tc))
        solfile.write('COMMENT : Length = {tc}\n'.format(tc=tc))
        solfile.write('TYPE : TOUR\n')
        solfile.write('DIMENSION : {}\n'.format(str(self.V)))
        solfile.write('TOUR_SECTION\n')
        for e in seq :
            solfile.write(str(e+1)) #1-based
            solfile.write('\n')
        solfile.write('-1\nEOF\n')
        solfile.close()


if __name__ == '__main__' :
    if len(sys.argv) != 1+1 :
        print('Usage is:',sys.argv[0],'<path to the tsp file>')
        exit(-1)
    
    path_to_tsp_file = sys.argv[1]
    tsp_filename = path_to_tsp_file.split('.tsp')[0]

    g = Graph(tsp_filename)
    g.kruskal()