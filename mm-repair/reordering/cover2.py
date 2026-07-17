import sys
import heapq

Description = """
The path PathCover+ approach to column reordering. 
"""

class Edge :
    def __init__(self, u, v, w, oldw=None) :
        self.u=u
        self.v=v 
        self.w=w
        self.oldw = w if oldw is None else oldw 

    def __lt__(self, other) :
        return self.w > other.w or (self.w==other.w and self.oldw>other.oldw)

    def __eq__(self, other) :
        return self.w == other.w and self.oldw == other.oldw

    def __gt__(self, other) :
        return self.w < other.w or (self.w==other.w and self.oldw<other.oldw)

    def __str__(self) :
        return '({},{}) {}'.format(self.u,self.v,self.w)

    def __repr__(self) :
        return self.__str__()

class Graph:
    def __init__(self, name):
        tspfile = open(name+'.tsp','r')
        line = tspfile.readline()
        while not line.startswith('DIMENSION : '):
            line = tspfile.readline()
        numv = line[len('DIMENSION : '):-1]

        self.name = name
        self.V = int(numv)
        self.adj = [dict() for _ in range(self.V)]
        self.Q = []
        self.Q0 = []

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

        heapq.heapify(self.Q)
 
    def add_edge(self, u, v, w):
        u,v = sorted((u,v))
        if w :
            self.adj[u][v] = w
            self.adj[v][u] = w
            e = (-w, -w, u, v, u, v)
            self.Q.append(e)
        else :
            self.Q0.append((u,v))

    def pop(self) :
        e = heapq.heappop(self.Q)
        return e

    def push(self, e) :
        heapq.heappush(self.Q , e )
        return

    def search(self, parent, i):
        if parent[i] == i:
            return i
        return self.search(parent, parent[i])
 
    def apply_union(self, parent, rank, x, y):
        xroot = self.search(parent, x)
        #assert(x == xroot)
        yroot = self.search(parent, y)
        #assert(y == yroot)
        if rank[xroot] < rank[yroot]:
            parent[xroot] = yroot
            return yroot,xroot
        elif rank[xroot] > rank[yroot]:
            parent[yroot] = xroot
            return xroot,yroot
        else:
            parent[yroot] = xroot
            rank[xroot] += 1
            return xroot,yroot

    def cover(self) :
        total_cost = 0
        res = []
        

        parent = [u for u in range(self.V)]
        rank = [0 for _ in range(self.V)]
        sels = [2] * self.V

        #We proceed considering edges by decreasing weight.
        #Each time we select an edge (u,v), we merge (u,v) into a new "supervertex" w.
        #All the edges that were incident to u and v are now incident to w.
        #The approach we follow is quite reminescent of the Sibeyn's algorithm
        #for the minimum spanning tree (MST) computation. For the Sibeyn's algorithm 
        #cf. e.g. Mehlhorn K., Sanders P., "Algorithms and Data Structures", chap. 11, Springer, 2008.

        #processing nodes with weight w>0
        while self.Q :
            _,w0,u,v,u0,v0 = self.pop()
            if u==v or sels[u0]==0 or sels[v0]==0 :
                continue #discard edge
            up, vp = self.search(parent,u), self.search(parent,v)
            if up==vp :
                continue #discard edge
            if up!=u or vp!=v :
                #update edge
                w = self.adj[up][vp]
                e = (-w,w0,up,vp,u0,v0)
                self.push(e)
            elif up != vp :
                adj_new = list(self.adj[u].items())
                adj_new.extend(self.adj[v].items())
                adj_new.sort(reverse=True) #heaviest first
                adj_new = {e[0]:e[1] for e in adj_new if (e[0] not in [up,vp])}
                relink_to,other = self.apply_union(parent,rank,up,vp)
                
                for n in self.adj[other] :
                    if n==relink_to :
                        continue
                    assert(other in self.adj[n])
                    w = adj_new[n]
                    assert(self.adj[other][n] >= w)
                    assert(self.adj[n][other] >= w)
                    del self.adj[n][other]
                    self.adj[n][relink_to] = w
                self.adj[other] = None #other removed
                self.adj[relink_to] = adj_new
                
                sels[u0] -= 1
                sels[v0] -= 1
                res.append((u0,v0))
        
        #processing nodes with weight w=0
        while self.Q0 :
            u0,v0 = self.Q0.pop()
            if sels[u0]==0 or sels[v0]==0 :
                continue
            u0p, v0p = self.search(parent,u0), self.search(parent,v0)
            if u0p != v0p :
                self.apply_union(parent,rank,u0p,v0p)
                sels[u0] -= 1
                sels[v0] -= 1
                res.append((u0,v0))

        #At this point we have computed the tour for all vertices.
        #We scan the set of selected edges and concatenate them as to form
        #an acyclic path.
        path_adj = [[] for _ in range(self.V)]
        for u,v in res :
            path_adj[u].append(v)
            path_adj[v].append(u)

        head,tail = None,None
        for i,l in enumerate(path_adj) :
            if len(l) == 1 :
                if head is None :
                    head = i
                else :
                    tail = i
            assert(l)
        assert(head is not None)
        assert(tail is not None)
        
        path = [head]
        prev,curr = None,head
        for _ in range(self.V-1) :
            next = path_adj[curr][0]
            if next==prev :
                next = path_adj[curr][1]
            path.append(next)
            prev,curr=path[-2],path[-1]
        
        assert(path[-1] == tail)

        #writing solution file
        solfile = open(self.name + '.cover2.solution', 'w')
        tc = str(total_cost)
        solfile.write('NAME : {name}.{tc}.tour\n'.format(name=self.name,tc=tc))
        solfile.write('COMMENT : Length = {tc}\n'.format(tc=tc))
        solfile.write('TYPE : TOUR\n')
        solfile.write('DIMENSION : {}\n'.format(str(self.V)))
        solfile.write('TOUR_SECTION\n')
        for e in path :
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
    g.cover()