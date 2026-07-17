# PageRank computation for boolean RePair-compressed matrices


The tools in this directory can be used to compute the pagerank vector given the adjacency matrix of a graph. 
The supported input format is the [Matrix Market Exchange Formats](https://math.nist.gov/MatrixMarket/formats.html). The adjacency matrix of some web graphs in this format are available at [SuiteSparse Matrix collection](https://sparse.tamu.edu/LAW) or in a [Kaggle dataset](https://www.kaggle.com/datasets/wolfram77/graphs-law). The graphs are courtesy of the [Laboratory for Web Algorithmics](http://law.di.unimi.it/index.php), from the Univerity of Milano.


For techical reasons, to compute PageRank the adjaceny matrix of the web graph has to be: trasposed, self-loops (if any) removed, and the positions of the nonzero elements have to be sorted in row-major order. This is done by the `mtx2rowm` tool (which internally calls `egde2egde.py`). The `mtx2rowm` tool also subtracts one by all node ids (mtx format is 1 based) and writes to a binary file the number of nonzero elements in each column of the transformed matrix. 


For example, to compute the pagerank of the [cnr-2000.mtx](https://www.kaggle.com/datasets/wolfram77/graphs-law?select=cnr-2000.mtx) graph the procedure is the following. The command
```bash
mtx2rowm cnr-2000.mtx
```
generates the files `cnr-2000.mtx.rowm` and `cnr-2000.mtx.ccount`. The former contains the positions of the 3128710 nonzeros in row-major order, while the latter is a binary file of length 1302228 contaning the number of elements per column stored in 32-bit format; the graph has 325557 nodes, so the file has size $1302228 = 325557 \cdot 4$. `mtx2rowm` also generates the file `cnr-2000.mtx.header` which contains some human readable information on the matrix (it can be safely deleted). 

The output of `mtx2rowm` already provides the detailed instructions for computing the PageRank in a typical setting:
```
Matrix size: 325557
Number of nonzeros 3128710

The (transposed) matrix in row major order ready for
PageRank computation is in file cnr-2000.mtx.rowm
The column count file required by repagerank is cnr-2000.mtx.ccount

Assuming you are in the pagerank directory (adjust paths otherwise)
to compress the matrix with matrepair use the command line:
   ../matrepair --bool -r cnr-2000.mtx.rowm 325557 325557
Then, to compute PageRank use the command line
   repagerank -e 1e-4 -v cnr-2000.mtx.rowm cnr-2000.mtx.ccount
```
The command `../matrepair --bool -r cnr-2000.mtx.rowm 325557 325557` compresses the `cnr-2000.mtx.rowm` matrix. It generates different compressed formats; the most compact one is ReANS which consists of the files `cnr-2000.mtx.rowm.val`, `cnr-2000.mtx.rowm.vc.R.iv` and `cnr-2000.mtx.rowm.vc.C.ansf1` that together take 2717846 bytes, corresponding to roughly 6.95 bits per nonzero entry.  

Finally, the command line
```bash
repagerank -e 1e-4 -v cnr-2000.mtx.rowm cnr-2000.mtx.ccount
```
uses the compressed matrix and the file `cnr-2000.mtx.ccount` to compute PageRank and report the three nodes with the highest rank. In our example this should produce the following output:
```
==== Command line:
 repagerank -e 1e-4 -v cnr-2000.mtx.rowm 325557 cnr-2000.mtx.ccount
Number of nodes: 325557
Number of dandling nodes: 86959
Number of arcs: 3128710
Converged after 52 iterations
Sum of ranks: 1.000000 (should be 1)
Top 3 ranks:
  60597 0.025900
  60595 0.025900
  247028 0.005233
Top: 60597 60595 247028
Elapsed time: 0 secs
```
Type `repagerank` without arguments to get the available command line options. The executables `csrvpagerank` and `re32pagerank` and `reivpagerank` work as `repagerank` (they share the same code) but assume the compressed matrix is in format CSRV, Re32, and ReIV respectively.  
