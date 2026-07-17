# Matrix multiplication on RePair compressed matrices


## Prerequisites 

* Python version 3.8 or later

* [sdsl-lite](https://github.com/simongog/sdsl-lite) (`git clone https://github.com/simongog/sdsl-lite.git` + 
`cd sdsl-lite` + `./install.sh`)

* [psutil](https://pypi.org/project/psutil/) (`pip install psutil`)


## Installation 

Clone/download then `make`


## Sample computation

Given the [covtype](https://www.kaggle.com/giovannimanzini/some-machine-learning-matrices?select=covtype) matrix in csv format (581012 rows, 54 columns) we compress it with the command:
```bash
matrepair -r somedir/covtype 581012 54
```
that creates in the directory `somedir` the files `covtype.val`, `covtype.vc`, `covtype.vc.R`, `covtype.vc.C`, `covtype.vc.R.iv`, `covtype.vc.C.iv`, and `covtype.vc.C.ansf.1`, representing the output of different grammar compression algorithms (see [Available Matrix Compression Formats](Available-Matrix-Compression-Formats)). 

Next, we create a vector containing 54 entries equal to 1.0 and store it to the file `x54.dbl`:
```bash
makevec.py x54.dbl 54 1
```
Now we can compute the matrix vector products $`y=Mx`$ and $`z^T = y^T M`$ with the command:
```bash
remm  -y y.dbl -z z.dbl somedir/covtype 581012 54 x54.dbl
```
The output vector `y.dbl` has length 581012 and contains the sum of the entries of each row:
```bash
od -An -t f8 y.dbl | head
                    10300                    10077
                    13195                    13219
                     9963                     9709
                    10409                    10298
                    10456                    10378
                    10629                    13361
                    12989                    10640
                     9626                     9513
                    10310                     9479
                     9475                     9531
```
The output vector `z.dbl` has length 54 and contains the sum of the entries of each column of $`A^T A`$:
```bash
od -An -t f8 z.dbl | head
           14549643830975             757939056719
              65655800755            1372019643035
             227139887795           13332507438427
            1032604199658            1089965982517
             697072520994           11059033547178
               2526285260                225296661
               1924930316                178077606
                 12992451                 43226269
                 24825698                 78839086
                  7002929                 33807845
```


---

## Input matrix format

By default matrepair assumes the input matrix is in textual `csv` format. This behavior can be changed using the following command line options:

* `--bool`  the input matrix has only 0/1 entries. The matrix is represented by a text file in which each line contains a pair of row and column indices denoting the position of a nonzero elements. The pairs must be ordered in row-major order without duplicates

* `--sparse`  the input matrix is represented by a text file in which each line contains a triplet consisting of row and column indices, denoting the position of a nonzero, followed by the nonzero value. The triplets must be ordered in row-major order and without duplicates. 

* `--f64` the input matrix is represented in dense format using one float64 per entry. The total file size  is $`8 \cdot rows \cdot columns`'$ bytes

* `--f32` the input matrix is represented in dense format using one float32 per entry. The total file size  is $`4 \cdot rows \cdot columns`$ bytes


* `--i32` the input matrix is represented in dense format using one int32 per entry. The total file size  is $`4 \cdot rows \cdot columns`$ bytes

---

## Parallel computation


In order to exploits parallelism in matrix operations, we need to partition the matrix in blocks of rows and compress them separately. The command 
```bash
matrepair -r -b3 somedir/covtype 581012 54
```
splits the input matrix into three blocks of rows which are compressed separately. The command therefore creates the files `covtype.val`, `covtype.3.0.vc`, `covtype.3.0.vc.R`, `covtype.3.0.vc.C`, `covtype.3.0.vc.R.iv`, `covtype.3.0.vc.C.iv`, and `covtype.3.0.vc.C.ansf.1`, and similarly `covtype.3.1.vc` ...  `covtype.3.1.vc.C.ansf.1`, `covtype.3.2.vc` ... `covtype.3.2.vc.C.ansf.1`. 

Next, we can compute the matrix vector products $`y=Mx`$ and $`z^T = y^T M`$ in parallel using the above row blocks with the command:
```bash
remm -b3 -y y.dbl -z z.dbl covtype 581012 54 x54.dbl
```

---

## Available Matrix Compression Formats

### CSRV
Compressed Sparse Row/Value format. The encoding consists of the files with extensions `.[if]val` and `.vc` where `.[if]val` can be either `.ival`, `.fval`, or `.val` according to the type of the input matrix entries (respectively int32, float32 or float64)


The `.[if]val` file contains the *distinct* nonzero entries of the input matrix represented as 8 byte doubles or 4 byte float32 or int32. The `.vc` file contains for each nonzero element *a = A[i][j]* an encoding of the pair *(id,j)* where *j* is the column index and *id* is the index of the position in the `.val` file containing the actual value *a*. 


### Re32
The `.vc` file of the CSRV format is grammar compressed with RePair. All symbols are represented by 32 bit unsigned integers. The encoding consists of the files with extensions `.[if]val`, `.vc.R`, and `vc.C`


### ReIV
The `.vc.C` and `vc.R` files of the Re32 format are represented as packed arrays with the minimum number of bits per entry using the `int_vector` class from SDSL-lite.
The encoding consists of the files with extensions `.[if]val`, `.vc.R.iv`, and `vc.C.iv`


### ReANS
The `.vc.C` file of the re32 format is compressed using the ANS-fold-1 algorithm. The `.vc.R` file is represented as a packed array. This is the option usually providing the best compression.
The encoding consists of the files with extensions `.[if]val`, `.vc.R.iv`, and `vc.C.ansf.1`


---

## Tools

### matrepair
Tool to compute the CSRV representation of a matrix and to grammar-compress it. By default assumes the input matrix be in `csv` format; use one of the options `--i32`, `--f32` or `--f64` to specify that the input matrix is in binary format. If `--i32` or `--f32` are used the matrix entries are stored  as `int32` or `float32` in the value file which is accordingly named with the extension `.ival` or `.fval`.
Boolean (0/1) matrices are supported by the option `--bool`, which assumes that the input matrix consists of a text file containing a list of row/column pairs indicating the positions of the 1s. The list of pairs must be in row major order without duplicates, with one pair per row.  
The option `-r` shows a nice report detailing running times and compression ratios for the different formats Re32, ReIV and ReAns

### remm
Tool to compute a series of left and right matrix-vector multiplications reporting the overall running time and peak memory usage. Takes as input a compressed matrix in ReANS format (files `.val`, `.vc.R.iv`, `.vc.C.ansf.1`) of size *RxC*, a vector *x* of size *C* stored in a binary file (in float64 format) and an integer parameter *n*; computes *n* times the operations *y=Mx*, *z=y^t M*, *x = z/|z|*. You can use the tool `makevec.py` below to create an input vector *x* of the proper size.

### csrvmm
Analogous to *remm* (uses the same code) except that the input matrix is in Compressed Sparse Row Value (CSRV) format (just the files `.val` and `.vc`, no grammar compression) 

### reivmm, re32mm 
Analogous to *remm* except that the input matrix is expected to be in the format ReIV or Re32.

### makevec.py
Tool to create a vector of a given length and write it to a file in binary format (default float64 bit doubles, float32 and int32 formats are also supported). The vector is specified giving a set of values which are repeated cyclically.


---

## Bulk testing 

The tool *mmtest.py* can be used to test compression and (parallel) matrix-vector multiplication on a set of different matrices. The matrices, and their number of rows and columns, are specified inside *mmtest.py* in the global variables `Files` and `Sizes`. The first variable is a list of file names while the second is a dictionary providing the number of rows and columns for each file (extra entries in `Sizes` are ignored). The default content of the variables `Files` and `Sizes` can be overriden by the options `--files` and `-sizes`.
By default the input matrix is assumed to contains doubles in csv format,
the options `--i32`, `--f32`, and `--f64` signal that the input is in binary
format with entries of type `int32`, `float32` and `float64`. Type `mmtest.py -h` for usage help. 


The command
```bash 
mmtest.py mg -d /data
```
converts the input matrices from the `/data` directory from the *csv* format to the dense uncompressed format and applies to the latter the compressors *gzip* and *xz* showing the absolute size and the percentage with respect to the dense uncompressed matrix. Used for generating Table 1 in the paper. 


The command
```bash 
mmtest.py mz -b 8 -d /data
```
computes the CSRV and grammar representations of the input matrices from the `/data` directory and show their size as percentage of the dense uncompressed matrices. Before computing the CSRV representation the input matrices are split into 8 row blocks. Used for generating Table 1 in the paper.


The command
```bash 
mmtest.py mm -b 8 -d /data -n num
```
executes *num* iterations of the matrix multiplication algorithms *csrvmm*, *re32mm*, *reivmm* and *reansmm* showing the average time per iteration and the peak memory usage. The command assumes that the input matrices have been already split into 8 row blocks and compressed as described above. Used for generating Table 2 in the paper.


---

## PageRank computation 

As an application of the boolean matrix compressed format we have implemented the PageRank algorithm. See the `pagerank` subfolder.


---

## Column reordering

For the column reordering algorithms, take a look at the `reordering` subfolder. Be warned that they only support input in `csv` format.


---

## Experimental: the DRV format

The DRV (Dense Row Values) format is suitable for dense matrices since it treats 0 as any other matrix entry. Each matrix entry (including 0) is represented by an int32 id. This format will usually provide a better compression, but the matrix operations will have cost proportional to the size of the uncompressed matrix. At the moment matrix vector multiplication is not supported for the DRV format. To bulk test compression in the DRV format use `mmtest.py` with option `md`, for example:
```bash 
mmtest.py md -d /data
```
The purpose of this format is to explore the maximum compression achievable using Repair.


---

## Internal tools 


### csvmat2csrv
Tool to compute the CSRV representation of a matrix. The input matrix is assumed to consists of `float64` (doubles) stored in `csv` format (one row per line). 
Used by *matrepair*. Outputs the `.vc` and `.val` files. 


### bin2csrv, bin2csrvf, bin2csrvi
Tools to compute the CSRV representation of a matrix stored in binary form. 
The three versions assume that the matrix entries are stored respectively as float64 (double), float32, int32 and write such entries in the '[if]val' file in the same format. 
Used by *matrepair*. Outputs the `.vc` and `.[if]val` files. 




### brepair/irepair0
Tool using the RePair algorithm to grammar-compress a sequence of integers; the integer 0 is never compressed (ie, used in the rhs of a rule). Used by *matrepair* to compress the `.vc` file producing the `.vc.R` (rules) and `.vc.C` (sequence) files.

### sdsl/encode.x
Tool to encode a sequence of 32-bit integers as a sdsl integer vector using the minimum number of bits per entry from (https://github.com/simongog/sdsl-lite)[sdsl-lite]. Used by *matrepair* to generate the `.iv` files.

### ans/encode.x
Tool to encode a sequence of 32-bit integers using the *ANSfold-1* encoder from (https://github.com/mpetri/ans-large-alphabet)[ans-large-alphabet]. Used by *matrepair* to generate the `.ansf.1` files.



### others/csvmat2bin.py
Tool to convert a matrix in csv format into binary int32/float32/float64 format (possibly removing some trailing or leading rows/columns). All matrix entries are represented so the outfile has size `rows*cols*sizeof(entry)`. Note that when using the int32 or float32 output formats some information will be lost if the input values are not of the right type.


### others/stripcsvmat.py

Tool to strip the initial or final columns and/or the initial rows from a csv matrix, producing a smaller csv matrix. Used to remove unwanted data when preparing the test dataset.


### others/mat2csrv.py
Tool to compute the CSRV representation of a matrix. The input matrix is assumed to be in `csv` format unless its name ends with the `.dbl` extension in that case it is assumed to be in dense format with a 8 byte double per entry. Outputs the `.vc` and `.val` files. Superseeded by `csvmat2csrv` and  `bin2csrv[if]`.

### Citation

If you use this software or the techniques described herein in your research,
please cite the paper describing it. The reference is omitted here because this
repository is under double-blind review; it is cited in the accompanying
manuscript and will be restored in the camera-ready artifact.
