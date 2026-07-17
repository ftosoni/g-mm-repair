# Testing multithreaded matrix multiplications
In this directory you’ll find the code for testing matrix-vector multiplications either gzip'ed and uncompressed matrices.

## Compile

1) Run `cmake -DCMAKE_BUILD_TYPE=Release .` or `cmake -DCMAKE_BUILD_TYPE=Debug .`
2) Run `make -j 6` (optionally specifying `VERBOSE=1`)

## Uncompressed representation

1) Given the covtype matrix in CSV format (581012 rows, 54 columns) we firstly obtain its block-wise dense representation with (assuming we want to use 8 threads)
```bash
python3 csv2repr.py <path to covtype>/covtype 581012 8 unco
```

2) As for matrix-to-vector multiplications, you need to have a vector, of course; you can, e.g., generate a vector containing 54 entries equal to `1.0` and store it to the file `x54.dbl` using the `makevec.py` tool in the parent directory.
```bash
makevec.py x54.dbl 54 1
```

3) Lastly, execute the `test_unco` program
```bash
./test_unco <path to covtype>/covtype <path to x54.dbl>/x54.dbl 581012 54 8
```

## Gzip'ed representation

1) As above, but the last parameter is `gzip` instead of `unco`.

2) As before, but you also need to gzip the input vector; for the above-mentioned example you should execute 
```bash
gzip <path to x54.dbl>/x54.dbl
```
thereby generating the file `x54.dbl.gz`.

3) Lastly, execute the `test_gzip` program
```bash
./test_gzip <path to covtype>/covtype <path to x54.dbl.gz>/x54.dbl.gz 581012 54 8
```
