# Default position of the lib and include files for the SDSL-lite library; 
# change this definition if they are somewhere else
LIB_DIR = ${HOME}/lib
INC_DIR = ${HOME}/include

# Compilation flags
CFLAGS=-Wall -std=c99 -g -O3
CC=gcc 
CXX_FLAGS=-std=c++17 -g -O3 -march=native

# main executables in this directory
CONV_EXECS=bin2csrv bin2csrvf bin2csrvi bin2csv bin2csvf csvmat2csrv 
EXECS=csrvmm re32mm remm reivmm reans32mm $(CONV_EXECS) $(PR_EXECS)
# executables in pagerank directory
PR_EXECS=pagerank/repagerank pagerank/repagerank_old  pagerank/reivpagerank pagerank/re32pagerank pagerank/csrvpagerank pagerank/reans32pagerank

# comment out this definition to get rid of malloc_count 
#    MALLOC_FLAGS=tools/malloc_count.c -DMALLOC_COUNT -ldl
# malloc_count dependencies
ifdef MALLOC_FLAGS
MALLOC_FILES=tools/malloc_count.c tools/malloc_count.h
endif

# targets not producing a file declared phony
.PHONY: all brepair ans sdsl clean

# main target
all: $(EXECS) brepair ansf sdsl

# general rules for the targets in this directory
%: %.c
	$(CC) $(CFLAGS) -o $@ $< 
%: %.cpp
	$(CXX) $(CXX_FLAGS) -o $@ $< 

# special rules for the bin2csrv variants
bin2csrvf: bin2csrv.cpp
	$(CXX) $(CXX_FLAGS) -o $@ $< -DTypecode=2
bin2csrvi: bin2csrv.cpp
	$(CXX) $(CXX_FLAGS) -o $@ $< -DTypecode=1
bin2csvf: bin2csv.cpp
	$(CXX) $(CXX_FLAGS) -o $@ $< -DTypecode=2


# left and right matrix vector multiplication, double entries
# one could create versions for float and int defining FLOAT_VAL or INT_VALS 

# vc/val vectors only, no Repair
csrvmm: remm.c csrmatrix.h vector.h $(MALLOC_FILES)
	$(CC) $(CFLAGS) -o $@ $< $(MALLOC_FLAGS) -DCSR_MATRIX -pthread

# 32 bit vectors for R and C
re32mm: remm.c rematrix.h vector.h tools/xerrors.h $(MALLOC_FILES)
	$(CC) $(CFLAGS) -o $@ $< $(MALLOC_FLAGS) -pthread

# sdsl IV for R and C
reivmm: remm.c rematrix.hpp vector.h $(MALLOC_FILES)
	$(CXX) $(CXX_FLAGS) -o $@ $< $(MALLOC_FLAGS) -I$(INC_DIR) -L$(LIB_DIR) -lsdsl -pthread  -DUSE_INTVEC

# sdsl IV for R, ANS for C
remm: remm.c rematrix.hpp vector.h ans/decode.hpp $(MALLOC_FILES)
	$(CXX) $(CXX_FLAGS) -o $@ $< $(MALLOC_FLAGS) -I$(INC_DIR) -L$(LIB_DIR) -lsdsl -pthread -DUSE_ANSIV 

# not reported in the VLDB paper: not an interesting time/space tradeoff for the dataset
# 32 bit ints for R, ANS for C
reans32mm: remm.c rematrix.h vector.h ans/decode.hpp $(MALLOC_FILES)
	$(CXX) $(CXX_FLAGS) -o $@ $< $(MALLOC_FLAGS) -DUSE_ANS -pthread


# pagerank computation for different matrix formats
# see corresponding reXXmm executables for the format of the input matrix
pagerank/repagerank: pagerank/repagerank.c rematrix.hpp vector.h ans/decode.hpp $(MALLOC_FILES)
	$(CXX) $(CXX_FLAGS) -o $@ $< $(MALLOC_FLAGS) -I$(INC_DIR) -L$(LIB_DIR) -lsdsl -pthread -DUSE_ANSIV 

pagerank/repagerank_old: pagerank/repagerank_old.c rematrix.hpp vector.h ans/decode.hpp $(MALLOC_FILES)
	$(CXX) $(CXX_FLAGS) -o $@ $< $(MALLOC_FLAGS) -I$(INC_DIR) -L$(LIB_DIR) -lsdsl -pthread -DUSE_ANSIV 

pagerank/reivpagerank: pagerank/repagerank.c rematrix.hpp vector.h ans/decode.hpp $(MALLOC_FILES)
	$(CXX) $(CXX_FLAGS) -o $@ $< $(MALLOC_FLAGS) -I$(INC_DIR) -L$(LIB_DIR) -lsdsl -pthread -DUSE_INTVEC

pagerank/re32pagerank: pagerank/repagerank.c rematrix.hpp vector.h ans/decode.hpp $(MALLOC_FILES)
	$(CC) $(CFLAGS) -o $@ $< $(MALLOC_FLAGS) -I$(INC_DIR) -L$(LIB_DIR) -lsdsl -pthread

pagerank/reans32pagerank: pagerank/repagerank.c rematrix.hpp vector.h ans/decode.hpp $(MALLOC_FILES)
	$(CXX) $(CXX_FLAGS) -o $@ $< $(MALLOC_FLAGS) -DUSE_ANS -pthread

pagerank/csrvpagerank: pagerank/repagerank.c rematrix.hpp vector.h ans/decode.hpp $(MALLOC_FILES)
	$(CC) $(CFLAGS) -o $@ $< $(MALLOC_FLAGS) -I$(INC_DIR) -L$(LIB_DIR) -lsdsl -pthread -DCSR_MATRIX


# directory containing the (balanced) irepar/idespair tools 
brepair:
	make -C brepair

# directory containing the ans-fold encoding-decoding tools 
ansf:
	make -C ans

# directory containing the integer vector encoding-decoding tools 
sdsl:
	make -C sdsl


clean:
	rm -f $(EXECS)
	make -C brepair clean
	make -C ans clean
	make -C sdsl clean
