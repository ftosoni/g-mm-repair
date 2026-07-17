# compilation flags
CXX_FLAGS=-std=c++17 -O3 -g
CFLAGS=-O3 -Wall -std=c99 -g
CC=gcc

# comment out this definition to get rid of malloc_count 
MALLOC_FLAGS=../tools/malloc_count.c -DMALLOC_COUNT -ldl

# executables in this directory
EXECS=remm recg csrmm csrcg ansremm decode

# malloc_count dedendencies for 
ifdef MALLOC_FLAGS
MALLOC_FILES=../tools/malloc_count.c ../tools/malloc_count.h
endif


# main executables 
EXECS=encode.x decode.x

# targets not producing a file declared phony
.PHONY: all clean

all: $(EXECS)

# test multiple version, requires boost
ansfold.x: ansfold.cpp cutil.hpp util.hpp methods.hpp
	$(CXX) $(CXX_FLAGS) -o $@ $< -lboost_system -lboost_regex -lboost_filesystem -lboost_program_options

encode.x: encode.cpp $(MALLOC_FILES)
	$(CXX) $(CXX_FLAGS) -o $@ $< $(MALLOC_FLAGS)

decode.x: decode.cpp decode.hpp $(MALLOC_FILES)
	$(CXX) $(CXX_FLAGS) -o $@ $< $(MALLOC_FLAGS)



%.o: %.c %.h
	$(CC) $(CFLAGS) -c -o $@ $<

clean:
	rm -f $(EXECS) *.o 
