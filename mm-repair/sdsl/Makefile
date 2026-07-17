# Default position of the lib and include files for the SDSL-lite library; 
# change this definition if they are somewhere else
LIB_DIR = ${HOME}/lib
INC_DIR = ${HOME}/include

# Compilation flags and targets
CXX_FLAGS=-std=c++11 -O3 -Wall -I$(INC_DIR) -L$(LIB_DIR)
CCLIB=-lsdsl
SOURCES=$(wildcard *.cpp)
EXECS=$(SOURCES:.cpp=.x)

# malloc_count flags and dependencies
# comment next line to eclude malloc_count
# MALLOC_FLAGS=../tools/malloc_count.c -DMALLOC_COUNT -ldl
ifdef MALLOC_FLAGS
MALLOC_FILES=../tools/malloc_count.c ../tools/malloc_count.h
endif

all: $(EXECS)
	        
%.x: %.cpp $(MALLOC_FILES)
	$(CXX) $(CXX_FLAGS) -o $@ $< $(MALLOC_FLAGS) $(CCLIB) 

clean:
	rm -f $(EXECS)
	rm -rf *.dSYM

