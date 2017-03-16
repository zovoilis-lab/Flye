ASSEMBLE := $(shell pwd)/src/assemble
POLISH := $(shell pwd)/src/polishing
REPEAT_GRAPH := $(shell pwd)/src/repeat_graph

export INCLUDE = $(shell pwd)/src/include
export LIBCUCKOO = $(shell pwd)/lib/libcuckoo
export BIN_DIR = $(shell pwd)/bin
export CXXFLAGS = -I${LIBCUCKOO} -I${INCLUDE}

.PHONY: clean all profile debug
.DEFAULT_GOAL := all

all: 
	make release -C src
profile:
	make profile -C src
debug:
	make debug -C src
clean: 
	make clean -C src
