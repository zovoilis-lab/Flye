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
	make release -C ${ASSEMBLE}
	make release -C ${POLISH}
profile:
	make profile -C src
	make profile -C ${POLISH}
	make profile -C ${ASSEMBLE}
debug:
	make debug -C src
	make debug -C ${ASSEMBLE}
	make debug -C ${POLISH}
clean: 
	make clean -C src
	make clean -C ${POLISH}
	make clean -C ${ASSEMBLE}
