ASSEMBLE := $(shell pwd)/assemble
POLISH := $(shell pwd)/polishing
KNOT := $(shell pwd)/knot_graph

export LIBCUCKOO = $(shell pwd)/libcuckoo
export INCLUDE = $(shell pwd)/include
export BIN_DIR = $(shell pwd)/bin
export CXXFLAGS = -I${LIBCUCKOO} -I${INCLUDE}

.PHONY: clean all profile debug
.DEFAULT_GOAL := all

all: 
	make release -C ${ASSEMBLE}
	make release -C ${POLISH}
	make release -C ${KNOT}
profile:
	make profile -C ${ASSEMBLE}
	make release -C ${POLISH}
	make release -C ${KNOT}
debug:
	make debug -C ${ASSEMBLE}
	make debug -C ${POLISH}
	make debug -C ${KNOT}
clean: 
	make clean -C ${POLISH}
	make clean -C ${ASSEMBLE}
	make clean -C ${KNOT}
