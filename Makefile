ASSEMBLE := $(shell pwd)/assemble
POLISH := $(shell pwd)/polishing

export LIBCUCKOO = $(shell pwd)/libcuckoo
export INCLUDE = $(shell pwd)/include
export BIN_DIR = $(shell pwd)/bin
export CXXFLAGS = -I${LIBCUCKOO} -I${INCLUDE}

.PHONY: clean all profile debug
.DEFAULT_GOAL := all

all: 
	make release -C ${ASSEMBLE}
	make release -C ${POLISH}
profile:
	make profile -C ${ASSEMBLE}
	make release -C ${POLISH}
debug:
	make debug -C ${ASSEMBLE}
	make debug -C ${POLISH}
clean: 
	make clean -C ${POLISH}
	make clean -C ${ASSEMBLE}
