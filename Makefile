export COMMON = -I$(shell pwd)/src/include
export LIBCUCKOO = -I$(shell pwd)/lib/libcuckoo
export INTERVAL_TREE = -I$(shell pwd)/lib/interval_tree
export LEMON = -I$(shell pwd)/lib/lemon-1.3.1/install/include
export BIN_DIR = $(shell pwd)/bin

export CXXFLAGS = ${LIBCUCKOO} ${INTERVAL_TREE} ${LEMON} ${COMMON} 
#export LDFLAGS = 

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
