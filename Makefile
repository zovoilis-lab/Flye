ROOT_DIR := $(shell dirname $(realpath $(lastword $(MAKEFILE_LIST))))

export COMMON = -I${ROOT_DIR}/src/include
export LIBCUCKOO = -I${ROOT_DIR}/lib/libcuckoo
export INTERVAL_TREE = -I${ROOT_DIR}/lib/interval_tree
export LEMON = -I${ROOT_DIR}/lib/lemon
export BIN_DIR = ${ROOT_DIR}/bin
export MINIMAP2_DIR = ${ROOT_DIR}/lib/minimap2
export GRAPHMAP_DIR = ${ROOT_DIR}/lib/graphmap
export GRAPHMAP_BIN = ${GRAPHMAP_DIR}/bin/Linux-x64

export CXXFLAGS = ${LIBCUCKOO} ${INTERVAL_TREE} ${LEMON} ${COMMON} 
export LDFLAGS = -lz

.PHONY: clean all profile debug minimap2 graphmap

.DEFAULT_GOAL := all


${BIN_DIR}/flye-minimap2:
	make -C ${MINIMAP2_DIR}
	cp ${MINIMAP2_DIR}/minimap2 ${BIN_DIR}/flye-minimap2

${BIN_DIR}/flye-graphmap:
	make modules -C ${GRAPHMAP_DIR}
	make -C ${GRAPHMAP_DIR}
	cp ${GRAPHMAP_BIN}/graphmap ${BIN_DIR}/flye-graphmap

minimap2: ${BIN_DIR}/flye-minimap2

graphmap: ${BIN_DIR}/flye-graphmap

init:
	git submodule init
	git submodule update

all: init minimap2 graphmap
	make release -C src
profile: init minimap2 graphmap init
	make profile -C src
debug: init minimap2 graphmap init
	make debug -C src
clean:
	make clean -C src
	make clean -C ${MINIMAP2_DIR}
	make clean -C ${GRAPHMAP_DIR}	
	rm ${BIN_DIR}/flye-minimap2
	rm ${BIN_DIR}/flye-graphmap
