ROOT_DIR := $(shell dirname $(realpath $(lastword $(MAKEFILE_LIST))))

export COMMON = -I${ROOT_DIR}/src/include
export LIBCUCKOO = -I${ROOT_DIR}/lib/libcuckoo
export INTERVAL_TREE = -I${ROOT_DIR}/lib/interval_tree
export LEMON = -I${ROOT_DIR}/lib/lemon
export BIN_DIR = ${ROOT_DIR}/bin
export MINIMAP2_DIR = ${ROOT_DIR}/lib/minimap2
export GRAPHMAP_DIR = ${ROOT_DIR}/lib/graphmap/bin/Linux-x64/

export CXXFLAGS = ${LIBCUCKOO} ${INTERVAL_TREE} ${LEMON} ${COMMON} 
#export LDFLAGS = 

.PHONY: clean all profile debug minimap2 graphmap

.DEFAULT_GOAL := all


${BIN_DIR}/abruijn-minimap2:
	make -C ${MINIMAP2_DIR}
	cp ${MINIMAP2_DIR}/minimap2 ${BIN_DIR}/abruijn-minimap2
	

${BIN_DIR}/abruijn-graphmap:
	make modules -C ${GRAPHMAP_DIR}
	make -C ${GRAPHMAP_DIR}
	cp ${GRAPHMAP_DIR}/graphmap ${BIN_DIR}/abruijn-graphmap

minimap2: ${BIN_DIR}/abruijn-minimap2

graphmap: ${BIN_DIR}/abruijn-graphmap

all: minimap2 graphmap
	make release -C src
profile: minimap2 graphmap
	make profile -C src
debug: minimap2 graphmap
	make debug -C src
clean:
	make clean -C src
	make clean -C ${MINIMAP2_DIR}
	make clean -C ${GRAPHMAP_DIR}	
	rm ${BIN_DIR}/abruijn-minimap2
	rm ${BIN_DIR}/abruijn-graphmap
