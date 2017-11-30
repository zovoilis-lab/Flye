ROOT_DIR := $(shell dirname $(realpath $(lastword $(MAKEFILE_LIST))))

export COMMON = -I${ROOT_DIR}/src/include
export LIBCUCKOO = -I${ROOT_DIR}/lib/libcuckoo
export INTERVAL_TREE = -I${ROOT_DIR}/lib/interval_tree
export LEMON = -I${ROOT_DIR}/lib/lemon
export BIN_DIR = ${ROOT_DIR}/bin
export MINIMAP2_DIR = ${ROOT_DIR}/lib/minimap2

export CXXFLAGS = ${LIBCUCKOO} ${INTERVAL_TREE} ${LEMON} ${COMMON} 
export LDFLAGS = -lz

.PHONY: clean all profile debug minimap2

.DEFAULT_GOAL := all


${BIN_DIR}/abruijn-minimap2:
	make -C ${MINIMAP2_DIR}
	cp ${MINIMAP2_DIR}/minimap2 ${BIN_DIR}/abruijn-minimap2

minimap2: ${BIN_DIR}/abruijn-minimap2

all: minimap2
	make release -C src
profile: minimap2
	make profile -C src
debug: minimap2
	make debug -C src
clean:
	make clean -C src
	make clean -C ${MINIMAP2_DIR}
	rm ${BIN_DIR}/abruijn-minimap2
