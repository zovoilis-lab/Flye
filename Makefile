ROOT_DIR := $(shell dirname $(realpath $(lastword $(MAKEFILE_LIST))))

export LIBCUCKOO = -I${ROOT_DIR}/lib/libcuckoo
export INTERVAL_TREE = -I${ROOT_DIR}/lib/interval_tree
export LEMON = -I${ROOT_DIR}/lib/lemon
export BIN_DIR = ${ROOT_DIR}/bin
export MINIMAP2_DIR = ${ROOT_DIR}/lib/minimap2
export SAMTOOLS_DIR = ${ROOT_DIR}/lib/samtools-1.10

export CXXFLAGS += ${LIBCUCKOO} ${INTERVAL_TREE} ${LEMON} -I${MINIMAP2_DIR}
export LDFLAGS += -lz -L${MINIMAP2_DIR} -lminimap2

.PHONY: clean all profile debug minimap2 samtools samtools_config

.DEFAULT_GOAL := all


${BIN_DIR}/flye-minimap2:
	make -C ${MINIMAP2_DIR}
	cp ${MINIMAP2_DIR}/minimap2 ${BIN_DIR}/flye-minimap2

minimap2: ${BIN_DIR}/flye-minimap2

samtools: ${BIN_DIR}/flye-samtools

samtools_config:
	cd ${SAMTOOLS_DIR} && ./configure --without-curses --disable-bz2 --disable-lzma

${BIN_DIR}/flye-samtools: samtools_config
	make samtools -C ${SAMTOOLS_DIR}
	cp ${SAMTOOLS_DIR}/samtools ${BIN_DIR}/flye-samtools

all: minimap2 samtools
	make release -C src -j 4
profile: minimap2 samtools
	make profile -C src -j 4
debug: minimap2 samtools
	make debug -C src -j 4
clean:
	make clean -C src
	make clean -C ${MINIMAP2_DIR}
	make clean -C ${SAMTOOLS_DIR}
	rm ${BIN_DIR}/flye-minimap2
	rm ${BIN_DIR}/flye-samtools
