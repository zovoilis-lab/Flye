LIBCPPLEX_DIR := $(shell pwd)/lib/cpplex

export COMMON = -I$(shell pwd)/src/include
export LIBCUCKOO = -I$(shell pwd)/lib/libcuckoo
export LIBCPPLEX = -I${LIBCPPLEX_DIR}/pilal/include -I${LIBCPPLEX_DIR}/simplex/include
export BIN_DIR = $(shell pwd)/bin

export CXXFLAGS = ${COMMON} ${LIBCUCKOO} ${LIBCPPLEX}
export LDFLAGS = -L${LIBCPPLEX_DIR} -lcpplex

.PHONY: clean all profile debug

.DEFAULT_GOAL := all

all: 
	make all -C ${LIBCPPLEX_DIR}
	make release -C src
profile:
	make all -C ${LIBCPPLEX_DIR}
	make profile -C src
debug:
	make all -C ${LIBCPPLEX_DIR}
	make debug -C src
clean: 
	make clean -C ${LIBCPPLEX_DIR}
	make clean -C src
