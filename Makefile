ASSEMBLE := $(shell pwd)/assemble
POLISH := $(shell pwd)/polishing

export LIBBF = $(shell pwd)/libbf
export BIN_DIR = $(shell pwd)/bin
export CXXFLAGS = -std=c++11 -I${LIBBF}
export LDFLAGS = -L${LIBBF} -lbf

.PHONY: clean all

all: 
	make all -C ${LIBBF}
	make all -C ${ASSEMBLE}
	make all -C ${POLISH}
debug:
	make debug -C ${LIBBF}
	make debug -C ${ASSEMBLE}
	make debug -C ${POLISH}
clean: 
	make clean -C ${LIBBF}
	make clean -C ${POLISH}
	make clean -C ${ASSEMBLE}
