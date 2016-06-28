ASSEMBLE := $(shell pwd)/assemble
POLISH := $(shell pwd)/polishing

export LIBCUCKOO = $(shell pwd)/libcuckoo
export BIN_DIR = $(shell pwd)/bin
export CXXFLAGS = -I${LIBCUCKOO}

.PHONY: clean all

profile:
	make profile -C ${ASSEMBLE}
	make all -C ${POLISH}
all: 
	make all -C ${ASSEMBLE}
	make all -C ${POLISH}
debug:
	make debug -C ${ASSEMBLE}
	make debug -C ${POLISH}
clean: 
	make clean -C ${POLISH}
	make clean -C ${ASSEMBLE}
