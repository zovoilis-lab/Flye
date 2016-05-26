CXXFLAGS += -Wall -std=c++11
LDFLAGS += -std=c++11
OUT_LIB := ${LIBBF}/libbf.a

.PHONY: all clean

OBJECTS = $(patsubst %.cc, %.o, $(shell find . -name *.cc))

debug: CXXFLAGS += -g -Og -D_GLIBCXX_DEBUG -D_DEBUG
debug: bf

all: CXXFLAGS += -O3 -DNDEBUG
all: bf


bf: ${OBJECTS}
	ar rs ${OUT_LIB} ${OBJECTS}

%.o: %.cc *.h
	${CXX} -c ${CXXFLAGS} $< -o $@

clean:
	rm ${OBJECTS}
	rm ${OUT_LIB}
