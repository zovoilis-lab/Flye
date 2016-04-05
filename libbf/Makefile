CXXFLAGS += -Wall
OUT_LIB := ${LIBBF}/libbf.a

.PHONY: all clean

OBJECTS = $(patsubst %.cc, %.o, $(shell find . -name *.cc))

all: CXXFLAGS += -O3 -DNDEBUG
all: 
	echo ${OBJECTS}
all: bf


bf: ${OBJECTS}
	ar rs ${OUT_LIB} ${OBJECTS}

%.o: %.cc *.h
	${CXX} -c ${CXXFLAGS} $< -o $@

clean:
	rm ${OBJECTS}
	rm ${OUT_LIB}
