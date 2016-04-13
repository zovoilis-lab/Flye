CXXFLAGS += -Wall
OUT_BIN := ${BIN_DIR}/polish

.PHONY: all clean debug

all: CXXFLAGS += -O3 -DNDEBUG
all: abruijn

debug: CXXFLAGS += -g -O0 -D_GLIBCXX_DEBUG -D_DEBUG -D_DEBUG_LOG -D_LOG
debug: abruijn

log: CXXFLAGS += -D_LOG -O2 -DNDEBUG
log: abruijn

objects := ${patsubst %.cpp,%.o,${wildcard *.cpp}}
abruijn: ${objects}
	${CXX} ${objects} -o ${OUT_BIN} ${LDFLAGS}

%.o: %.cpp *.h
	${CXX} -c ${CXXFLAGS} $< -o $@

clean:
	rm *.o
	rm ${OUT_BIN}
