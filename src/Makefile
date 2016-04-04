CXXFLAGS += -Wall
OUT_BIN := ${BIN_DIR}/abruijn

.PHONY: all clean debug

all: CXXFLAGS += -O2 -DNDEBUG
all: abruijn

debug: CXXFLAGS += -g -D_GLIBCXX_DEBUG -D_DEBUG -D_DEBUG_LOG -D_LOG
debug: abruijn

log: CXXFLAGS += -D_LOG -O2 -DNDEBUG
log: abruijn

objects := ${patsubst %.cpp,%.o,${wildcard *.cpp}}
abruijn: ${objects}
	${CXX} ${LDFLAGS} ${objects} -o ${OUT_BIN}

%.o: %.cpp *.h
	${CXX} -c ${CXXFLAGS} $< -o $@

clean:
	rm *.o
	rm ${OUT_BIN}
