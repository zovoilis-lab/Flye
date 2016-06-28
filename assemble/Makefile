CXXFLAGS += -Wall -pthread -std=c++11
LDFLAGS += -pthread -std=c++11
OUT_BIN := ${BIN_DIR}/abruijn-assemble

.PHONY: all clean debug

profile: CXXFLAGS += -pg
profile: LDFLAGS += -pg
profile: all

all: CXXFLAGS += -D_LOG -O3 -DNDEBUG
all: abruijn

debug: CXXFLAGS += -g -Og -D_GLIBCXX_DEBUG -D_DEBUG -D_DEBUG_LOG -D_LOG
debug: abruijn

objects := ${patsubst %.cpp,%.o,${wildcard *.cpp}}
abruijn: ${objects}
	${CXX} ${objects} -o ${OUT_BIN} ${LDFLAGS}

%.o: %.cpp *.h
	${CXX} -c ${CXXFLAGS} $< -o $@

clean:
	rm *.o
	rm ${OUT_BIN}
