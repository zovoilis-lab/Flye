.PHONY: all clean debug profile

CXXFLAGS += -Wall -pthread -std=c++11
LDFLAGS += -pthread -std=c++11


profile: CXXFLAGS += -pg
profile: LDFLAGS += -pg
profile: release

release: CXXFLAGS += -D_LOG -O3 -DNDEBUG
release: abruijn-repeat

debug: CXXFLAGS += -g -O0 -D_GLIBCXX_DEBUG -D_DEBUG -D_DEBUG_LOG -D_LOG
debug: abruijn-repeat


#sequence module
sequence_obj := ${patsubst %.cpp,%.o,${wildcard sequence/*.cpp}}

sequence/%.o: sequence/%.cpp sequence/*.h
	${CXX} -c ${CXXFLAGS} $< -o $@

#abruijn-repeat module
REPEAT_BIN := ${BIN_DIR}/abruijn-repeat
repeat_obj := ${patsubst %.cpp,%.o,${wildcard repeat_graph/*.cpp}}
abruijn-repeat: ${repeat_obj} ${sequence_obj}
	${CXX} ${repeat_obj} ${sequence_obj} -o ${REPEAT_BIN} ${LDFLAGS}

repeat_graph/%.o: repeat_graph/%.cpp repeat_graph/*.h sequence/*.h
	${CXX} -c ${CXXFLAGS} $< -o $@


clean:
	rm ${REPEAT_BIN}
	rm ${repeat_obj}
	rm ${sequence_obj}
