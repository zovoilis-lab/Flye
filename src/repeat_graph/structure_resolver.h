//(c) 2017 by Authors
//This file is a part of ABruijn program.
//Released under the BSD license (see LICENSE file)

#pragma once

#include "repeat_graph.h"


class StructureResolver
{
public:
	StructureResolver(RepeatGraph& graph, const SequenceContainer& asmSeqs,
				   	  const SequenceContainer& readSeqs):
		_graph(graph), _asmSeqs(asmSeqs), _readSeqs(readSeqs) {}

	void unrollLoops();

private:
	RepeatGraph& _graph;
	const SequenceContainer& _asmSeqs;
	const SequenceContainer& _readSeqs;
};


