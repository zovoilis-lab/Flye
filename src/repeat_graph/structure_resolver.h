//(c) 2017 by Authors
//This file is a part of ABruijn program.
//Released under the BSD license (see LICENSE file)

#pragma once

#include "repeat_graph.h"
#include "read_aligner.h"

class StructureResolver
{
public:
	StructureResolver(RepeatGraph& graph, const SequenceContainer& asmSeqs,
				   	  const SequenceContainer& readSeqs, 
					  const ReadAligner& aligner):
		_graph(graph), _asmSeqs(asmSeqs), _readSeqs(readSeqs), 
		_aligner(aligner) {}

	//void unrollLoops();
	void scaffold();

private:
	//typedef std::pair<GraphEdge*, GraphEdge*> ScaffoldPair;
	struct ScaffoldInfo
	{
		GraphEdge* leftUnique;
		GraphEdge* rightUnique;
		std::unordered_set<GraphEdge*> repetitiveEdges;
	};

	std::vector<ScaffoldInfo> getScaffoldPairs();
	void untangle(ScaffoldInfo);

	RepeatGraph& _graph;
	const SequenceContainer& _asmSeqs;
	const SequenceContainer& _readSeqs;
	const ReadAligner& _aligner;
};


