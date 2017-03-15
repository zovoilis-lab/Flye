//(c) 2016 by Authors
//This file is a part of ABruijn program.
//Released under the BSD license (see LICENSE file)

#pragma once

#include "repeat_graph.h"

class RepeatResolver
{
public:
	RepeatResolver(RepeatGraph& graph, const SequenceContainer& asmSeqs,
				   const SequenceContainer& readSeqs): 
		_graph(graph), _asmSeqs(asmSeqs), _readSeqs(readSeqs) {}

	void resolveRepeats();

private:
	struct EdgeAlignment
	{
		OverlapRange overlap;
		GraphEdge* edge;
		SequenceSegment* segment;
	};

	struct Connection
	{
		GraphPath path;
		int32_t leftOverlap;
		int32_t rightOverlap;
	};

	void resolveConnections(const std::vector<Connection>& conns);
	size_t separatePath(const GraphPath& path, size_t startId);
	std::vector<EdgeAlignment> 
		chainReadAlignments(const SequenceContainer& edgeSeqs,
							std::vector<EdgeAlignment> ovlps);

	const int _readJump = 1500;
	const int _readOverhang = 500;

	RepeatGraph& _graph;

	const SequenceContainer& _asmSeqs;
	const SequenceContainer& _readSeqs;
};
