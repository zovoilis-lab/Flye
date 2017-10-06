//(c) 2016 by Authors
//This file is a part of ABruijn program.
//Released under the BSD license (see LICENSE file)

#pragma once

#include "repeat_graph.h"
#include "multiplicity_inferer.h"

class RepeatResolver
{
public:
	RepeatResolver(RepeatGraph& graph, const SequenceContainer& asmSeqs,
				   const SequenceContainer& readSeqs, 
				   const MultiplicityInferer& multInf): 
		_graph(graph), _asmSeqs(asmSeqs), _readSeqs(readSeqs), 
		_multInf(multInf) {}

	void alignReads();
	void findRepeats();
	void resolveRepeats();

	const std::vector<GraphAlignment>& getReadsAlignment() const
	{
		return _readAlignments;
	}

private:
	struct Connection
	{
		GraphPath path;
		SequenceSegment readSequence;
	};

	void clearResolvedRepeats();
	void removeUnsupportedEdges();
	std::vector<Connection> getConnections();
	int  resolveConnections(const std::vector<Connection>& conns);
	void separatePath(const GraphPath& path, SequenceSegment segment,
					  FastaRecord::Id startId);
	std::vector<GraphAlignment> 
		chainReadAlignments(const SequenceContainer& edgeSeqs,
							const std::vector<EdgeAlignment>& ovlps) const;
	void updateAlignments();

	std::vector<GraphAlignment> _readAlignments;

	RepeatGraph& _graph;
	const SequenceContainer&   _asmSeqs;
	const SequenceContainer&   _readSeqs;
	const MultiplicityInferer& _multInf;
};
