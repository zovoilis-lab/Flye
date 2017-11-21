//(c) 2016-2017 by Authors
//This file is a part of ABruijn program.
//Released under the BSD license (see LICENSE file)

#pragma once

#include "repeat_graph.h"
#include "read_aligner.h"
#include "graph_processing.h"

struct Contig
{
};

class ContigExtender
{
public:
	ContigExtender(RepeatGraph& graph, const ReadAligner& aligner,
				   const SequenceContainer& asmSeqs, 
				   const SequenceContainer& readSeqs):
		_graph(graph), _aligner(aligner), 
		_asmSeqs(asmSeqs), _readSeqs(readSeqs) {}

	void generateUnbranchingPaths();
	//void generatePathSequence();

	const std::vector<UnbranchingPath>& getUnbranchingPaths() 
		{return _unbranchingPaths;}
private:


	std::vector<UnbranchingPath> _unbranchingPaths;
	
	RepeatGraph& _graph;
	const ReadAligner& _aligner;
	const SequenceContainer& _asmSeqs;
	const SequenceContainer& _readSeqs;
};
