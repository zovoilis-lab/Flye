//(c) 2016-2017 by Authors
//This file is a part of ABruijn program.
//Released under the BSD license (see LICENSE file)

#pragma once

#include "repeat_graph.h"
#include "read_aligner.h"
#include "graph_processing.h"



class ContigExtender
{
public:
	ContigExtender(RepeatGraph& graph, const ReadAligner& aligner,
				   const SequenceContainer& asmSeqs, 
				   const SequenceContainer& readSeqs):
		_graph(graph), _aligner(aligner), 
		_asmSeqs(asmSeqs), _readSeqs(readSeqs) {}

	void generateUnbranchingPaths();
	void generateContigs();
	void outputContigs(const std::string& filename);
	//std::vector<UnbranchingPath> getContigPaths();

	const std::vector<UnbranchingPath>& getUnbranchingPaths() 
		{return _unbranchingPaths;}
private:
	struct Contig
	{
		Contig(const UnbranchingPath& corePath):
			graphEdges(corePath), graphPaths({&corePath})
		{}

		UnbranchingPath graphEdges;
		std::vector<const UnbranchingPath*> graphPaths;
		DnaSequence sequence;
	};

	std::vector<UnbranchingPath*> asUPaths(const GraphPath& path);

	std::vector<UnbranchingPath> _unbranchingPaths;
	std::unordered_map<GraphEdge*, UnbranchingPath*> _edgeToPath;
	std::vector<Contig> _contigs;
	
	RepeatGraph& _graph;
	const ReadAligner& _aligner;
	const SequenceContainer& _asmSeqs;
	const SequenceContainer& _readSeqs;
};
