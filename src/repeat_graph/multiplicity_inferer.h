//(c) 2016 by Authors
//This file is a part of ABruijn program.
//Released under the BSD license (see LICENSE file)

#pragma once

#include "repeat_graph.h"
#include "read_aligner.h"

//A simple class that assigns edges multiplicity based on the coverage
//and copmutes the mean coverage of all edges
class MultiplicityInferer
{
public:
	MultiplicityInferer(RepeatGraph& graph, ReadAligner& aligner,
						const SequenceContainer& asmSeqs, 
						const SequenceContainer& readSeqs):
		_graph(graph), _aligner(aligner), _asmSeqs(asmSeqs), 
		_readSeqs(readSeqs), _uniqueCovThreshold(0), _meanCoverage(0) {}

	void estimateCoverage();
	int  getMeanCoverage() const {return _meanCoverage;}
	void removeUnsupportedEdges();
	void removeUnsupportedConnections();

	void collapseHeterozygousLoops();
	void collapseHeterozygousBulges();
	void trimTips();

	//coverage threshold for an edge to be considered "unique"
	int  getUniqueCovThreshold() const 	{return _uniqueCovThreshold;}

private:

	RepeatGraph& _graph;
	ReadAligner& _aligner;
	const SequenceContainer& _asmSeqs;
	const SequenceContainer& _readSeqs;
	int _uniqueCovThreshold; 
	int _meanCoverage;
};
