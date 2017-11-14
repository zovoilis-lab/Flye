//(c) 2016 by Authors
//This file is a part of ABruijn program.
//Released under the BSD license (see LICENSE file)

#pragma once

#include "repeat_graph.h"

//A simple class that assigns edges multiplicity based on the coverage
//and copmutes the mean coverage of all edges
class MultiplicityInferer
{
public:
	MultiplicityInferer(RepeatGraph& graph):
		_graph(graph), _uniqueCovThreshold(0), _meanCoverage(0) {}

	void estimateCoverage(const std::vector<GraphAlignment>& readAln);
	int  getMeanCoverage() const {return _meanCoverage;}

	//coverage threshold for an edge to be considered "unique"
	int  getUniqueCovThreshold() const 	{return _uniqueCovThreshold;}

private:

	RepeatGraph& _graph;
	int _uniqueCovThreshold; 
	int _meanCoverage;
};
