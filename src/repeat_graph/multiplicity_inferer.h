//(c) 2016 by Authors
//This file is a part of ABruijn program.
//Released under the BSD license (see LICENSE file)

#pragma once

#include "repeat_graph.h"
#include "repeat_resolver.h"


class MultiplicityInferer
{
public:
	MultiplicityInferer(RepeatGraph& graph):
		_graph(graph) {}

	void fixEdgesMultiplicity(const std::vector<GraphAlignment>& readAln);

private:
	void estimateByCoverage(const std::vector<GraphAlignment>& readAln);
	void balanceGraph();

	RepeatGraph& _graph;
};
