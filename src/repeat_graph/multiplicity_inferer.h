//(c) 2016 by Authors
//This file is a part of ABruijn program.
//Released under the BSD license (see LICENSE file)

#pragma once

#include "repeat_graph.h"


class MultiplicityInferer
{
public:
	MultiplicityInferer(RepeatGraph& graph):
		_graph(graph) {}

	void fixEdgesMultiplicity();

private:
	void estimateByCoverage();
	void balanceGraph();

	RepeatGraph& _graph;
};
