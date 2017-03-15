//(c) 2016 by Authors
//This file is a part of ABruijn program.
//Released under the BSD license (see LICENSE file)

#pragma once

#include "repeat_graph.h"

class GraphProcessor
{
public:
	GraphProcessor(RepeatGraph& graph):
		_graph(graph) {}

	void simplify()
	{
		this->trimTips();
		this->unrollLoops();
		this->condenceEdges();
		this->updateEdgesMultiplicity();
	}

	void trimTips();
	void unrollLoops();
	void condenceEdges();
	void updateEdgesMultiplicity();

private:
	RepeatGraph& _graph;

	std::unordered_set<FastaRecord::Id> _outdatedEdges;
	const int _tipThreshold = 1500;
};
