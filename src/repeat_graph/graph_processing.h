//(c) 2016 by Authors
//This file is a part of ABruijn program.
//Released under the BSD license (see LICENSE file)

#pragma once

#include "repeat_graph.h"
#include "repeat_resolver.h"

struct UnbranchingPath
{
	UnbranchingPath(const GraphPath& path, 
					FastaRecord::Id id = FastaRecord::ID_NONE,
		   			bool circular = false, int length = 0, int meanCoverage = 0):
		   	path(path), id(id), circular(circular), repetitive(false), 
			length(length), meanCoverage(meanCoverage) {}

	std::string name() const
	{
		std::string nameTag = circular ? "circular" : "linear";
		return nameTag + "_" + std::to_string(id.signedId());
	}

	std::string nameUnsigned() const
	{
		std::string nameTag = circular ? "circular" : "linear";
		std::string idTag = id.strand() ? std::to_string(id.signedId()) : 
										  std::to_string(id.rc().signedId());
		return nameTag + "_" + idTag;
	}

	std::string edgesStr() const
	{
		if (path.empty()) return "";

		std::string contentsStr;
		for (auto& edge : path)
		{
			contentsStr += std::to_string(edge->edgeId.signedId()) + " -> ";
		}

		contentsStr.erase(contentsStr.size() - 4);
		return contentsStr;
	}

	GraphPath path;
	FastaRecord::Id id;
	std::string sequence;
	bool circular;
	bool repetitive;
	int length;
	int meanCoverage;
};


class GraphProcessor
{
public:
	GraphProcessor(RepeatGraph& graph, const SequenceContainer& asmSeqs,
				   const SequenceContainer& readSeqs):
		_graph(graph), _asmSeqs(asmSeqs), _readSeqs(readSeqs),
		_tipThreshold(Parameters::get().minimumOverlap) {}

	void condence();
	void fixChimericJunctions();
	std::vector<UnbranchingPath> getUnbranchingPaths();

private:
	void trimTips();
	void condenceEdges();
	void updateEdgesMultiplicity();
	void collapseBulges();

	RepeatGraph& _graph;
	const SequenceContainer& _asmSeqs;
	const SequenceContainer& _readSeqs;
	const int _tipThreshold;
};
