//(c) 2016 by Authors
//This file is a part of ABruijn program.
//Released under the BSD license (see LICENSE file)

#pragma once

#include "sequence_container.h"
#include "overlap.h"

struct GluePoint
{
	GluePoint(size_t id = 0, FastaRecord::Id seqId = FastaRecord::ID_NONE,
			  int32_t position = 0):
		pointId(id), seqId(seqId), position(position) {}

	size_t 	pointId;
	FastaRecord::Id seqId;
	int32_t	position;
};


struct RepeatCluster
{
	RepeatCluster(FastaRecord::Id seqId = FastaRecord::ID_NONE,
				  size_t clusterId = 0, int32_t start = 0,
				  int32_t end = 0):
		seqId(seqId), clusterId(clusterId), start(start), end(end) {}

	FastaRecord::Id seqId;
	size_t clusterId;
	int32_t start;
	int32_t end;
};


struct GraphEdge
{
	GraphEdge(GluePoint gpLeft, GluePoint gpRight, bool rep, 
			  FastaRecord::Id edgeId, size_t clusterId):
		gpLeft(gpLeft), gpRight(gpRight), 
		repetitive(rep), edgeId(edgeId), clusterId(clusterId) {}

	GluePoint gpLeft;
	GluePoint gpRight;
	bool repetitive;

	FastaRecord::Id edgeId;
	size_t clusterId;
};


class RepeatGraph
{
public:
	RepeatGraph(const SequenceContainer& asmSeqs,
				const SequenceContainer& readSeqs):
		_asmSeqs(asmSeqs), _readSeqs(readSeqs)
	{}

	void build();
	void resolveRepeats();
	void outputDot(const std::string& filename, bool collapseRepeats);

private:
	void getRepeatClusters(const OverlapContainer& ovlps);
	void buildGraph(const OverlapContainer& ovlps);
	void initializeEdges();

	const int _maxSeparation = 1500;
	const SequenceContainer& _asmSeqs;
	const SequenceContainer& _readSeqs;

	std::unordered_map<FastaRecord::Id, 
					   std::vector<GluePoint>> _gluePoints;
	std::unordered_map<FastaRecord::Id, 
					   std::vector<GraphEdge>> _graphEdges;
	std::unordered_map<FastaRecord::Id, 
					   std::vector<RepeatCluster>> _repeatClusters;
	
};
