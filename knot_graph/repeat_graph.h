//(c) 2016 by Authors
//This file is a part of ABruijn program.
//Released under the BSD license (see LICENSE file)

#pragma once

#include "sequence_container.h"
#include "overlap.h"

struct SequenceSegment
{
	SequenceSegment(FastaRecord::Id seqId, int32_t start, int32_t end):
		seqId(seqId), start(start), end(end) {}

	FastaRecord::Id seqId;
	int32_t start;
	int32_t end;
};

struct GraphNode;

struct GraphEdge
{
	GraphEdge(GraphNode* nodeLeft, GraphNode* nodeRight, 
			  FastaRecord::Id edgeId):
		nodeLeft(nodeLeft), nodeRight(nodeRight), 
		edgeId(edgeId), multiplicity(0)
		{}

	bool isRepetitive() {return multiplicity > 1;}
	void addSequence(FastaRecord::Id id, int32_t start, int32_t end)
	{
		seqSegments.emplace_back(id, start, end);
		++multiplicity;
	}
	int32_t length()
	{
		int32_t maxLen = std::numeric_limits<int32_t>::min();
		for (auto& seqSeg : seqSegments)
		{
			maxLen = std::max(maxLen, seqSeg.end - seqSeg.start);
		}
		return maxLen;
	}

	GraphNode* nodeLeft;
	GraphNode* nodeRight;

	FastaRecord::Id edgeId;
	std::vector<SequenceSegment> seqSegments;
	int multiplicity;
};

struct GraphNode
{
	bool isBifurcation() {return outEdges.size() > 1 || inEdges.size() > 1;}

	std::vector<GraphEdge*> inEdges;
	std::vector<GraphEdge*> outEdges;
};


class RepeatGraph
{
public:
	RepeatGraph(const SequenceContainer& asmSeqs,
				const SequenceContainer& readSeqs):
		_asmSeqs(asmSeqs), _readSeqs(readSeqs), _nextEdgeId(0)
	{}

	void build();
	void resolveRepeats();
	void outputDot(const std::string& filename, bool collapseRepeats);

private:
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

	struct Connection
	{
		GraphEdge* edgeIn;
		GraphEdge* edgeOut;
	};

	struct EdgeAlignment
	{
		OverlapRange overlap;
		GraphEdge* edge;
		SequenceSegment* segment;
	};

	void getGluepoints(const OverlapContainer& ovlps);
	void getRepeatClusters(const OverlapContainer& ovlps);
	void initializeEdges();
	void resolveConnections(const std::vector<Connection>& conns);
	std::vector<Connection> 
		chainReadAlignments(const SequenceContainer& edgeSeqs,
							std::vector<EdgeAlignment> ovlps);
	bool isRepetitive(GluePoint gpLeft, GluePoint gpRight);
	void fixTips();

	////////////////////////////////////////////////////////////////

	const int _maxSeparation = 1500;
	const SequenceContainer& _asmSeqs;
	const SequenceContainer& _readSeqs;


	std::unordered_map<FastaRecord::Id, 
					   std::vector<GluePoint>> _gluePoints;
	std::list<GraphNode> _graphNodes;
	std::list<GraphEdge> _graphEdges;

	//std::unordered_map<FastaRecord::Id, 
	//				   std::vector<GraphEdge>> _graphEdges;
	std::unordered_map<FastaRecord::Id, 
					   std::vector<RepeatCluster>> _repeatClusters;
	size_t _nextEdgeId;
};
