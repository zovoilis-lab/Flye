//(c) 2016 by Authors
//This file is a part of ABruijn program.
//Released under the BSD license (see LICENSE file)

#pragma once

#include "../sequence/sequence_container.h"
#include "../sequence/overlap.h"
#include "config.h"
#include "utils.h"
#include <list>

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
		edgeId(edgeId), multiplicity(0), selfComplement(false), 
		readSequence(false)
		{}

	bool isRepetitive() {return multiplicity > 1;}
	void addSequence(FastaRecord::Id id, int32_t start, int32_t end)
	{
		seqSegments.emplace_back(id, start, end);
		++multiplicity;
	}
	int32_t length()
	{
		if (seqSegments.empty()) return 0;

		int64_t sumLen = 0;
		for (auto& seqSeg : seqSegments)
		{
			sumLen += seqSeg.end - seqSeg.start;
		}
		return sumLen / seqSegments.size();
	}

	GraphNode* nodeLeft;
	GraphNode* nodeRight;

	FastaRecord::Id edgeId;
	std::vector<SequenceSegment> seqSegments;
	int multiplicity;
	bool selfComplement;
	bool readSequence;
};

struct GraphNode
{
	bool isBifurcation() {return outEdges.size() != 1 || inEdges.size() != 1;}
	std::vector<GraphNode*> neighbors()
	{
		std::unordered_set<GraphNode*> result;
		for (auto& edge : inEdges) result.insert(edge->nodeLeft);
		for (auto& edge : outEdges) result.insert(edge->nodeRight);

		return std::vector<GraphNode*>(result.begin(), result.end());
	}

	std::vector<GraphEdge*> inEdges;
	std::vector<GraphEdge*> outEdges;
};

typedef std::vector<GraphEdge*> GraphPath;

class RepeatGraph
{
public:
	RepeatGraph(const SequenceContainer& asmSeqs):
		 _nextEdgeId(0), _asmSeqs(asmSeqs)
	{}

	void build();
	void outputDot(const std::string& filename);
	GraphPath complementPath(const GraphPath& path);

	//nodes
	GraphNode* addNode()
	{
		GraphNode* node = new GraphNode();
		_graphNodes.insert(node);
		return node;
	}

	class IterNodes
	{
	public:
		IterNodes(RepeatGraph& graph): _graph(graph) {}

		std::unordered_set<GraphNode*>::iterator begin() 
			{return _graph._graphNodes.begin();}
		std::unordered_set<GraphNode*>::iterator end() 
			{return _graph._graphNodes.end();}
	
	private:
		RepeatGraph& _graph;
	};
	IterNodes iterNodes() {return IterNodes(*this);}
	//
	
	//edges
	GraphEdge* addEdge(GraphEdge&& edge)
	{
		GraphEdge* newEdge = new GraphEdge(edge);
		_graphEdges.insert(newEdge);
		newEdge->nodeLeft->outEdges.push_back(newEdge);
		newEdge->nodeRight->inEdges.push_back(newEdge);
		return newEdge;
	}
	class IterEdges
	{
	public:
		IterEdges(RepeatGraph& graph): _graph(graph) {}

		std::unordered_set<GraphEdge*>::iterator begin() 
			{return _graph._graphEdges.begin();}
		std::unordered_set<GraphEdge*>::iterator end() 
			{return _graph._graphEdges.end();}
	
	private:
		RepeatGraph& _graph;
	};
	IterEdges iterEdges() {return IterEdges(*this);}
	void removeEdge(GraphEdge* edge)
	{
		vecRemove(edge->nodeRight->inEdges, edge);
		vecRemove(edge->nodeLeft->outEdges, edge);
		_graphEdges.erase(edge);
		delete edge;
	}
	//

	size_t _nextEdgeId;	//TODO: temporary

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

	void getGluepoints(const OverlapContainer& ovlps);
	void getRepeatClusters(const OverlapContainer& ovlps);
	void initializeEdges();
	bool isRepetitive(GluePoint gpLeft, GluePoint gpRight);
	
	const SequenceContainer& _asmSeqs;
	const int _maxSeparation = Constants::maxSeparation;

	std::unordered_map<FastaRecord::Id, 
					   std::vector<GluePoint>> _gluePoints;
	std::unordered_map<FastaRecord::Id, 
					   std::vector<RepeatCluster>> _repeatClusters;

	std::unordered_set<GraphNode*> _graphNodes;
	std::unordered_set<GraphEdge*> _graphEdges;
};
