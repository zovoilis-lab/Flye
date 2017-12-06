//(c) 2016 by Authors
//This file is a part of ABruijn program.
//Released under the BSD license (see LICENSE file)

#pragma once

#include <list>

#include "../sequence/sequence_container.h"
#include "../sequence/overlap.h"
#include "../common/config.h"
#include "../common/utils.h"

struct SequenceSegment
{
	SequenceSegment(FastaRecord::Id seqId = FastaRecord::ID_NONE, 
					int32_t seqLen = 0, int32_t start = 0, 
					int32_t end = 0):
		seqId(seqId), seqLen(seqLen), start(start), 
		end(end), readSequence(false) {}

	SequenceSegment complement() const
	{
		SequenceSegment other(*this);
		other.seqId = seqId.rc();
		other.start = seqLen - end - 1;
		other.end = seqLen - start - 1;
		return other;
	}

	int32_t length() const {return end - start;}

	bool operator==(const SequenceSegment& other)
	{
		return seqId == other.seqId && start == other.start && end == other.end;
	}

	FastaRecord::Id seqId;
	int32_t seqLen;
	int32_t start;
	int32_t end;
	bool readSequence;
};

struct GraphNode;

struct GraphEdge
{
	GraphEdge(GraphNode* nodeLeft, GraphNode* nodeRight, 
			  FastaRecord::Id edgeId = FastaRecord::ID_NONE):
		nodeLeft(nodeLeft), nodeRight(nodeRight), 
		edgeId(edgeId), multiplicity(0), repetitive(false), 
		selfComplement(false), resolved(false), meanCoverage(0) {}

	bool isRepetitive() const 
		{return repetitive;}

	bool isLooped() const 
		{return nodeLeft == nodeRight;}

	bool isTip() const;

	void addSequence(FastaRecord::Id id, int32_t length, 
					 int32_t start, int32_t end)
	{
		seqSegments.emplace_back(id, length, start, end);
		++multiplicity;
	}

	int32_t length() const
	{
		if (seqSegments.empty()) return 0;

		int64_t sumLen = 0;
		for (auto& seqSeg : seqSegments)
		{
			sumLen += seqSeg.end - seqSeg.start;
		}
		return sumLen / seqSegments.size();
	}

	std::unordered_set<GraphEdge*> adjacentEdges();

	GraphNode* nodeLeft;
	GraphNode* nodeRight;

	FastaRecord::Id edgeId;
	std::vector<SequenceSegment> seqSegments;

	int  multiplicity;
	bool repetitive;
	bool selfComplement;
	bool resolved;
	int  meanCoverage;
};

struct GraphNode
{
	bool isBifurcation() const
		{return outEdges.size() != 1 || inEdges.size() != 1;}

	std::unordered_set<GraphNode*> neighbors() const
	{
		std::unordered_set<GraphNode*> result;
		for (auto& edge : inEdges) 
		{
			if (edge->nodeLeft != this) result.insert(edge->nodeLeft);
		}
		for (auto& edge : outEdges) 
		{
			if (edge->nodeRight != this) result.insert(edge->nodeRight);
		}

		return result;
	}

	bool isEnd() const
	{
		int inDegree = 0;
		for (auto& edge : inEdges)
		{
			if (!edge->isLooped()) ++inDegree;
		}
		int outDegree = 0;
		for (auto& edge : outEdges)
		{
			if (!edge->isLooped()) ++outDegree;
		}
		return (inDegree == 1 && outDegree == 0) || 
			   (inDegree == 0 && outDegree == 1);
	}

	bool isTelomere() const
	{
		int numIn = 0;
		int numOut = 0;
		for (auto& edge: inEdges)
		{
			if (!edge->isLooped()) ++numIn;
		}
		for (auto& edge: outEdges)
		{
			if (!edge->isLooped()) ++numOut;
		}
		if ((bool)numIn != (bool)numOut)
		{
			return true;
		}
		return false;
	}


	bool isResolved() const
	{
		int inDegree = 0;
		for (auto& edge : inEdges)
		{
			if (!edge->isLooped()) ++inDegree;
		}
		int outDegree = 0;
		for (auto& edge : outEdges)
		{
			if (!edge->isLooped()) ++outDegree;
		}
		return inDegree == 1 && outDegree == 1;
	}

	std::vector<GraphEdge*> inEdges;
	std::vector<GraphEdge*> outEdges;
};

typedef std::vector<GraphEdge*> GraphPath;

struct EdgeAlignment
{
	OverlapRange overlap;
	GraphEdge* edge;
	SequenceSegment segment;
};
typedef std::vector<EdgeAlignment> GraphAlignment;

class RepeatGraph
{
public:
	RepeatGraph(const SequenceContainer& asmSeqs):
		 _nextEdgeId(0), _asmSeqs(asmSeqs)
	{}

	void build();
	GraphPath  complementPath(const GraphPath& path) const;
	GraphEdge* complementEdge(GraphEdge* edge) const;
	GraphNode* complementNode(GraphNode* node) const;

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
		
		_idToEdge[newEdge->edgeId] = newEdge;
		if (newEdge->selfComplement)
		{
			_idToEdge[newEdge->edgeId.rc()] = newEdge;
		}
		return newEdge;
	}
	bool hasEdge(GraphEdge* edge)
	{
		return _graphEdges.count(edge);
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
		_idToEdge.erase(edge->edgeId);
		delete edge;
	}
	void removeNode(GraphNode* node)
	{
		std::unordered_set<GraphEdge*> toRemove;
		for (auto& edge : node->outEdges) 
		{
			vecRemove(edge->nodeRight->inEdges, edge);
			toRemove.insert(edge);
		}
		for (auto& edge : node->inEdges) 
		{
			vecRemove(edge->nodeLeft->outEdges, edge);
			toRemove.insert(edge);
		}
		for (auto& edge : toRemove)
		{
			_graphEdges.erase(edge);
			delete edge;
		}
		_graphNodes.erase(node);
		delete node;
	}
	//
	FastaRecord::Id newEdgeId()
	{
		size_t curId = _nextEdgeId;
		_nextEdgeId += 2;
		return FastaRecord::Id(curId);
	}

private:
	size_t _nextEdgeId;

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
	void initializeEdges(const OverlapContainer& asmOverlaps);
	void collapseTandems();
	void logEdges();
	
	const SequenceContainer& _asmSeqs;
	const int _maxSeparation = Config::get("max_separation");

	std::unordered_map<FastaRecord::Id, 
					   std::vector<GluePoint>> _gluePoints;

	std::unordered_set<GraphNode*> _graphNodes;
	std::unordered_set<GraphEdge*> _graphEdges;
	std::unordered_map<FastaRecord::Id, GraphEdge*> _idToEdge;

};
