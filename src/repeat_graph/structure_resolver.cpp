//(c) 2017 by Authors
//This file is a part of ABruijn program.
//Released under the BSD license (see LICENSE file)

#include "structure_resolver.h"

void StructureResolver::unrollLoops()
{
	
	std::unordered_set<GraphEdge*> toUnroll;
	for (GraphEdge* edge : _graph.iterEdges())
	{
		if (!edge->edgeId.strand() || edge->selfComplement ||
			!edge->isLooped()) continue;

		if (edge->nodeLeft->inEdges.size() == 2 &&
			edge->nodeLeft->outEdges.size() == 2 &&
			edge->nodeLeft->neighbors().size() == 2)
		{
			toUnroll.insert(edge);
		}
	}
	
	auto unroll = [this](GraphEdge* leftEdge, GraphEdge* loop, 
						 GraphEdge* rightEdge,
						 SequenceSegment segment,
						 std::vector<FastaRecord::Id> ids)
	{
		int mult = loop->multiplicity;
		_graph.removeEdge(loop);
		if (mult == 0) return;

		GraphNode* lastNode = leftEdge->nodeRight;
		for (int i = 0; i < mult; ++i)
		{
			GraphNode* newNode = _graph.addNode();
			GraphEdge* newEdge = _graph.addEdge(GraphEdge(lastNode, newNode,
													  	  ids[i]));
			newEdge->seqSegments.push_back(segment);
			lastNode = newNode;
		}
		vecRemove(rightEdge->nodeLeft->outEdges, rightEdge);
		rightEdge->nodeLeft = lastNode;
		lastNode->outEdges.push_back(rightEdge);
	};

	for (auto& edge : toUnroll)
	{
		GraphEdge* leftEdge = !edge->nodeLeft->inEdges[0]->isLooped() ? 
					edge->nodeLeft->inEdges[0] : edge->nodeLeft->inEdges[1];
		GraphEdge* rightEdge = !edge->nodeLeft->outEdges[0]->isLooped() ? 
					edge->nodeLeft->outEdges[0] : edge->nodeLeft->outEdges[1];
		if (leftEdge->edgeId == rightEdge->edgeId.rc() ||
			_graph.complementEdge(leftEdge) == rightEdge) continue;
		//if (leftEdge->repetitive || rightEdge->repetitive) continue;
		
		std::sort(edge->seqSegments.begin(), edge->seqSegments.end(),
				  [](const SequenceSegment& seg1, const SequenceSegment& seg2)
				  	{return seg1.length() < seg2.length();});
		auto segment = edge->seqSegments[edge->seqSegments.size() / 2];
		std::vector<FastaRecord::Id> fwdIds;
		std::vector<FastaRecord::Id> revIds;
		for (int i = 0; i < edge->multiplicity; ++i)
		{
			fwdIds.emplace_back(_graph.newEdgeId());
			revIds.emplace_back(fwdIds.back().rc());
		}
		std::reverse(revIds.begin(), revIds.end());

		Logger::get().debug() << "Unrolling " << edge->edgeId.signedId()
			<< " of multiplicity " << edge->multiplicity;

		GraphEdge* complLoop = _graph.complementEdge(edge);
		unroll(leftEdge, edge, rightEdge, segment, fwdIds);
		unroll(_graph.complementEdge(rightEdge), complLoop,
			   _graph.complementEdge(leftEdge), segment.complement(), revIds);
	}
}
