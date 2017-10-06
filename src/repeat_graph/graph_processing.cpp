//(c) 2016 by Authors
//This file is a part of ABruijn program.
//Released under the BSD license (see LICENSE file)

#include <deque>

#include "graph_processing.h"
#include "../sequence/contig_generator.h"
#include "../common/logger.h"
#include "../common/config.h"
#include "../common/utils.h"

void GraphProcessor::condence()
{
	//this->trimTips();
	this->condenceEdges();
	this->fixChimericJunctions();
	this->collapseBulges();
	this->condenceEdges();
	this->trimTips();
}

void GraphProcessor::fixChimericJunctions()
{
	//a very specific case: 1 in - 1 out
	std::unordered_set<GraphNode*> simpleCases;
	for (auto& node : _graph.iterNodes())
	{
		if (!node->isBifurcation() &&
			(node->inEdges.front()->edgeId ==
			 node->outEdges.front()->edgeId.rc()))
		{
			simpleCases.insert(node);
		}
	}
	for (auto& node : simpleCases)
	{
		GraphNode* newNode = _graph.addNode();
		GraphEdge* cutEdge = node->outEdges.front();
		newNode->outEdges.push_back(cutEdge);
		cutEdge->nodeLeft = newNode;
		node->outEdges.clear();
	}

	//more common case: 2 in - 2 out
	std::unordered_set<GraphNode*> complexCases;
	for (auto& node : _graph.iterNodes())
	{
		if (node->inEdges.size() != 2 ||
			node->outEdges.size() != 2) continue;
		auto& inEdges = node->inEdges;
		auto& outEdges = node->outEdges;
		if (inEdges[0]->edgeId.rc() != outEdges[0]->edgeId)
		{
			//match INs with OUTs
			std::swap(inEdges[0], inEdges[1]);
		}

		if (inEdges[0]->edgeId.rc() == outEdges[0]->edgeId &&
			inEdges[1]->edgeId.rc() == outEdges[1]->edgeId)
		{
			complexCases.insert(node);
		}
	}
	for(auto& node : complexCases)
	{
		GraphNode* newNode = _graph.addNode();
		node->inEdges[1]->nodeRight = newNode;
		node->outEdges[0]->nodeLeft = newNode;
		newNode->inEdges.push_back(node->inEdges[1]);
		newNode->outEdges.push_back(node->outEdges[0]);

		node->inEdges.pop_back();
		node->outEdges.erase(node->outEdges.begin());
	}

	Logger::get().debug() << "Removed " 
		<< simpleCases.size() + complexCases.size()
		<< " chimeric junctions";
}

void GraphProcessor::collapseBulges()
{
	const int MAX_BUBBLE = Parameters::get().minimumOverlap;
	std::unordered_set<std::pair<GraphNode*, GraphNode*>,
					   pairhash> toFix;
	for (auto& edge : _graph.iterEdges())
	{
		//fixing simple bubbles like this: -<=>-
		if (edge->nodeLeft->outEdges.size()  != 2 ||
			edge->nodeLeft->inEdges.size()   != 1 ||
			edge->nodeRight->inEdges.size()  != 2 ||
			edge->nodeRight->outEdges.size() != 1) continue;

		GraphEdge* otherEdge = edge->nodeLeft->outEdges[0];
		if (otherEdge == edge) otherEdge = edge->nodeLeft->outEdges[1];
		if (otherEdge->nodeRight != edge->nodeRight) continue;
		if (edge->edgeId == otherEdge->edgeId.rc()) continue;

		if (edge->length() > MAX_BUBBLE || 
			otherEdge->length() > MAX_BUBBLE) continue;
		//if (edge->nodeLeft->inEdges.front()->length() < MAX_BUBBLE ||
		//	edge->nodeRight->outEdges.front()->length() < MAX_BUBBLE) continue;

		toFix.emplace(edge->nodeLeft, edge->nodeRight);
	}

	for (auto& nodes : toFix)
	{
		GraphEdge* edgeOne = nodes.first->outEdges[0];
		GraphEdge* edgeTwo = nodes.first->outEdges[1];
		if (abs(edgeOne->edgeId.signedId()) > abs(edgeTwo->edgeId.signedId()))
		{
			std::swap(edgeOne, edgeTwo);
		}
		for (auto& seg : edgeTwo->seqSegments)
		{
			edgeOne->seqSegments.push_back(seg);
		}
		_graph.removeEdge(edgeTwo);
	}
	Logger::get().debug() << "Collapsed " << toFix.size() / 2 << " bulges";
}

void GraphProcessor::trimTips()
{
	std::unordered_set<GraphNode*> toRemove;
	for (GraphEdge* tipEdge : _graph.iterEdges())
	{
		if (tipEdge->length() < _tipThreshold && tipEdge->isTip())
		{
			int leftDegree = tipEdge->nodeLeft->inEdges.size();
			toRemove.insert(leftDegree == 0 ? tipEdge->nodeLeft : 
							tipEdge->nodeRight);
		}
	}

	Logger::get().debug() << toRemove.size() << " tips removed";
	for (auto& edge : toRemove)	_graph.removeNode(edge);
}

void GraphProcessor::condenceEdges()
{
	int edgesRemoved = 0;
	int edgesAdded = 0;

	auto collapseEdges = [] (const GraphPath& edges)
	{
		std::vector<GraphEdge> newEdges;
		std::list<SequenceSegment> growingSeqs(edges.front()->seqSegments.begin(),
											   edges.front()->seqSegments.end());
		assert(edges.size() > 1);
		size_t prevStart = 0;
		for (size_t i = 1; i < edges.size(); ++i)
		{
			auto prevSeqs = growingSeqs;
			for (auto prevSeg = growingSeqs.begin(); 
				 prevSeg != growingSeqs.end(); )
			{
				bool continued = false;
				for (auto& nextSeg : edges[i]->seqSegments)
				{
					if (prevSeg->seqId == nextSeg.seqId &&
						prevSeg->end == nextSeg.start)
					{
						continued = true;
						prevSeg->end = nextSeg.end;
					}
				}
				if (!continued)
				{
					prevSeg = growingSeqs.erase(prevSeg);
				}
				else
				{
					++prevSeg;
				}
			}

			if (growingSeqs.empty())
			{
				newEdges.emplace_back(edges[prevStart]->nodeLeft, 
									  edges[i - 1]->nodeRight);
				std::copy(prevSeqs.begin(), prevSeqs.end(),
				  		  std::back_inserter(newEdges.back().seqSegments));

				std::copy(edges[i]->seqSegments.begin(), 
						  edges[i]->seqSegments.end(), 
						  std::back_inserter(growingSeqs));
				prevStart = i;
			}
		}

		newEdges.emplace_back(edges[prevStart]->nodeLeft, 
							  edges.back()->nodeRight);
		std::copy(growingSeqs.begin(), growingSeqs.end(),
				  std::back_inserter(newEdges.back().seqSegments));

		return newEdges;
	};

	std::vector<GraphPath> toCollapse;
	std::unordered_set<FastaRecord::Id> usedDirections;
	for (auto& node : _graph.iterNodes())
	{
		if (!node->isBifurcation()) continue;

		for (auto& direction : node->outEdges)
		{
			if (usedDirections.count(direction->edgeId)) continue;
			usedDirections.insert(direction->edgeId);

			GraphNode* curNode = direction->nodeRight;
			GraphPath traversed;
			traversed.push_back(direction);

			std::unordered_set<FastaRecord::Id> complementEdges;
			complementEdges.insert(direction->edgeId.rc());

			while (!curNode->isBifurcation() &&
				   !curNode->outEdges.empty() &&
				   !complementEdges.count(curNode->outEdges.front()->edgeId))
			{
				traversed.push_back(curNode->outEdges.front());
				complementEdges.insert(traversed.back()->edgeId.rc());
				curNode = curNode->outEdges.front()->nodeRight;
			}
			usedDirections.insert(traversed.back()->edgeId.rc());
			
			if (traversed.size() > 1)
			{
				toCollapse.emplace_back(std::move(traversed));
			}
		}
	}

	for (auto& unbranchingPath : toCollapse)
	{
		std::string collapsedStr;
		for (auto& edge : unbranchingPath)
		{
			collapsedStr += std::to_string(edge->edgeId.signedId()) + " -> ";
		}
		GraphPath complPath = _graph.complementPath(unbranchingPath);
		auto newEdges = collapseEdges(unbranchingPath);
		if (newEdges.size() == unbranchingPath.size()) continue;

		std::string addedStr;
		for (auto& edge : newEdges)
		{
			//do not add collapsed short loops
			//TODO: kinda hacky..
			//if (edge.length() < Parameters::get().minimumOverlap &&
			//	edge.isLooped()) continue;	

			GraphEdge addFwd = edge;
			addFwd.edgeId = _graph.newEdgeId();

			GraphEdge addRev(_graph.complementNode(edge.nodeRight),
							 _graph.complementNode(edge.nodeLeft),
							 addFwd.edgeId.rc());

			//complementary segments
			for (auto seqSeg : addFwd.seqSegments)
			{
				addRev.seqSegments.push_back(seqSeg.complement());
			}

			GraphEdge* addedEdge = _graph.addEdge(std::move(addFwd));
			_graph.addEdge(std::move(addRev));

			addedStr += std::to_string(addedEdge->edgeId.signedId()) + " -> ";
		}

		std::unordered_set<GraphEdge*> toRemove;
		for (auto& edge : unbranchingPath) toRemove.insert(edge);
		for (auto& edge : complPath) toRemove.insert(edge);
		for (auto& edge : toRemove) _graph.removeEdge(edge);

		edgesRemoved += unbranchingPath.size();
		edgesAdded += newEdges.size();

		//if (collapsedStr.size() > 4) collapsedStr.erase(collapsedStr.size() - 4);
		//if (addedStr.size() > 4) addedStr.erase(addedStr.size() - 4);
		collapsedStr.erase(collapsedStr.size() - 4);
		addedStr.erase(addedStr.size() - 4);
		Logger::get().debug() << "Collapsed: " << collapsedStr 
			<< " to " << addedStr;
	}

	Logger::get().debug() << "Removed " << edgesRemoved << " edges";
	Logger::get().debug() << "Added " << edgesAdded << " edges";
}
