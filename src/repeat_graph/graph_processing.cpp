//(c) 2016 by Authors
//This file is a part of ABruijn program.
//Released under the BSD license (see LICENSE file)

#include <deque>

#include "graph_processing.h"
#include "logger.h"
#include "utils.h"


void GraphProcessor::unrollLoops()
{
	auto unrollEdge = [](GraphEdge& loopEdge)
	{
		GraphEdge* prevEdge = nullptr;
		for (auto& inEdge : loopEdge.nodeLeft->inEdges)
		{
			if (inEdge != &loopEdge) prevEdge = inEdge;
		}
		GraphEdge* nextEdge = nullptr;
		for (auto& outEdge : loopEdge.nodeLeft->outEdges)
		{
			if (outEdge != &loopEdge) nextEdge = outEdge;
		}

		/*
		Logger::get().debug() << "In";
		for (auto seg : prevEdge->seqSegments)
		{
			Logger::get().debug() << seg.seqId << " " << seg.start << " " << seg.end;
		}
		Logger::get().debug() << "Out";
		for (auto seg : nextEdge->seqSegments)
		{
			Logger::get().debug() << seg.seqId << " " << seg.start << " " << seg.end;
		}*/

		auto growingSeqs = prevEdge->seqSegments;
		std::vector<SequenceSegment> updatedSeqs;
		for (auto& seq : growingSeqs)
		{
			while (true)
			{
				bool updated = false;
				for (auto& otherSeq : loopEdge.seqSegments)
				{
					if (seq.seqId != otherSeq.seqId ||
						seq.end != otherSeq.start) continue;

					updated = true;
					seq.end = otherSeq.end;
					break;
				}

				bool complete = false;
				for (auto& otherSeq : nextEdge->seqSegments)
				{
					if (seq.seqId != otherSeq.seqId ||
						seq.end != otherSeq.start) continue;

					complete = true;
					break;
				}

				if (complete)
				{
					updatedSeqs.push_back(seq);
					break;
				}
				if (!updated) break;
			}
		}
		if (!updatedSeqs.empty())
		{
			prevEdge->seqSegments = updatedSeqs;
			return true;
		}
		else
		{
			Logger::get().debug() << "Can't unroll";
			return false;
		}
	};

	for (auto itEdge = _graph._graphEdges.begin(); 
		 itEdge != _graph._graphEdges.end();)
	{
		if (itEdge->nodeLeft == itEdge->nodeRight &&
			itEdge->nodeLeft->inEdges.size() == 2 &&
			itEdge->nodeLeft->outEdges.size() == 2 &&
			itEdge->nodeLeft->neighbors().size() == 3) 	//including itself
		{
			if (unrollEdge(*itEdge))
			{
				Logger::get().debug() << "Unroll " << itEdge->edgeId.signedId();
				vecRemove(itEdge->nodeLeft->inEdges, &(*itEdge));
				vecRemove(itEdge->nodeLeft->outEdges, &(*itEdge));
				itEdge = _graph._graphEdges.erase(itEdge);
				continue;
			}
		}
		++itEdge;
	}
}

void GraphProcessor::trimTips()
{

	int trimmed = 0;
	for (auto itEdge = _graph._graphEdges.begin(); 
		 itEdge != _graph._graphEdges.end();)
	{
		int prevDegree = itEdge->nodeLeft->inEdges.size();
		int nextDegree = itEdge->nodeRight->outEdges.size();
		if (itEdge->length() < _tipThreshold && 
			(prevDegree == 0 || nextDegree == 0))
		{
			GraphNode* toFix = (prevDegree != 0) ? itEdge->nodeLeft : 
								itEdge->nodeRight;
			//remove the edge
			vecRemove(itEdge->nodeRight->inEdges, &(*itEdge));
			vecRemove(itEdge->nodeLeft->outEdges, &(*itEdge));
			itEdge = _graph._graphEdges.erase(itEdge);

			++trimmed;

			//label all adjacent repetitive edges
			std::unordered_set<GraphNode*> visited;
			std::deque<GraphNode*> queue;
			queue.push_back(toFix);
			while (!queue.empty())
			{
				GraphNode* node = queue.back();
				queue.pop_back();
				if (visited.count(node)) continue;
				visited.insert(node);

				for (auto& edge : node->inEdges) 
				{
					if (edge->isRepetitive() &&
						edge->nodeLeft != edge->nodeRight)
					{
						_outdatedEdges.insert(edge->edgeId);
						queue.push_back(edge->nodeLeft);
					}
				}
				for (auto& edge : node->outEdges) 
				{
					if (edge->isRepetitive() &&
						edge->nodeLeft != edge->nodeRight)
					{
						_outdatedEdges.insert(edge->edgeId);
						queue.push_back(edge->nodeRight);
					}
				}
			}
			///
		}
		else
		{
			++itEdge;
		}
	}

	Logger::get().debug() << trimmed << " tips removed";
}

void GraphProcessor::condenceEdges()
{
	auto collapseEdges = [this](const GraphPath& edges)
	{
		std::list<SequenceSegment> growingSeqs(edges.front()->seqSegments.begin(),
											   edges.front()->seqSegments.end());
		assert(edges.size() > 1);
		for (size_t i = 1; i < edges.size(); ++i)
		{
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
				Logger::get().debug() << "Can't collape edge";
				return;
			}
		}

		///New edge
		_graph._graphEdges.emplace_back(edges.front()->nodeLeft,
								 edges.back()->nodeRight,
								 FastaRecord::Id(_graph._nextEdgeId++));
		std::copy(growingSeqs.begin(), growingSeqs.end(),
				  std::back_inserter(_graph._graphEdges.back().seqSegments));
		_graph._graphEdges.back().multiplicity = 
								std::max((int)growingSeqs.size(), 2);
		_graph._graphEdges.back().nodeLeft->outEdges
								.push_back(&_graph._graphEdges.back());
		_graph._graphEdges.back().nodeRight->inEdges
								.push_back(&_graph._graphEdges.back());
		Logger::get().debug() << "Added " 
						<< _graph._graphEdges.back().edgeId.signedId();
		///
		_outdatedEdges.insert(_graph._graphEdges.back().edgeId);
		
		std::unordered_set<GraphEdge*> toRemove(edges.begin(), edges.end());
		for (auto itEdge = _graph._graphEdges.begin(); 
			 itEdge != _graph._graphEdges.end(); )
		{
			if (toRemove.count(&(*itEdge)))
			{
				Logger::get().debug() << "Removed " << itEdge->edgeId.signedId();
				vecRemove(itEdge->nodeRight->inEdges, &(*itEdge));
				vecRemove(itEdge->nodeLeft->outEdges, &(*itEdge));
				_outdatedEdges.erase(itEdge->edgeId);

				itEdge = _graph._graphEdges.erase(itEdge);
			}
			else
			{
				++itEdge;
			}
		}
	};

	std::vector<GraphPath> toCollapse;
	std::unordered_set<FastaRecord::Id> usedDirections;
	for (auto& node : _graph._graphNodes)
	{
		if (!node.isBifurcation()) continue;

		for (auto& direction : node.outEdges)
		{
			if (usedDirections.count(direction->edgeId)) continue;
			usedDirections.insert(direction->edgeId);

			GraphNode* curNode = direction->nodeRight;
			GraphPath traversed;
			traversed.push_back(direction);
			while (!curNode->isBifurcation() &&
				   !curNode->outEdges.empty())
			{
				traversed.push_back(curNode->outEdges.front());
				curNode = curNode->outEdges.front()->nodeRight;
			}
			usedDirections.insert(traversed.back()->edgeId.rc());
			
			if (traversed.size() > 1)
			{
				toCollapse.emplace_back(std::move(traversed));
			}
		}
	}

	for (auto& edges : toCollapse)
	{
		GraphPath complEdges = _graph.complementPath(edges);
		collapseEdges(edges);
		collapseEdges(complEdges);
	}
}

void GraphProcessor::updateEdgesMultiplicity()
{
	bool anyChanges = true;
	while (anyChanges)
	{
		anyChanges = false;
		for (auto& node : _graph._graphNodes)
		{
			int numUnknown = 0;
			bool isOut = false;
			int inDegree = 0;
			int outDegree = 0;
			GraphEdge* unknownEdge = nullptr;

			for (auto& outEdge : node.outEdges)
			{
				if (outEdge->nodeLeft == outEdge->nodeRight) continue;

				if (!_outdatedEdges.count(outEdge->edgeId))
				{
					outDegree += outEdge->multiplicity;
				}
				else
				{
					++numUnknown;
					unknownEdge = outEdge;
					isOut = true;
				}
			}
			for (auto& inEdge : node.inEdges)
			{
				if (inEdge->nodeLeft == inEdge->nodeRight) continue;

				if (!_outdatedEdges.count(inEdge->edgeId))
				{
					inDegree += inEdge->multiplicity;
				}
				else
				{
					++numUnknown;
					unknownEdge = inEdge;
					isOut = false;
				}
			}
			int newMultiplicity = isOut ? inDegree - outDegree : 
										  outDegree - inDegree;

			if (numUnknown == 1)
			{
				Logger::get().debug() << "Updated: " 
					<< unknownEdge->edgeId.signedId() << " " 
					<< unknownEdge->multiplicity 
					<< " " << newMultiplicity;

				GraphEdge* complEdge = _graph.complementPath({unknownEdge}).front();
				unknownEdge->multiplicity = newMultiplicity;
				complEdge->multiplicity = newMultiplicity;
				_outdatedEdges.erase(unknownEdge->edgeId);
				_outdatedEdges.erase(complEdge->edgeId);
				anyChanges = true;
			}
		}
	}

	Logger::get().debug() << _outdatedEdges.size() << " edges were not updated!";
	for (auto edge : _outdatedEdges)
	{
		Logger::get().debug() << edge.signedId();
	}
}

