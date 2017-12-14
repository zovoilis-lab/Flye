//(c) 2017 by Authors
//This file is a part of ABruijn program.
//Released under the BSD license (see LICENSE file)

#include "structure_resolver.h"

/*void Scaffolder::unrollLoops()
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
		int mult = std::max(loop->multiplicity, 1);
		_graph.removeEdge(loop);
		//if (mult == 0) return;

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

		int mult = std::max(edge->multiplicity, 1);
		for (int i = 0; i < mult; ++i)
		{
			fwdIds.emplace_back(_graph.newEdgeId());
			revIds.emplace_back(fwdIds.back().rc());
		}
		std::reverse(revIds.begin(), revIds.end());

		Logger::get().debug() << "Unrolling " << edge->edgeId.signedId()
			<< " of multiplicity " << mult;

		GraphEdge* complLoop = _graph.complementEdge(edge);
		unroll(leftEdge, edge, rightEdge, segment, fwdIds);
		unroll(_graph.complementEdge(rightEdge), complLoop,
			   _graph.complementEdge(leftEdge), segment.complement(), revIds);
	}
}*/

void StructureResolver::scaffold()
{
	auto connections = this->getScaffoldPairs();
	for (auto& conn : connections)
	{

		this->untangle(conn);
	}
}

void StructureResolver::untangle(ScaffoldInfo scfInfo)
{
	auto bestExtendingRead = [this] (GraphEdge* edge)
	{
		int32_t maxExtension = 0;
		GraphAlignment bestAlignment;
		for (auto& path : _aligner.getAlignments())
		{
			for (size_t i = 0; i < path.size(); ++i)
			{
				if (path[i].edge == edge && i < path.size() - 1)
				{
					int32_t alnLen = path.back().overlap.curEnd - 
									 path[i + 1].overlap.curBegin;
					if (alnLen > maxExtension)
					{
						maxExtension = alnLen;
						bestAlignment.clear();
						std::copy(path.begin() + i + 1, path.end(),
								  std::back_inserter(bestAlignment));
					}
					break;
				}
			}
		}
		return bestAlignment;
	};

	auto leftAln = bestExtendingRead(scfInfo.leftUnique);
	auto rightAln = bestExtendingRead(_graph.complementEdge(scfInfo.rightUnique));
	if (leftAln.empty() || rightAln.empty()) return;

	auto separate = [this] (GraphEdge* leftEdge, GraphEdge* rightEdge,
							const std::vector<SequenceSegment>& newSegments,
							const std::vector<FastaRecord::Id>& newIds)
	{
		GraphNode* leftNode = _graph.addNode();
		vecRemove(leftEdge->nodeRight->inEdges, leftEdge);
		leftEdge->nodeRight = leftNode;
		leftNode->inEdges.push_back(leftEdge);

		GraphNode* rightNode = leftNode;
		for (size_t i = 0; i < newSegments.size(); ++i)
		{
			rightNode = _graph.addNode();
			GraphEdge* newEdge = _graph.addEdge(GraphEdge(leftNode, rightNode,
														  newIds[i]));
			newEdge->seqSegments.push_back(newSegments[i]);
			leftNode = rightNode;
			Logger::get().debug() << "\tAdded: " << newIds[i].signedId();
		}

		vecRemove(rightEdge->nodeLeft->outEdges, rightEdge);
		rightEdge->nodeLeft = rightNode;
		rightNode->outEdges.push_back(rightEdge);
	};

	/*int32_t remainedGap = 0;
	for (auto& edge : scfInfo.repetitiveEdges)
	{
		for (auto& seg : edge->seqSegments) remainedGap += seg.length();
	}*/
	//Logger::get().debug() << "\tGraph gap: " << remainedGap;

	Logger::get().debug() << "Scaffold: " << scfInfo.leftUnique->edgeId.signedId()
		<< " -> " << scfInfo.rightUnique->edgeId.signedId();
	for (auto& edge : scfInfo.repetitiveEdges)
	{
		Logger::get().debug() << "\tRepeat " << edge->edgeId.signedId();
	}

	int32_t spannedGap = 0;
	spannedGap += leftAln.back().overlap.curEnd - 
					leftAln.front().overlap.curBegin;
	spannedGap += rightAln.back().overlap.curEnd - 
					rightAln.front().overlap.curBegin;
	Logger::get().debug() << "\tSpanned gap: " << spannedGap;

	std::vector<SequenceSegment> newSegments;
	std::vector<FastaRecord::Id> newIds;
	//left side
	SequenceSegment leftSeg(leftAln.front().overlap.curId, 
							leftAln.front().overlap.curLen,
							leftAln.front().overlap.curBegin,
							leftAln.back().overlap.curEnd);
	leftSeg.segType = SequenceSegment::Read;
	newSegments.push_back(leftSeg);
	newIds.push_back(_graph.newEdgeId());

	//remaining gap
	SequenceSegment gap;
	//gap.segType = SequenceSegment::Gap;
	newSegments.push_back(gap);
	newIds.push_back(_graph.newEdgeId());
	
	//right side
	SequenceSegment rightSeg(leftAln.front().overlap.curId, 
							 leftAln.front().overlap.curLen,
							 leftAln.front().overlap.curBegin,
							 leftAln.back().overlap.curEnd);
	rightSeg.segType = SequenceSegment::Read;
	newSegments.push_back(leftSeg.complement());
	newIds.push_back(_graph.newEdgeId());

	//keeping graph symmetric
	std::vector<SequenceSegment> complSegments;
	for (auto& seg : newSegments) complSegments.push_back(seg.complement());
	std::reverse(complSegments.begin(), complSegments.end());
	std::vector<FastaRecord::Id> complIds;
	for (auto& id : newIds) complIds.push_back(id.rc());
	std::reverse(complIds.begin(), complIds.end());

	separate(scfInfo.leftUnique, scfInfo.rightUnique, newSegments, newIds);
	separate(_graph.complementEdge(scfInfo.rightUnique), 
			 _graph.complementEdge(scfInfo.leftUnique), 
			 complSegments, complIds);
}

std::vector<StructureResolver::ScaffoldInfo> StructureResolver::getScaffoldPairs()
{
	std::vector<ScaffoldInfo> connectionPairs;

	for (auto& edge : _graph.iterEdges())
	{
		if (edge->repetitive) continue;

		std::vector<GraphEdge*> dfsStack;
		std::unordered_set<GraphEdge*> visited;
		std::unordered_set<GraphEdge*> reachableUnique;
		std::unordered_set<GraphEdge*> traversedRepeats;

		dfsStack.push_back(edge);
		while(!dfsStack.empty())
		{
			auto curEdge = dfsStack.back(); 
			dfsStack.pop_back();
			if (visited.count(curEdge)) continue;

			visited.insert(curEdge);
			visited.insert(_graph.complementEdge(curEdge));
			for (auto& adjEdge: curEdge->nodeRight->outEdges)
			{
				if (adjEdge->isRepetitive() && !visited.count(adjEdge))
				{
					dfsStack.push_back(adjEdge);
					traversedRepeats.insert(adjEdge);
				}
				else if (!adjEdge->isRepetitive() && adjEdge != edge)
				{
					reachableUnique.insert(adjEdge);
				}
			}
		}

		if (reachableUnique.size() == 1)
		{
			GraphEdge* outEdge = *reachableUnique.begin();
			if ((edge->nodeRight->isBifurcation() || 
					outEdge->nodeLeft->isBifurcation()) &&
				edge->edgeId != outEdge->edgeId.rc() &&
				abs(edge->edgeId.signedId()) < abs(outEdge->edgeId.signedId()))
			{
				connectionPairs.push_back({edge, outEdge, 
										   traversedRepeats});
			}
		}
	}

	return connectionPairs;
}
