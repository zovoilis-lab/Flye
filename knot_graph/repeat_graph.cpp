//(c) 2016 by Authors
//This file is a part of ABruijn program.
//Released under the BSD license (see LICENSE file)


#include "repeat_graph.h"
#include "overlap.h"
#include "vertex_index.h"
#include "config.h"
#include "disjoint_set.h"
#include "bipartie_mincost.h"

#include <deque>

namespace
{
	template<class T>
	void vecRemove(std::vector<T>& v, T val)
	{
		v.erase(std::remove(v.begin(), v.end(), val), v.end()); 
	}

	std::vector<OverlapRange> filterOvlp(const std::vector<OverlapRange>& ovlps)
	{
		std::vector<OverlapRange> filtered;
		for (auto& ovlp : ovlps)
		{
			bool found = false;
			for (auto& otherOvlp : filtered)
			{
				if (otherOvlp.curBegin == ovlp.curBegin &&
					otherOvlp.curEnd == ovlp.curEnd &&
					otherOvlp.extBegin == ovlp.extBegin &&
					otherOvlp.extEnd == ovlp.extEnd)
				{
					found = true;
					break;
				}
			}
			if (!found) filtered.push_back(ovlp);
		}
		return filtered;
	}

	struct pairhash 
	{
	public:
		template <typename T, typename U>
		std::size_t operator()(const std::pair<T, U> &x) const
		{
			return std::hash<T>()(x.first) ^ std::hash<U>()(x.second);
		}
	};

	struct Point2d
	{
		Point2d(FastaRecord::Id curId = FastaRecord::ID_NONE, int32_t curPos = 0, 
			  FastaRecord::Id extId = FastaRecord::ID_NONE, int32_t extPos = 0):
			curId(curId), curPos(curPos), extId(extId), extPos(extPos) {}
		
		FastaRecord::Id curId;
		int32_t curPos;
		FastaRecord::Id extId;
		int32_t extPos;
	};

	struct Point1d
	{
		Point1d(FastaRecord::Id seqId = FastaRecord::ID_NONE, int32_t pos = 0):
			seqId(seqId), pos(pos) {}
		
		FastaRecord::Id seqId;
		int32_t pos;
	};
}

void RepeatGraph::unrollLoops()
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

	for (auto itEdge = _graphEdges.begin(); itEdge != _graphEdges.end();)
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
				itEdge = _graphEdges.erase(itEdge);
				continue;
			}
		}
		++itEdge;
	}
}




void RepeatGraph::build()
{
	//getting overlaps
	VertexIndex asmIndex(_asmSeqs);
	asmIndex.countKmers(1);
	asmIndex.buildIndex(1, 500, 10);

	OverlapDetector asmOverlapper(_asmSeqs, asmIndex, 
								  Constants::maximumJump, 
								  Parameters::minimumOverlap,
								  0);
	OverlapContainer asmOverlaps(asmOverlapper, _asmSeqs);
	asmOverlaps.findAllOverlaps();

	this->getRepeatClusters(asmOverlaps);
	this->getGluepoints(asmOverlaps);
	this->initializeEdges();
}

void RepeatGraph::simplify()
{
	this->trimTips();
	this->unrollLoops();
	this->condenceEdges();
	this->updateEdgesMultiplicity();
}

void RepeatGraph::getRepeatClusters(const OverlapContainer& asmOverlaps)
{
	//forming overlap-based clusters
	typedef SetNode<OverlapRange> DjsOverlap;
	std::unordered_map<FastaRecord::Id, 
					   std::list<DjsOverlap>> overlapClusters;
	for (auto& ovlpHash : asmOverlaps.getOverlapIndex())
	{
		for (auto& ovlp : ovlpHash.second)
		{
			overlapClusters[ovlp.curId].emplace_back(ovlp);
			overlapClusters[ovlp.extId].emplace_back(ovlp.reverse());

			//reverse complement
			int32_t curLen = _asmSeqs.seqLen(ovlp.curId);
			int32_t extLen = _asmSeqs.seqLen(ovlp.extId);
			OverlapRange complOvlp = ovlp.complement(curLen, extLen);
			//

			overlapClusters[complOvlp.curId].emplace_back(complOvlp);
			overlapClusters[complOvlp.extId].emplace_back(complOvlp.reverse());
		}
	}

	for (auto& seqClusterPair : overlapClusters)
	{
		for (auto& ovlp1 : seqClusterPair.second)
		{
			for (auto& ovlp2 : seqClusterPair.second)
			{
				if (ovlp1.data.curIntersect(ovlp2.data) > -_maxSeparation)
				{
					unionSet(&ovlp1, &ovlp2);
				}
			}
		}
		std::unordered_map<DjsOverlap*, 
						   std::vector<DjsOverlap*>> intersectClusters;
		for (auto& ovlp : seqClusterPair.second)
		{
			intersectClusters[findSet(&ovlp)].push_back(&ovlp);
		}

		for (auto& clustHash : intersectClusters)
		{
			int32_t clustStart = std::numeric_limits<int32_t>::max();
			int32_t clustEnd = std::numeric_limits<int32_t>::min();
			for (auto& ovlp : clustHash.second)
			{
				clustStart = std::min(clustStart, ovlp->data.curBegin);
				clustEnd = std::max(clustEnd, ovlp->data.curEnd);
			}

			_repeatClusters[seqClusterPair.first]
					.emplace_back(seqClusterPair.first, 0,
								  clustStart, clustEnd);
		}
	}
}

void RepeatGraph::getGluepoints(const OverlapContainer& asmOverlaps)
{
	//cluster interval endpoints
	typedef SetNode<Point2d> SetPoint2d;
	std::list<SetPoint2d> endpoints;
	size_t pointId = 0;

	for (auto& seqOvlps : asmOverlaps.getOverlapIndex())
	{
		for (auto& ovlp : seqOvlps.second)
		{
			endpoints.emplace_back(Point2d(ovlp.curId, ovlp.curBegin, 
										 ovlp.extId, ovlp.extBegin));
			endpoints.emplace_back(Point2d(ovlp.curId, ovlp.curEnd, 
										 ovlp.extId, ovlp.extEnd));
		}
	}

	for (auto& p1 : endpoints)
	{
		for (auto& p2 : endpoints)
		{
			if (p1.data.curId == p2.data.curId &&
				abs(p1.data.curPos - p2.data.curPos) < _maxSeparation)
			{
				unionSet(&p1, &p2);
			}
		}
	}

	std::unordered_map<SetPoint2d*, std::vector<SetPoint2d*>> clusters;
	for (auto& endpoint : endpoints)
	{
		clusters[findSet(&endpoint)].push_back(&endpoint);
	}

	typedef SetNode<Point1d> SetPoint1d;
	std::list<SetPoint1d> glueSets;
	std::unordered_map<SetPoint1d*, SetPoint1d*> complements;

	for (auto& clustEndpoints : clusters)
	{
		FastaRecord::Id clustSeq = clustEndpoints.second.front()->data.curId;
		if (!clustSeq.strand()) continue;	//only for forward strands

		int64_t sum = 0;
		for (auto& ep : clustEndpoints.second)
		{
			sum += ep->data.curPos;
		}
		int32_t clusterXpos = sum / clustEndpoints.second.size();

		std::vector<Point1d> clusterPoints;
		clusterPoints.emplace_back(clustSeq, clusterXpos);

		std::list<SetPoint2d> extCoords;
		for (auto& ep : clustEndpoints.second)
		{
			extCoords.emplace_back(ep->data);
		}
		
		//projections on other overlaps
		for (auto& ovlp : asmOverlaps.getOverlapIndex().at(clustSeq))
		{
			if (ovlp.curBegin <= clusterXpos && clusterXpos <= ovlp.curEnd)
			{
				//TODO: projection with k-mers / alignment
				float lengthRatio = (float)ovlp.extRange() / ovlp.curRange();
				int32_t projectedPos = ovlp.extBegin + 
								float(clusterXpos - ovlp.curBegin) * lengthRatio;
				projectedPos = std::max(ovlp.extBegin, 
										std::min(projectedPos, ovlp.extEnd));

				extCoords.emplace_back(Point2d(clustSeq, clusterXpos,
									   		   ovlp.extId, projectedPos));
			}
		}

		//cluster them
		for (auto& p1 : extCoords)
		{
			for (auto& p2 : extCoords)
			{
				if (p1.data.extId == p2.data.extId &&
					abs(p1.data.extPos - p2.data.extPos) < _maxSeparation)
				{
					unionSet(&p1, &p2);
				}
			}
		}
		std::unordered_map<SetPoint2d*, std::vector<SetPoint2d*>> extClusters;
		for (auto& endpoint : extCoords)
		{
			extClusters[findSet(&endpoint)].push_back(&endpoint);
		}

		//now, get coordinates for each cluster
		for (auto& extClust : extClusters)
		{
			int64_t sum = 0;
			for (auto& ep : extClust.second)
			{
				sum += ep->data.extPos;
			}
			int32_t clusterYpos = sum / extClust.second.size();
			FastaRecord::Id extSeq = extClust.second.front()->data.extId;
			clusterPoints.emplace_back(extSeq, clusterYpos);
		}

		//merge with the previous clusters
		std::vector<SetPoint1d*> toMerge;
		for (auto& clustPt : clusterPoints)
		{
			int32_t seqLen = _asmSeqs.seqLen(clustPt.seqId);
			Point1d complPt(clustPt.seqId.rc(), seqLen - clustPt.pos - 1);

			bool used = false;
			for (auto& glueNode : glueSets)
			{
				if (glueNode.data.seqId == clustPt.seqId &&
					abs(glueNode.data.pos - clustPt.pos) < _maxSeparation)
				{
					used = true;
					toMerge.push_back(&glueNode);
				}
			}
			if (!used)
			{
				glueSets.emplace_back(clustPt);
				auto fwdPtr = &glueSets.back();
				glueSets.emplace_back(complPt);
				auto revPtr = &glueSets.back();

				complements[fwdPtr] = revPtr;
				complements[revPtr] = fwdPtr;
				toMerge.push_back(fwdPtr);
			}
		}
		for (size_t i = 0; i < toMerge.size() - 1; ++i)
		{
			unionSet(toMerge[i], toMerge[i + 1]);
			unionSet(complements[toMerge[i]], 
					 complements[toMerge[i + 1]]);
		}
	}

	std::unordered_map<SetPoint1d*, size_t> setToId;
	for (auto& setNode : glueSets)
	{
		if (!setToId.count(findSet(&setNode)))
		{
			setToId[findSet(&setNode)] = pointId++;
		}
		_gluePoints[setNode.data.seqId].emplace_back(setToId[findSet(&setNode)],
													 setNode.data.seqId,
													 setNode.data.pos);
	}

	for (auto& seqPoints : _gluePoints)
	{
		std::sort(seqPoints.second.begin(), seqPoints.second.end(),
				  [](const GluePoint& pt1, const GluePoint& pt2)
				  {return pt1.position < pt2.position;});

		//flanking points
		seqPoints.second.emplace(seqPoints.second.begin(), pointId++, 
								 seqPoints.first, 0);
		seqPoints.second.emplace_back(pointId++, seqPoints.first, 
									  _asmSeqs.seqLen(seqPoints.first) - 1);
	}
}

bool RepeatGraph::isRepetitive(GluePoint gpLeft, GluePoint gpRight)
{
	float maxRate = 0.0f;
	int32_t maxOverlap = 0;
	
	for (auto& cluster : _repeatClusters[gpLeft.seqId])
	{
		int32_t overlap = std::min(gpRight.position, cluster.end) -
						  std::max(gpLeft.position, cluster.start);
		float rate = float(overlap) / (gpRight.position - gpLeft.position);
		maxRate = std::max(maxRate, rate);
		maxOverlap = std::max(overlap, maxOverlap);

	}
	return maxRate > 0.5;
}

void RepeatGraph::initializeEdges()
{
	typedef std::pair<GraphNode*, GraphNode*> NodePair;
	std::unordered_map<NodePair, GraphEdge*, pairhash> repeatEdges;

	std::unordered_map<size_t, GraphNode*> nodeIndex;
	auto idToNode = [&nodeIndex, this](size_t nodeId)
	{
		if (!nodeIndex.count(nodeId))
		{
			_graphNodes.emplace_back();
			nodeIndex[nodeId] = &_graphNodes.back();
		}
		return nodeIndex[nodeId];
	};

	auto addUnique = [&idToNode, this](GluePoint gpLeft, 
									   GluePoint gpRight)
	{
		GraphNode* leftNode = idToNode(gpLeft.pointId);
		GraphNode* rightNode = idToNode(gpRight.pointId);

		_graphEdges.emplace_back(leftNode, rightNode, 
								 FastaRecord::Id(_nextEdgeId));
		leftNode->outEdges.push_back(&_graphEdges.back());
		rightNode->inEdges.push_back(&_graphEdges.back());
		++_nextEdgeId;

		_graphEdges.back().addSequence(gpLeft.seqId, gpLeft.position, 
									   gpRight.position);
	};

	auto addRepeat = [&idToNode, this, &repeatEdges]
		(GluePoint gpLeft, GluePoint gpRight, bool selfComplement)
	{
		GraphNode* leftNode = idToNode(gpLeft.pointId);
		GraphNode* rightNode = idToNode(gpRight.pointId);

		if (!repeatEdges.count({leftNode, rightNode}))
		{
			_graphEdges.emplace_back(leftNode, rightNode, 
									 FastaRecord::Id(_nextEdgeId));
			leftNode->outEdges.push_back(&_graphEdges.back());
			rightNode->inEdges.push_back(&_graphEdges.back());
			++_nextEdgeId;
			if (selfComplement) ++_nextEdgeId;

			repeatEdges[std::make_pair(leftNode, rightNode)] = 
												&_graphEdges.back();
		}

		GraphEdge* edge = repeatEdges[std::make_pair(leftNode, rightNode)];
		edge->selfComplement = selfComplement;
		edge->addSequence(gpLeft.seqId, gpLeft.position, gpRight.position);
	};

	for (auto& seqEdgesPair : _gluePoints)
	{
		if (!seqEdgesPair.first.strand()) continue;
		FastaRecord::Id complId = seqEdgesPair.first.rc();

		if (seqEdgesPair.second.size() != _gluePoints[complId].size())
		{
			throw std::runtime_error("Graph is not symmetric");
		}

		for (size_t i = 0; i < seqEdgesPair.second.size() - 1; ++i)
		{
			GluePoint gpLeft = seqEdgesPair.second[i];
			GluePoint gpRight = seqEdgesPair.second[i + 1];

			size_t complPos = seqEdgesPair.second.size() - i - 2;
			GluePoint complLeft = _gluePoints[complId][complPos];
			GluePoint complRight = _gluePoints[complId][complPos + 1];

			bool repetitive = this->isRepetitive(gpLeft, gpRight);
			if (this->isRepetitive(complLeft, complRight) != repetitive)
			{
				throw std::runtime_error("Complementary repeats error");
			}

			bool selfComplement = (gpLeft.pointId == gpRight.pointId &&
								   gpRight.pointId == complLeft.pointId &&
								   complLeft.pointId == complRight.pointId);
			if (!repetitive)
			{
				addUnique(gpLeft, gpRight);
				addUnique(complLeft, complRight);
			}
			else
			{
				addRepeat(gpLeft, gpRight, selfComplement);
				addRepeat(complLeft, complRight, selfComplement);
			}
			
			///
			FastaRecord::Id edgeId = FastaRecord::ID_NONE;
			if(!repetitive)
			{
				edgeId = FastaRecord::Id(_nextEdgeId - 2);
			}
			else
			{
				auto edge = repeatEdges[std::make_pair(idToNode(gpLeft.pointId), 
													idToNode(gpRight.pointId))];
				edgeId = edge->edgeId;
			}
			std::string unique = !repetitive ? "*" : " ";
			Logger::get().debug() << unique << "\t" << edgeId.signedId() << "\t" 
								  << gpLeft.seqId << "\t"
								  << gpLeft.position << "\t" 
								  << gpRight.position << "\t"
								  << gpRight.position - gpLeft.position;
			///
		}
	}
}

std::vector<RepeatGraph::EdgeAlignment>
	RepeatGraph::chainReadAlignments(const SequenceContainer& edgeSeqs,
								 	 std::vector<EdgeAlignment> ovlps)
{
	std::sort(ovlps.begin(), ovlps.end(),
			  [](const EdgeAlignment& e1, const EdgeAlignment& e2)
			  	{return e1.overlap.curBegin < e2.overlap.curBegin;});

	typedef std::vector<EdgeAlignment*> Chain;

	std::list<Chain> activeChains;
	for (auto& edgeAlignment : ovlps)
	{
		std::list<Chain> newChains;
		int32_t maxSpan = 0;
		Chain* maxChain = nullptr;
		for (auto& chain : activeChains)
		{
			OverlapRange& nextOvlp = edgeAlignment.overlap;
			OverlapRange& prevOvlp = chain.back()->overlap;

			int32_t readDiff = nextOvlp.curBegin - prevOvlp.curEnd;
			int32_t graphDiff = nextOvlp.extBegin +
							edgeSeqs.seqLen(prevOvlp.extId) - prevOvlp.extEnd;

			if (_readJump > readDiff && readDiff > 0 &&
				_readJump > graphDiff && graphDiff > 0 &&
				abs(readDiff - graphDiff) < _readJump / Constants::farJumpRate &&
				chain.back()->edge->nodeRight == 
				edgeAlignment.edge->nodeLeft)
			{
				int32_t readSpan = nextOvlp.curEnd -
								   chain.front()->overlap.curBegin;
				if (readSpan > maxSpan)
				{
					maxSpan = readSpan;
					maxChain = &chain;
				}
			}
		}
		
		if (maxChain)
		{
			newChains.push_back(*maxChain);
			maxChain->push_back(&edgeAlignment);
		}

		activeChains.splice(activeChains.end(), newChains);
		if (edgeAlignment.overlap.curBegin < _readOverhang)
		{
			activeChains.push_back({&edgeAlignment});
		}
	}

	int32_t maxSpan = 0;
	Chain* maxChain = nullptr;
	std::unordered_set<FastaRecord::Id> inEdges;
	std::unordered_set<FastaRecord::Id> outEdges;
	for (auto& chain : activeChains)
	{
		//check right overhang
		int32_t overhang = _readSeqs.seqLen(chain.back()->overlap.curId) -
							chain.back()->overlap.curEnd;
		if (overhang > _readOverhang) continue;

		int32_t readSpan = chain.back()->overlap.curEnd - 
						   chain.front()->overlap.curBegin;
		if (readSpan > maxSpan)
		{
			maxSpan = readSpan;
			maxChain = &chain;
		}

		inEdges.insert(chain.front()->edge->edgeId);
		outEdges.insert(chain.back()->edge->edgeId);
	}

	std::vector<EdgeAlignment> result;
	if (maxChain)
	{
		//chech number of non-repetitive edges
		int numUnique = 0;
		for (auto& edge : *maxChain) 
		{
			if (!edge->edge->isRepetitive()) ++numUnique;
		}
		if (numUnique < 2) return {};
		if (inEdges.size() != 1 || outEdges.size() != 1) return {};
		
		//check length consistency
		int32_t readSpan = maxChain->back()->overlap.curEnd - 
						   maxChain->front()->overlap.curBegin;
		int32_t graphSpan = maxChain->front()->overlap.extRange();
		for (size_t i = 1; i < maxChain->size(); ++i)
		{
			graphSpan += (*maxChain)[i]->overlap.extEnd +
						 edgeSeqs.seqLen((*maxChain)[i - 1]->overlap.extId) - 
						 (*maxChain)[i - 1]->overlap.extEnd;	
		}
		float lengthDiff = abs(readSpan - graphSpan);
		float meanLength = (readSpan + graphSpan) / 2.0f;
		if (lengthDiff > meanLength / Constants::overlapDivergenceRate)
		{
			return {};
		}

		for (auto& aln : *maxChain) result.push_back(*aln);
	}

	return result;
}

void RepeatGraph::trimTips()
{
	const int THRESHOLD = 1500;

	int trimmed = 0;
	for (auto itEdge = _graphEdges.begin(); itEdge != _graphEdges.end();)
	{
		int prevDegree = itEdge->nodeLeft->inEdges.size();
		int nextDegree = itEdge->nodeRight->outEdges.size();
		if (itEdge->length() < THRESHOLD && 
			(prevDegree == 0 || nextDegree == 0))
		{
			GraphNode* toFix = (prevDegree != 0) ? itEdge->nodeLeft : 
								itEdge->nodeRight;
			//remove the edge
			vecRemove(itEdge->nodeRight->inEdges, &(*itEdge));
			vecRemove(itEdge->nodeLeft->outEdges, &(*itEdge));
			itEdge = _graphEdges.erase(itEdge);

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

void RepeatGraph::condenceEdges()
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
		_graphEdges.emplace_back(edges.front()->nodeLeft,
								 edges.back()->nodeRight,
								 FastaRecord::Id(_nextEdgeId++));
		std::copy(growingSeqs.begin(), growingSeqs.end(),
				  std::back_inserter(_graphEdges.back().seqSegments));
		_graphEdges.back().multiplicity = std::max((int)growingSeqs.size(), 2);
		_graphEdges.back().nodeLeft->outEdges.push_back(&_graphEdges.back());
		_graphEdges.back().nodeRight->inEdges.push_back(&_graphEdges.back());
		Logger::get().debug() << "Added " << _graphEdges.back().edgeId.signedId();
		///
		_outdatedEdges.insert(_graphEdges.back().edgeId);
		
		std::unordered_set<GraphEdge*> toRemove(edges.begin(), edges.end());
		for (auto itEdge = _graphEdges.begin(); itEdge != _graphEdges.end(); )
		{
			if (toRemove.count(&(*itEdge)))
			{
				Logger::get().debug() << "Removed " << itEdge->edgeId.signedId();
				vecRemove(itEdge->nodeRight->inEdges, &(*itEdge));
				vecRemove(itEdge->nodeLeft->outEdges, &(*itEdge));
				_outdatedEdges.erase(itEdge->edgeId);

				itEdge = _graphEdges.erase(itEdge);
			}
			else
			{
				++itEdge;
			}
		}
	};

	std::vector<GraphPath> toCollapse;
	std::unordered_set<FastaRecord::Id> usedDirections;
	for (auto& node : _graphNodes)
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
		GraphPath complEdges = this->complementPath(edges);
		collapseEdges(edges);
		collapseEdges(complEdges);
	}
}

void RepeatGraph::updateEdgesMultiplicity()
{
	bool anyChanges = true;
	while (anyChanges)
	{
		anyChanges = false;
		for (auto& node : _graphNodes)
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

				GraphEdge* complEdge = this->complementPath({unknownEdge}).front();
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

size_t RepeatGraph::separatePath(const GraphPath& graphPath, size_t startId)
{
	/*
	Logger::get().debug() << "---";
	for (auto edge : graphPath)
	{
		Logger::get().debug() << edge->edgeId << " " << edge->multiplicity;
	}*/

	//first edge
	_graphNodes.emplace_back();
	vecRemove(graphPath.front()->nodeRight->inEdges, graphPath.front());
	graphPath.front()->nodeRight = &_graphNodes.back();
	GraphNode* prevNode = &_graphNodes.back();
	prevNode->inEdges.push_back(graphPath.front());

	//repetitive edges in the middle
	size_t edgesAdded = 0;
	for (size_t i = 1; i < graphPath.size() - 1; ++i)
	{
		--graphPath[i]->multiplicity;

		_graphNodes.emplace_back();
		GraphNode* nextNode = &_graphNodes.back();

		_graphEdges.emplace_back(prevNode, nextNode, 
								 FastaRecord::Id(startId));
		_graphEdges.back().seqSegments = graphPath[i]->seqSegments;
		_graphEdges.back().multiplicity = 1;
		startId += 2;
		edgesAdded += 1;

		prevNode->outEdges.push_back(&_graphEdges.back());
		nextNode->inEdges.push_back(&_graphEdges.back());
		prevNode = nextNode;
	}

	//last edge
	vecRemove(graphPath.back()->nodeLeft->outEdges, graphPath.back());
	graphPath.back()->nodeLeft = prevNode;
	prevNode->outEdges.push_back(graphPath.back());

	return edgesAdded;
}

RepeatGraph::GraphPath RepeatGraph::complementPath(const GraphPath& path)
{
	std::unordered_map<FastaRecord::Id, GraphEdge*> idToEdge;
	for (auto& edge : _graphEdges) 
	{
		idToEdge[edge.edgeId] = &edge;
		if (edge.selfComplement)
		{
			idToEdge[edge.edgeId.rc()] = &edge;
		}
	}

	GraphPath complEdges;
	for (auto itEdge = path.rbegin(); itEdge != path.rend(); ++itEdge)
	{
		complEdges.push_back(idToEdge.at((*itEdge)->edgeId.rc()));
	}

	return complEdges;
}

void RepeatGraph::resolveConnections(const std::vector<Connection>& connections)
{
	///////////
	std::unordered_map<GraphEdge*, std::unordered_map<GraphEdge*, int>> stats;
	for (auto& conn : connections)
	{
		++stats[conn.path.front()][conn.path.back()];
	}
	for (auto& leftEdge : stats)
	{
		Logger::get().debug() << "For " << leftEdge.first->edgeId << " "
			<< leftEdge.first->seqSegments.front().seqId << " "
			<< leftEdge.first->seqSegments.front().end;
		
		for (auto& rightEdge : leftEdge.second)
		{
			Logger::get().debug() << "\t" << rightEdge.first->edgeId << " "
				<< rightEdge.first->seqSegments.front().seqId << " "
				<< rightEdge.first->seqSegments.front().start << " " 
				<< rightEdge.second;
		}
		Logger::get().debug() << "";
	}
	///////////
	std::unordered_map<GraphEdge*, int> leftCoverage;
	std::unordered_map<GraphEdge*, int> rightCoverage;
	
	//create bipartie graph matrix
	std::unordered_map<GraphEdge*, size_t> leftEdgesId;
	std::unordered_map<size_t, GraphEdge*> leftIdToEdge;
	size_t nextLeftId = 0;
	std::unordered_map<GraphEdge*, size_t> rightEdgesId;
	std::unordered_map<size_t, GraphEdge*> rightIdToEdge;
	size_t nextRightId = 0;

	for (auto& conn : connections)
	{
		GraphEdge* leftEdge = conn.path.front();
		GraphEdge* rightEdge = conn.path.back();
		++leftCoverage[leftEdge];
		++rightCoverage[rightEdge];

		if (!leftEdgesId.count(leftEdge))
		{
			leftEdgesId[leftEdge] = nextLeftId;
			leftIdToEdge[nextLeftId++] = leftEdge;
		}
		if (!rightEdgesId.count(rightEdge))
		{
			rightEdgesId[rightEdge] = nextRightId;
			rightIdToEdge[nextRightId++] = rightEdge;
		}
	}

	size_t numNodes = std::max(leftEdgesId.size(), rightEdgesId.size());
	BipartieTable table;
	table.assign(numNodes, std::vector<double>(numNodes, 0));
	for (auto& conn : connections)
	{
		GraphEdge* leftEdge = conn.path.front();
		GraphEdge* rightEdge = conn.path.back();
		if (leftEdge->edgeId == rightEdge->edgeId ||
			leftEdge->edgeId == rightEdge->edgeId.rc()) continue;

		//solving min cost mathcing
		--table[leftEdgesId[leftEdge]][rightEdgesId[rightEdge]];
	}
	auto edges = bipartieMincost(table);
	typedef std::pair<size_t, size_t> MatchPair;
	std::vector<MatchPair> matchingPairs;
	for (size_t i = 0; i < edges.size(); ++i)
	{
		matchingPairs.emplace_back(i, edges[i]);
	}

	const float MIN_SUPPORT = 0.5f;
	std::unordered_set<FastaRecord::Id> usedEdges;
	std::vector<GraphPath> uniquePaths;
	int totalLinks = 0;
	for (auto match : matchingPairs)
	{
		GraphEdge* leftEdge = leftIdToEdge[match.first];
		GraphEdge* rightEdge = rightIdToEdge[match.second];

		int support = -table[match.first][match.second];
		float confidence = 2.0f * support / (leftCoverage[leftEdge] + 
											 rightCoverage[rightEdge]);
		if (!support) continue;
		if (usedEdges.count(leftEdge->edgeId)) continue;
		usedEdges.insert(rightEdge->edgeId.rc());

		Logger::get().debug() << "\tConnection " 
			<< leftEdge->seqSegments.front().seqId
			<< "\t" << leftEdge->seqSegments.front().end << "\t"
			<< rightEdge->seqSegments.front().seqId
			<< "\t" << rightEdge->seqSegments.front().start
			<< "\t" << support << "\t" << confidence;

		if (support < MIN_SUPPORT) continue;

		totalLinks += 2;
		for (auto& conn : connections)
		{
			if (conn.path.front() == leftEdge && 
				conn.path.back() == rightEdge)
			{
				uniquePaths.push_back(conn.path);
				break;
			}
		}
	}

	for (auto& path : uniquePaths)
	{
		GraphPath complPath = this->complementPath(path);
		size_t addedFwd = this->separatePath(path, _nextEdgeId);
		size_t addedRev = this->separatePath(complPath, _nextEdgeId + 1);
		assert(addedFwd == addedRev);
		_nextEdgeId += addedFwd + addedRev;
	}

	Logger::get().debug() << "Edges: " << totalLinks << " links: "
						  << connections.size();
}


void RepeatGraph::resolveRepeats()
{
	//create database
	std::unordered_map<FastaRecord::Id, 
					   std::pair<GraphEdge*, SequenceSegment*>> idToSegment;
	SequenceContainer pathsContainer;

	size_t nextSeqId = 0;
	for (auto& edge : _graphEdges)
	{
		for (auto& segment : edge.seqSegments)
		{
			size_t len = segment.end - segment.start;
			std::string sequence = _asmSeqs.getSeq(segment.seqId)
												.substr(segment.start, len);
			pathsContainer.addSequence(FastaRecord(sequence, "",
												   FastaRecord::Id(nextSeqId)));
			idToSegment[FastaRecord::Id(nextSeqId)] = {&edge, &segment};
			++nextSeqId;
		}
	}

	//index it and align reads
	VertexIndex pathsIndex(pathsContainer);
	pathsIndex.countKmers(1);
	pathsIndex.buildIndex(1, 5000, 1);
	OverlapDetector readsOverlapper(pathsContainer, pathsIndex, 
									_readJump, _maxSeparation, 0);
	OverlapContainer readsContainer(readsOverlapper, _readSeqs);
	readsContainer.findAllOverlaps();

	//get connections
	std::vector<Connection> readConnections;
	for (auto& seqOverlaps : readsContainer.getOverlapIndex())
	{
		std::vector<EdgeAlignment> alignments;
		for (auto& ovlp : filterOvlp(seqOverlaps.second))
		{
			alignments.push_back({ovlp, idToSegment[ovlp.extId].first,
								  idToSegment[ovlp.extId].second});
		}

		auto readPath = this->chainReadAlignments(pathsContainer, alignments);

		/*
		if (!readPath.empty())
		{
			Logger::get().debug() << _readSeqs.seqName(seqOverlaps.second.front().curId);
		}*/
		//std::vector<GraphPath> readConnections;
		GraphPath currentPath;
		int32_t prevFlank = 0;
		for (auto& aln : readPath)
		{
			if (currentPath.empty()) 
			{
				if (aln.edge->isRepetitive()) continue;
				prevFlank = aln.overlap.curRange();
			}
			/*
			if (!aln.edge->isRepetitive())
			{
				Logger::get().debug() << aln.edge->edgeId << "\t" 
									  << aln.edge->seqSegments.front().seqId << "\t"
									  << aln.edge->seqSegments.front().start << "\t"
									  << aln.edge->seqSegments.front().end << "\t"
									  << aln.overlap.curRange();
			}*/

			currentPath.push_back(aln.edge);
			if (!aln.edge->isRepetitive() && currentPath.size() > 1)
			{
				int32_t curFlank = aln.overlap.curRange();
				readConnections.push_back({currentPath, prevFlank, curFlank});
				readConnections.push_back({this->complementPath(currentPath),
										   curFlank, prevFlank});

				currentPath.clear();
				currentPath.push_back(aln.edge);
				prevFlank = curFlank;
			}
		}
	}

	this->resolveConnections(readConnections);

	//this->unrollLoops();
	//this->condenceEdges();
	//this->updateEdgesMultiplicity();
}


void RepeatGraph::outputDot(const std::string& filename, 
							bool collapseUnbranching)
{
	std::ofstream fout(filename);
	fout << "digraph {\n";
	
	///re-enumerating helper functions
	std::unordered_map<GraphNode*, int> nodeIds;
	int nextNodeId = 0;
	auto nodeToId = [&nodeIds, &nextNodeId](GraphNode* node)
	{
		if (!nodeIds.count(node))
		{
			nodeIds[node] = nextNodeId++;
		}
		return nodeIds[node];
	};

	std::unordered_map<FastaRecord::Id, size_t> edgeIds;
	size_t nextEdgeId = 0;
	auto pathToId = [&edgeIds, &nextEdgeId](GraphPath path)
	{
		if (!edgeIds.count(path.front()->edgeId))
		{
			edgeIds[path.front()->edgeId] = nextEdgeId;
			edgeIds[path.back()->edgeId.rc()] = nextEdgeId + 1;
			nextEdgeId += 2;
		}
		return FastaRecord::Id(edgeIds[path.front()->edgeId]);
	};
	
	const std::string COLORS[] = {"red", "darkgreen", "blue", "goldenrod", 
								  "cadetblue", "darkorchid", "aquamarine1", 
								  "darkgoldenrod1", "deepskyblue1", 
								  "darkolivegreen3"};
	std::unordered_map<FastaRecord::Id, size_t> colorIds;
	size_t nextColorId = 0;
	auto idToColor = [&colorIds, &nextColorId, &COLORS](FastaRecord::Id id)
	{
		if (!id.strand()) id = id.rc();
		if (!colorIds.count(id))
		{
			colorIds[id] = nextColorId;
			nextColorId = (nextColorId + 1) % 10;
		}
		return COLORS[colorIds[id]];
	};
	/////////////

	for (auto& node : _graphNodes)
	{
		if (!node.isBifurcation()) continue;

		for (auto& direction : node.outEdges)
		{

			GraphNode* curNode = direction->nodeRight;
			GraphPath traversed;
			traversed.push_back(direction);

			while (!curNode->isBifurcation() &&
				   !curNode->outEdges.empty())
			{
				traversed.push_back(curNode->outEdges.front());
				curNode = curNode->outEdges.front()->nodeRight;
			}
			
			bool resolvedRepeat = true;
			for (auto& edge : traversed)
			{
				if (edge->multiplicity != 0) resolvedRepeat = false;
			}
			if (resolvedRepeat) continue;

			///
			if (collapseUnbranching)
			{
				int32_t edgeLength = 0;
				for (auto& edge : traversed) edgeLength += edge->length();

				if (traversed.front()->isRepetitive())
				{
					FastaRecord::Id edgeId = pathToId(traversed);
					std::string color = idToColor(edgeId);
					//std::string color = _outdatedEdges.count(edgeId) ? "red" : "green";

					fout << "\"" << nodeToId(traversed.front()->nodeLeft) 
						 << "\" -> \"" << nodeToId(traversed.back()->nodeRight)
						 << "\" [label = \"" << edgeId.signedId() << 
						 " " << edgeLength << " (" 
						 << traversed.front()->multiplicity << ")\", color = \"" 
						 << color << "\" " << " penwidth = 3] ;\n";
				}
				else
				{
					fout << "\"" << nodeToId(traversed.front()->nodeLeft) 
						 << "\" -> \"" << nodeToId(traversed.back()->nodeRight)
						 << "\" [label = \"" << pathToId(traversed).signedId() 
						 << " " << edgeLength << "\", color = \"black\"] ;\n";
				}
			}
			else
			{
				for (auto& edge : traversed)
				{
					if (edge->isRepetitive())
					{
						FastaRecord::Id edgeId = edge->edgeId;
						std::string color = idToColor(edgeId);
						//std::string color = _outdatedEdges.count(edgeId) ? "red" : "green";

						fout << "\"" << nodeToId(edge->nodeLeft) 
							 << "\" -> \"" << nodeToId(edge->nodeRight)
							 << "\" [label = \"" << edgeId.signedId() 
							 << " " << edge->length() << " ("
							 << edge->multiplicity << ")" << "\", color = \"" 
							 << color << "\" " << " penwidth = 3] ;\n";
					}
					else
					{
						fout << "\"" << nodeToId(edge->nodeLeft) 
							 << "\" -> \"" << nodeToId(edge->nodeRight)
							 << "\" [label = \"" << edge->edgeId.signedId() << " "
							 << edge->length() << "\", color = \"black\"] ;\n";
					}
				}
			}
		}
	}

	fout << "}\n";
}
