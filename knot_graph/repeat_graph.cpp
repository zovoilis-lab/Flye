#include "repeat_graph.h"
#include "overlap.h"
#include "vertex_index.h"
#include "config.h"
#include "disjoint_set.h"

struct Point
{
	Point(FastaRecord::Id curId = FastaRecord::ID_NONE, int32_t curPos = 0, 
		  FastaRecord::Id extId = FastaRecord::ID_NONE, int32_t extPos = 0):
		curId(curId), curPos(curPos), extId(extId), extPos(extPos) {}
	
	FastaRecord::Id curId;
	int32_t curPos;
	FastaRecord::Id extId;
	int32_t extPos;
};

struct SeqPoint
{
	SeqPoint(FastaRecord::Id seqId = FastaRecord::ID_NONE, int32_t pos = 0):
		seqId(seqId), pos(pos) {}
	
	FastaRecord::Id seqId;
	int32_t pos;
};

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
	this->fixTips();
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
			//unionSet(&overlapClusters[ovlp.curId].back(),
			//		 &overlapClusters[ovlp.extId].back());

			//reverse complement
			int32_t curLen = _asmSeqs.seqLen(ovlp.curId);
			int32_t extLen = _asmSeqs.seqLen(ovlp.extId);
			OverlapRange complOvlp = ovlp.complement(curLen, extLen);
			//

			overlapClusters[complOvlp.curId].emplace_back(complOvlp);
			overlapClusters[complOvlp.extId].emplace_back(complOvlp.reverse());
			//unionSet(&overlapClusters[complOvlp.curId].back(),
			//		 &overlapClusters[complOvlp.extId].back());
			//
		}
	}
	
	/*
	for (auto& ovlpHash : overlapClusters)
	{
		for (auto& ovlp1 : ovlpHash.second)
		{
			for (auto& ovlp2 : ovlpHash.second)
			{
				if (ovlp1.data.curIntersect(ovlp2.data) > -_maxSeparation)
				{
					unionSet(&ovlp1, &ovlp2);
				}
			}
		}
	}

	//getting cluster assignments
	std::unordered_map<DjsOverlap*, size_t> knotMappings;
	std::unordered_map<DjsOverlap*, size_t> reprToKnot;
	size_t nextKnotId = 0;

	for (auto& ovlpHash : overlapClusters)
	{
		for (auto& ovlp : ovlpHash.second)
		{
			if (!reprToKnot.count(findSet(&ovlp)))
			{
				reprToKnot[findSet(&ovlp)] = nextKnotId++;
			}
			knotMappings[&ovlp] = reprToKnot[findSet(&ovlp)];
		}
	}*/

	for (auto& seqClusterPair : overlapClusters)
	{
		//forming intersection-based clusters
		/*for (auto& ovlp : seqClusterPair.second)
		{
			ovlp.parent = &ovlp;
			ovlp.rank = 0;
		}*/

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

			//auto knotId = knotMappings[clustHash.second.front()];
			_repeatClusters[seqClusterPair.first]
					.emplace_back(seqClusterPair.first, 0,
								  clustStart, clustEnd);
		}
	}
}

void RepeatGraph::getGluepoints(const OverlapContainer& asmOverlaps)
{
	//cluster interval endpoints
	typedef SetNode<Point> Endpoint;
	std::list<Endpoint> endpoints;
	size_t pointId = 0;

	for (auto& seqOvlps : asmOverlaps.getOverlapIndex())
	{
		for (auto& ovlp : seqOvlps.second)
		{
			endpoints.emplace_back(Point(ovlp.curId, ovlp.curBegin, 
										 ovlp.extId, ovlp.extBegin));
			endpoints.emplace_back(Point(ovlp.curId, ovlp.curEnd, 
										 ovlp.extId, ovlp.extEnd));
		}

		_gluePoints[seqOvlps.first].emplace_back(pointId, seqOvlps.first, 0);
		_gluePoints[seqOvlps.first].emplace_back(pointId + 1, seqOvlps.first, 
											 _asmSeqs.seqLen(seqOvlps.first));
		pointId += 2;
	}

	for (auto& p1 : endpoints)
	{
		for (auto& p2 : endpoints)
		{
			if (p1.data.curId == p2.data.curId &&
				abs(p1.data.curPos - p2.data.curPos) <= _maxSeparation)
			{
				unionSet(&p1, &p2);
			}
		}
	}

	std::unordered_map<Endpoint*, std::vector<Endpoint*>> clusters;
	for (auto& endpoint : endpoints)
	{
		clusters[findSet(&endpoint)].push_back(&endpoint);
	}

	typedef SetNode<SeqPoint> GlueSet;
	std::list<GlueSet> glueSets;

	///used later
	auto propagateClusters = [&glueSets, this]
		(const std::vector<GluePoint>& clusterPoints)
	{
		//add all cluster-specific points into the shared set
		std::vector<GlueSet*> toMerge;
		for (auto& clustPt : clusterPoints)
		{
			bool used = false;
			for (auto& glueNode : glueSets)
			{
				if (glueNode.data.seqId == clustPt.seqId &&
					abs(glueNode.data.pos - clustPt.position) < _maxSeparation)
				{
					used = true;
					toMerge.push_back(&glueNode);
				}
			}
			if (!used)
			{
				glueSets.emplace_back(SeqPoint(clustPt.seqId, clustPt.position));
				toMerge.push_back(&glueSets.back());
			}
		}
		for (size_t i = 0; i < toMerge.size() - 1; ++i)
		{
			unionSet(toMerge[i], toMerge[i + 1]);
		}
	};
	///

	for (auto& clustEndpoints : clusters)
	{
		int64_t sum = 0;
		for (auto& ep : clustEndpoints.second)
		{
			sum += ep->data.curPos;
		}
		int32_t clusterXpos = sum / clustEndpoints.second.size();
		FastaRecord::Id clustSeq = clustEndpoints.second.front()->data.curId;

		std::vector<GluePoint> clusterPoints;
		clusterPoints.emplace_back(0, clustSeq, clusterXpos);

		std::list<Endpoint> extCoords;
		for (auto& ep : clustEndpoints.second)
		{
			extCoords.emplace_back(ep->data);
		}
		
		for (auto& ovlp : asmOverlaps.getOverlapIndex().at(clustSeq))
		{
			if (ovlp.curBegin <= clusterXpos && clusterXpos <= ovlp.curEnd)
			{
				int32_t projectedPos = clusterXpos - ovlp.curBegin + ovlp.extBegin;
				extCoords.emplace_back(Point(clustSeq, clusterXpos,
											 ovlp.extId, projectedPos));
			}
		}

		//cluster them
		for (auto& p1 : extCoords)
		{
			for (auto& p2 : extCoords)
			{
				if (p1.data.extId == p2.data.extId &&
					abs(p1.data.extPos - p2.data.extPos) <= _maxSeparation)
				{
					unionSet(&p1, &p2);
				}
			}
		}
		std::unordered_map<Endpoint*, std::vector<Endpoint*>> extClusters;
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
			clusterPoints.emplace_back(0, extSeq, clusterYpos);
		}

		//reverse complements
		std::vector<GluePoint> complClusterPoints;
		for (auto& gp : clusterPoints)
		{
			int32_t seqLen = _asmSeqs.seqLen(gp.seqId);
			complClusterPoints.emplace_back(0, gp.seqId.rc(), 
											seqLen - gp.position);
		}
		
		propagateClusters(clusterPoints);
		propagateClusters(complClusterPoints);
	}

	std::unordered_map<GlueSet*, size_t> setToId;
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
	}
}

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
}

bool RepeatGraph::isRepetitive(GluePoint gpLeft, GluePoint gpRight)
{
	for (auto& cluster : _repeatClusters[gpLeft.seqId])
	{
		int32_t overlap = std::min(gpRight.position, cluster.end) -
						  std::max(gpLeft.position, cluster.start);
		float rate = float(overlap) / (gpRight.position - gpLeft.position);

		if (rate > 0.5)
		{
			return true;
		}
	}
	return false;
}

void RepeatGraph::initializeEdges()
{
	std::unordered_map<size_t, GraphNode*> idToNode;
	typedef std::pair<GraphNode*, GraphNode*> NodePair;
	std::unordered_map<NodePair, GraphEdge*, pairhash> repeatEdges;

	auto addUnique = [&idToNode, this](GluePoint gpLeft, 
												GluePoint gpRight)
	{
		if (!idToNode.count(gpLeft.pointId))
		{
			_graphNodes.emplace_back();
			idToNode[gpLeft.pointId] = &_graphNodes.back();
		}
		if (!idToNode.count(gpRight.pointId))
		{
			_graphNodes.emplace_back();
			idToNode[gpRight.pointId] = &_graphNodes.back();
		}
		GraphNode* leftNode = idToNode[gpLeft.pointId];
		GraphNode* rightNode = idToNode[gpRight.pointId];

		_graphEdges.emplace_back(leftNode, rightNode, 
								 FastaRecord::Id(_nextEdgeId++));
		leftNode->outEdges.push_back(&_graphEdges.back());
		rightNode->inEdges.push_back(&_graphEdges.back());

		_graphEdges.back().addSequence(gpLeft.seqId, gpLeft.position, 
									   gpRight.position);
	};

	auto addRepeat = [&idToNode, this, &repeatEdges]
		(GluePoint gpLeft, GluePoint gpRight)
	{
		if (!idToNode.count(gpLeft.pointId))
		{
			_graphNodes.emplace_back();
			idToNode[gpLeft.pointId] = &_graphNodes.back();
		}
		if (!idToNode.count(gpRight.pointId))
		{
			_graphNodes.emplace_back();
			idToNode[gpRight.pointId] = &_graphNodes.back();
		}
		GraphNode* leftNode = idToNode[gpLeft.pointId];
		GraphNode* rightNode = idToNode[gpRight.pointId];

		if (!repeatEdges.count({leftNode, rightNode}))
		{
			_graphEdges.emplace_back(leftNode, rightNode, 
									 FastaRecord::Id(_nextEdgeId++));
			leftNode->outEdges.push_back(&_graphEdges.back());
			rightNode->inEdges.push_back(&_graphEdges.back());

			repeatEdges[std::make_pair(leftNode, rightNode)] = 
												&_graphEdges.back();
		}

		repeatEdges[std::make_pair(leftNode, rightNode)]
			->addSequence(gpLeft.seqId, gpLeft.position, gpRight.position);
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

			if (!repetitive)
			{
				addUnique(gpLeft, gpRight);
				addUnique(complLeft, complRight);
			}
			else
			{
				addRepeat(gpLeft, gpRight);
				addRepeat(complLeft, complRight);
			}
		}
	}
	
	/*
	for (auto& seqEdges : _graphEdges)
	{
		for (auto& edge : seqEdges.second)
		{
			Logger::get().debug() << edge.edgeId << "\t" 
								  << edge.gpLeft.seqId << "\t"
								  << edge.gpLeft.position << "\t" 
								  << edge.gpRight.position << "\t"
								  << edge.gpRight.position - edge.gpLeft.position
								  << "\t" << edge.repetitive;
		}
	}*/
}

RepeatGraph::GraphPath
RepeatGraph::chainReadAlignments(const SequenceContainer& edgeSeqs,
								 std::vector<EdgeAlignment> ovlps)
{
	const int32_t JUMP = 1500;

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

			if (JUMP > readDiff && readDiff > 0 &&
				JUMP > graphDiff && graphDiff > 0 &&
				abs(readDiff - graphDiff) < JUMP / 10 &&
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
		if (edgeAlignment.overlap.curBegin < Constants::maximumOverhang)
		{
			activeChains.push_back({&edgeAlignment});
		}
	}

	int32_t maxSpan = 0;
	Chain* maxChain = nullptr;
	for (auto& chain : activeChains)
	{
		//check right overhang
		int32_t overhang = _readSeqs.seqLen(chain.back()->overlap.curId) -
							chain.back()->overlap.curEnd;
		if (overhang < Constants::maximumOverhang) continue;

		int32_t readSpan = chain.back()->overlap.curEnd - 
						   chain.front()->overlap.curBegin;
		if (readSpan > maxSpan)
		{
			maxSpan = readSpan;
			maxChain = &chain;
		}
	}

	std::vector<GraphEdge*> result;
	if (maxChain && maxChain->size() > 1)
	{
		//chech number of non-repetitive edges
		int numUnique = 0;
		for (auto& edge : *maxChain) 
		{
			if (!edge->edge->isRepetitive()) ++numUnique;
		}
		if (numUnique < 2) return {};
		
		//check length consistency
		int32_t readSpan = maxChain->back()->overlap.curEnd - 
						   maxChain->front()->overlap.curBegin;
		int32_t graphSpan = maxChain->front()->overlap.curRange();
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

		for (auto& aln : *maxChain) result.push_back(aln->edge);
	}

	return result;
}

void RepeatGraph::fixTips()
{
	const int THRESHOLD = 1500;
	std::unordered_set<FastaRecord::Id> suspicious;
	std::unordered_map<FastaRecord::Id, GraphEdge*> idToEdge;
	for (auto& edge : _graphEdges) idToEdge[edge.edgeId] = &edge;

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
			Logger::get().debug() << "Chop!";

			for (auto& edge : toFix->inEdges) suspicious.insert(edge->edgeId);
			for (auto& edge : toFix->outEdges) suspicious.insert(edge->edgeId);
		}
		else
		{
			++itEdge;
		}
	}

	auto collapseEdges = [this, &suspicious](const GraphPath& edges)
	{
		Logger::get().debug() << "-----";
		for (auto& edge : edges) 
		{
			Logger::get().debug() << edge->edgeId << " " << edge->seqSegments.size();
			for (auto& seg : edge->seqSegments)
			{
				Logger::get().debug() << "\t" << seg.seqId << " " << seg.start << " " << seg.end;
			}
		}

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
			
		for (auto& seg : growingSeqs)
		{
			Logger::get().debug() << seg.seqId << " " << seg.start << " " << seg.end;
		}

		///New edge
		_graphEdges.emplace_back(edges.front()->nodeLeft,
								 edges.back()->nodeRight,
								 FastaRecord::Id(_nextEdgeId++));
		std::copy(growingSeqs.begin(), growingSeqs.end(),
				  std::back_inserter(_graphEdges.back().seqSegments));
		_graphEdges.back().multiplicity = growingSeqs.size();
		_graphEdges.back().nodeLeft->outEdges.push_back(&_graphEdges.back());
		_graphEdges.back().nodeRight->inEdges.push_back(&_graphEdges.back());
		Logger::get().debug() << "Added " << _graphEdges.back().edgeId;
		///
		suspicious.insert(_graphEdges.back().edgeId);
		
		std::unordered_set<GraphEdge*> toRemove(edges.begin(), edges.end());
		for (auto itEdge = _graphEdges.begin(); itEdge != _graphEdges.end(); )
		{
			if (toRemove.count(&(*itEdge)))
			{
				Logger::get().debug() << "Removed " << itEdge->edgeId;
				vecRemove(itEdge->nodeRight->inEdges, &(*itEdge));
				vecRemove(itEdge->nodeLeft->outEdges, &(*itEdge));
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

	bool anyChanges = true;
	while (anyChanges)
	{
		anyChanges = false;
		for (auto& edge : _graphEdges)
		{
			if (!suspicious.count(edge.edgeId)) continue;
			if (edge.nodeLeft == edge.nodeRight) continue;

			//check beginning
			bool allReliable = true;
			int inDegree = 0;
			int outDegree = 0;
			for (auto& otherEdge : edge.nodeLeft->outEdges)
			{
				if (otherEdge->edgeId == edge.edgeId) continue;
				if (otherEdge->nodeLeft == otherEdge->nodeRight) continue;

				if (!suspicious.count(otherEdge->edgeId))
				{
					outDegree += otherEdge->multiplicity;
				}
				else
				{
					allReliable = false;
				}
			}
			for (auto& otherEdge : edge.nodeLeft->inEdges)
			{
				if (otherEdge->nodeLeft == otherEdge->nodeRight) continue;

				if (!suspicious.count(otherEdge->edgeId))
				{
					inDegree += otherEdge->multiplicity;
				}
				else
				{
					allReliable = false;
				}
			}

			if (allReliable)
			{
				Logger::get().debug() << "Updated: "<< edge.edgeId << " " 
									  << edge.multiplicity 
									  << " " << inDegree - outDegree;
				edge.multiplicity = inDegree - outDegree;
				suspicious.erase(edge.edgeId);
				anyChanges = true;
			}
		}
	}
}

size_t RepeatGraph::separatePath(const GraphPath& graphPath, size_t startId)
{
	//first edge
	_graphNodes.emplace_back();
	vecRemove(graphPath.front()->nodeRight->inEdges, graphPath.front());
	graphPath.front()->nodeRight = &_graphNodes.back();
	GraphNode* prevNode = &_graphNodes.back();

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

		prevNode = nextNode;
	}

	//last edge
	vecRemove(graphPath.back()->nodeLeft->outEdges, graphPath.back());
	graphPath.back()->nodeLeft = prevNode;

	return edgesAdded;
}

RepeatGraph::GraphPath RepeatGraph::complementPath(const GraphPath& path)
{
	std::unordered_map<FastaRecord::Id, GraphEdge*> idToEdge;
	for (auto& edge : _graphEdges) idToEdge[edge.edgeId] = &edge;

	GraphPath complEdges;
	for (auto itEdge = path.rbegin(); itEdge != path.rend(); ++itEdge)
	{
		complEdges.push_back(idToEdge[(*itEdge)->edgeId.rc()]);
	}

	return complEdges;
}

void RepeatGraph::resolveConnections(const std::vector<GraphPath>& connections)
{
	const float MIN_CONFIDENCE = 0.7f;

	typedef std::pair<GraphEdge*, GraphEdge*> EdgePair;
	std::unordered_map<EdgePair, int, pairhash> edges;

	for (auto& conn : connections)
	{
		++edges[EdgePair(conn.front(), conn.back())];
	}
	std::vector<EdgePair> sortedConnections;
	for (auto& edgeCount : edges)
	{
		sortedConnections.push_back(edgeCount.first);
	}
	std::sort(sortedConnections.begin(), sortedConnections.end(),
			  [&edges](const EdgePair& e1, const EdgePair& e2)
			  {return edges[e1] > edges[e2];});

	std::unordered_set<FastaRecord::Id> usedInEdges;
	std::unordered_set<FastaRecord::Id> usedOutEdges;
	int totalLinks = 0;
	std::vector<GraphPath> uniquePaths;

	for (auto& edgePair : sortedConnections)
	{
		//
		int32_t coverage = 0;
		for (auto& conn : sortedConnections)
		{
			if (conn.first == edgePair.first &&
				!usedOutEdges.count(conn.second->edgeId))
			{
				coverage += edges[conn];
			}
			if (conn.second == edgePair.second &&
				!usedInEdges.count(conn.first->edgeId))
			{
				coverage += edges[conn];
			}
		}
		float confidence = 2.0f * edges[edgePair] / coverage;
		//

		if (!usedInEdges.count(edgePair.first->edgeId) &&
			!usedOutEdges.count(edgePair.second->edgeId) &&
			confidence > MIN_CONFIDENCE)
		{

			Logger::get().debug() << "\tConnection " 
				<< edgePair.first->seqSegments.front().seqId
				<< "\t" << edgePair.first->seqSegments.front().end << "\t"
				<< edgePair.second->seqSegments.front().seqId
				<< "\t" << edgePair.second->seqSegments.front().start
				<< "\t" << edges[edgePair]
				<< "\t" << confidence;

			usedInEdges.insert(edgePair.first->edgeId);
			usedInEdges.insert(edgePair.second->edgeId.rc());
			usedOutEdges.insert(edgePair.second->edgeId);
			usedOutEdges.insert(edgePair.first->edgeId.rc());
			totalLinks += 2;

			//TODO: get a path of median length instead, not just the first one?
			for (auto& path : connections)
			{
				if (path.front() == edgePair.first && 
					path.back() == edgePair.second)
				{
					uniquePaths.push_back(path);
					break;
				}
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
									500, 1500, 0);
	OverlapContainer readsContainer(readsOverlapper, _readSeqs);
	readsContainer.findAllOverlaps();

	//get connections
	std::vector<GraphPath> allConnections;
	for (auto& seqOverlaps : readsContainer.getOverlapIndex())
	{
		std::vector<EdgeAlignment> alignments;
		for (auto& ovlp : filterOvlp(seqOverlaps.second))
		{
			alignments.push_back({ovlp, idToSegment[ovlp.extId].first,
								  idToSegment[ovlp.extId].second});
		}

		auto readPath = this->chainReadAlignments(pathsContainer, alignments);
		_readsAlignments[seqOverlaps.first] = readPath;

		std::vector<GraphPath> readConnections;
		GraphPath currentPath;

		//Logger::get().debug() << _readSeqs.seqName(ovlps.front().overlap.curId);
		for (auto& edge : readPath)
		{
			/*Logger::get().debug() << edge->edge->gpLeft.seqId << "\t" 
								  << edge->edge->gpLeft.position << "\t"
								  << edge->edge->gpRight.position << "\t"
								  << edge->edge->repetitive << "\t"
								  << edge->edge->edgeId << "\t"
								  << edge->overlap.curBegin << "\t" 
								  << edge->overlap.curEnd;*/
			currentPath.push_back(edge);
			if (!edge->isRepetitive())
			{
				if (currentPath.size() > 1)
				{
					readConnections.push_back(currentPath);
					currentPath.clear();
					currentPath.push_back(edge);
				}
			}
		}

		for (auto& conn : readConnections)
		{
			allConnections.push_back(conn);
			allConnections.push_back(this->complementPath(conn));
		}
	}

	this->resolveConnections(allConnections);
}


void RepeatGraph::outputDot(const std::string& filename, bool collapseRepeats)
{
	const std::string COLORS[] = {"red", "darkgreen", "blue", "goldenrod", 
								  "cadetblue",
								  "darkorchid", "aquamarine1", "darkgoldenrod1",
								  "deepskyblue1", "darkolivegreen3"};
	size_t nextColor = 0;

	std::ofstream fout(filename);
	fout << "digraph {\n";
	
	std::unordered_map<GraphNode*, int> nodeIds;
	int nextId = 0;
	auto nodeToId = [&nodeIds, &nextId](GraphNode* node)
	{
		if (!nodeIds.count(node))
		{
			nodeIds[node] = nextId++;
		}
		return nodeIds[node];
	};

	for (auto& edge : _graphEdges)
	{
		if (!edge.isRepetitive())
		{
			fout << "\"" << nodeToId(edge.nodeLeft) << "\" -> \"" 
				 << nodeToId(edge.nodeRight)
				 << "\" [label = \"" << edge.edgeId << " "
				 << edge.length() << "\", color = \"black\"] ;\n";
		}
		else
		{
			std::string color = COLORS[nextColor];
			nextColor = (nextColor + 1) % 10;

			fout << "\"" << nodeToId(edge.nodeLeft) << "\" -> \"" 
				 << nodeToId(edge.nodeRight)
				 << "\" [label = \"" << edge.edgeId << " ("
				 << edge.multiplicity << ") "
				 << edge.length() << "\", color = \"" << color << "\" "
				 << " penwidth = 3] ;\n";
		}
	}

	fout << "}\n";
}
