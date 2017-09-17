//(c) 2016 by Authors
//This file is a part of ABruijn program.
//Released under the BSD license (see LICENSE file)

#include <deque>
#include <iomanip>

#include "../sequence/overlap.h"
#include "../sequence/vertex_index.h"
#include "../common/config.h"
#include "../common/disjoint_set.h"
#include "repeat_graph.h"


namespace
{
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
	
	template<typename T>
	T median(std::vector<T>& vec)
	{
		std::sort(vec.begin(), vec.end());
		//NOTE: there's a bug in libstdc++ nth_element, 
		//that sometimes leads to a segfault
		//std::nth_element(vec.begin(), vec.begin() + vec.size() / 2, 
		//				 vec.end());
		return vec[vec.size() / 2];
	}
}

bool GraphEdge::isTip() const
{
	return nodeLeft->inEdges.empty() || nodeRight->outEdges.empty();
}

void RepeatGraph::build()
{
	//getting overlaps
	VertexIndex asmIndex(_asmSeqs);
	asmIndex.countKmers(1);
	asmIndex.buildIndex(1, Constants::repeatGraphMaxKmer, 
						Constants::repeatGraphKmerSample);

	OverlapDetector asmOverlapper(_asmSeqs, asmIndex, 
								  Constants::maximumJump, 
								  Parameters::get().minimumOverlap,
								  /*no overhang*/ 0, /*keep alignment*/ true);
	OverlapContainer asmOverlaps(asmOverlapper, _asmSeqs, false);
	asmOverlaps.findAllOverlaps();

	this->getGluepoints(asmOverlaps);
	this->collapseTandems();
	this->initializeEdges(asmOverlaps);
}


void RepeatGraph::getGluepoints(const OverlapContainer& asmOverlaps)
{
	//Warning - The most complex function ever, a bit dirty too :(
	//Note, that there are two kinds of clustering here:
	//based on X and Y coordinates (position + seqId) of the points,
	//that might often define different subsets
	
	Logger::get().debug() << "Computing gluepoints";
	typedef SetNode<Point2d> SetPoint2d;
	std::unordered_map<FastaRecord::Id, 
					   std::vector<SetPoint2d*>> endpoints;

	//first, extract endpoints from all overlaps.
	//each point has X and Y coordinates (curSeq and extSeq)
	for (auto& seqOvlps : asmOverlaps.getOverlapIndex())
	{
		for (auto& ovlp : seqOvlps.second)
		{
			endpoints[ovlp.curId]
				.push_back(new SetPoint2d(Point2d(ovlp.curId, ovlp.curBegin,
										  ovlp.extId, ovlp.extBegin)));
			endpoints[ovlp.curId]
				.push_back(new SetPoint2d(Point2d(ovlp.curId, ovlp.curEnd,
										  ovlp.extId, ovlp.extEnd)));
		}
	}

	//for each contig, cluster gluepoints that are close to each other
	//only cosider X coordinates for now
	for (auto& seqPoints : endpoints)
	{
		std::sort(seqPoints.second.begin(), seqPoints.second.end(),
				  [](const SetPoint2d* p1, const SetPoint2d* p2)
				  {return p1->data.curPos < p2->data.curPos;});

		for (size_t i = 0; i < seqPoints.second.size() - 1; ++i)
		{
			auto* p1 = seqPoints.second[i];
			auto* p2 = seqPoints.second[i + 1];
			if (abs(p1->data.curPos - p2->data.curPos) < _maxSeparation)
			{
				unionSet(p1, p2);
			}

		}
	}
	std::unordered_map<SetPoint2d*, std::vector<SetPoint2d*>> clusters;
	for (auto seqPoints : endpoints)
	{
		for (auto& endpoint : seqPoints.second)
		{
			clusters[findSet(endpoint)].push_back(endpoint);
		}
	}

	typedef SetNode<Point1d> SetPoint1d;
	std::unordered_map<FastaRecord::Id, std::list<SetPoint1d>> tempGluepoints;
	std::unordered_map<SetPoint1d*, SetPoint1d*> complements;

	//we will now split each cluster based on it's Y coordinates
	//and project these subgroups to the corresponding sequences
	for (auto& clustEndpoints : clusters)
	{
		//first, simply add projections for each point from the cluster
		FastaRecord::Id clustSeq = clustEndpoints.second.front()->data.curId;
		if (!clustSeq.strand()) continue;	//only for forward strands

		std::vector<int32_t> positions;
		for (auto& ep : clustEndpoints.second) 
		{
			positions.push_back(ep->data.curPos);
		}
		int32_t clusterXpos = median(positions);

		std::vector<Point1d> clusterPoints;
		clusterPoints.emplace_back(clustSeq, clusterXpos);

		std::list<SetPoint2d> extCoords;
		for (auto& ep : clustEndpoints.second)
		{
			extCoords.emplace_back(ep->data);
		}
		
		//Important part: we need also add extra projections
		//for gluepoints that lie within other existing overlaps
		//(handles situations with 'repeat hierarchy', when some
		//repeats are parts of the other bigger repeats)
		auto startCmp = [] (const OverlapRange& ovlp, int32_t pos)
						{return ovlp.curBegin < pos;};
		auto& allOvlp = asmOverlaps.getOverlapIndex().at(clustSeq);
		auto leftPos = std::lower_bound(allOvlp.begin(), allOvlp.end(), 
										clusterXpos - 50000, startCmp);
		auto rightPos = std::lower_bound(allOvlp.begin(), allOvlp.end(), 
										 clusterXpos, startCmp);

		for (auto& ovlp = leftPos; ovlp != rightPos; ++ovlp)
		{
			if (clusterXpos - ovlp->curBegin >= 0 && 
				ovlp->curEnd - clusterXpos >= 0)
			{
				int32_t projectedPos = ovlp->project(clusterXpos);
				extCoords.emplace_back(Point2d(clustSeq, clusterXpos,
									   		   ovlp->extId, projectedPos));
			}
		}

		//Finally, cluster the projected points based on Y coordinates
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
			std::vector<int32_t> positions;
			for (auto& ep : extClust.second) 
			{
				positions.push_back(ep->data.extPos);
			}
			int32_t clusterYpos = median(positions);


			FastaRecord::Id extSeq = extClust.second.front()->data.extId;
			clusterPoints.emplace_back(extSeq, clusterYpos);
		}

		//We should now consider how newly generaetd clusters
		//are integrated with the existing ones, we might need to
		//merge some of them together
		std::vector<SetPoint1d*> toMerge;
		for (auto& clustPt : clusterPoints)
		{
			int32_t seqLen = _asmSeqs.seqLen(clustPt.seqId);
			Point1d complPt(clustPt.seqId.rc(), seqLen - clustPt.pos - 1);

			auto& seqGluepoints = tempGluepoints[clustPt.seqId];
			auto& complGluepoints = tempGluepoints[clustPt.seqId.rc()];
			for (auto& glueNode : seqGluepoints)
			{
				if (abs(glueNode.data.pos - clustPt.pos) < _maxSeparation)
				{
					toMerge.push_back(&glueNode);
				}
			}

			seqGluepoints.emplace_back(clustPt);
			auto fwdPtr = &seqGluepoints.back();
			complGluepoints.emplace_back(complPt);
			auto revPtr = &complGluepoints.back();

			complements[fwdPtr] = revPtr;
			complements[revPtr] = fwdPtr;
			toMerge.push_back(fwdPtr);
		}
		for (size_t i = 0; i < toMerge.size() - 1; ++i)
		{
			unionSet(toMerge[i], toMerge[i + 1]);
			unionSet(complements[toMerge[i]], 
					 complements[toMerge[i + 1]]);
		}
	}

	//Generating final gluepoints, we might need additionally
	//group clusters that are too close
	size_t pointId = 0;
	std::unordered_map<SetPoint1d*, size_t> setToId;
	auto addConsensusPoint = [&setToId, this, &pointId]
		(const std::vector<SetPoint1d*>& group)
	{
		SetPoint1d* reprPoint = group.front();
		if (!setToId.count(findSet(reprPoint)))
		{
			setToId[findSet(reprPoint)] = pointId++;
		}

		std::vector<int32_t> positions;
		for (auto& ep : group) 
		{
			positions.push_back(ep->data.pos);
		}
		int32_t clusterXpos = median(positions);

		bool tandemRepeat = group.back()->data.pos - 
							group.front()->data.pos > _maxSeparation;
		
		if (tandemRepeat)
		{
			_gluePoints[reprPoint->data.seqId]
				.emplace_back(setToId[findSet(reprPoint)],
							  reprPoint->data.seqId, group.front()->data.pos);

			_gluePoints[reprPoint->data.seqId]
				.emplace_back(setToId[findSet(reprPoint)],
							  reprPoint->data.seqId, group.back()->data.pos);
		}
		else
		{
			_gluePoints[reprPoint->data.seqId]
				.emplace_back(setToId[findSet(reprPoint)],
							  reprPoint->data.seqId, clusterXpos);
		}

	};
	int numGluepoints = 0;
	for (auto& seqGluepoints : tempGluepoints)
	{
		std::vector<SetPoint1d*> sortedSets;
		for (auto& gp : seqGluepoints.second) sortedSets.push_back(&gp);
		std::sort(sortedSets.begin(), sortedSets.end(),
				  [](const SetPoint1d* pt1, const SetPoint1d* pt2)
				  {return pt1->data.pos < pt2->data.pos;});
		numGluepoints += sortedSets.size();

		std::vector<SetPoint1d*> currentGroup;
		for (auto& gp : sortedSets)
		{
			if (currentGroup.empty() || 
				gp->data.pos - currentGroup.back()->data.pos < _maxSeparation)
			{
				currentGroup.push_back(gp);
			}
			else
			{
				addConsensusPoint(currentGroup);
				currentGroup.clear();
				currentGroup.push_back(gp);
			}
		}
		if (!currentGroup.empty())
		{
			addConsensusPoint(currentGroup);
		}
	}
	//add flanking points
	for (auto& seqRec : _asmSeqs.getIndex())
	{
		auto& seqPoints = _gluePoints[seqRec.first];
		seqPoints.emplace(seqPoints.begin(), pointId++, 
						  seqRec.first, 0);
		seqPoints.emplace_back(pointId++, seqRec.first, 
							   _asmSeqs.seqLen(seqRec.first) - 1);
		numGluepoints += 2;
	}
	Logger::get().debug() << "Created " << numGluepoints << " gluepoints";

	//free memory
	for (auto& seqEndpoints : endpoints)
	{
		for (auto ep : seqEndpoints.second) delete ep;
	}
}

//Cleaning up some messy tandem repeats that are actually
//artifatcs of the alignment
void RepeatGraph::collapseTandems()
{
	std::unordered_map<size_t, std::unordered_set<size_t>> tandemLefts;
	std::unordered_map<size_t, std::unordered_set<size_t>> tandemRights;

	for (auto& seqPoints : _gluePoints)
	{
		size_t leftId = 0;
		size_t rightId = 0;
		while (rightId < seqPoints.second.size())
		{
			while (rightId < seqPoints.second.size() && 
				   seqPoints.second[leftId].pointId == 
				   		seqPoints.second[rightId].pointId)
			{
				++rightId;
			}
			if (rightId - leftId > 1)
			{
				size_t tandemId = seqPoints.second[leftId].pointId;
				if (rightId < seqPoints.second.size())
				{
					tandemRights[tandemId]
						.insert(seqPoints.second[rightId].pointId);
				}
				if (leftId > 0)
				{
					tandemLefts[tandemId]
						.insert(seqPoints.second[leftId - 1].pointId);
				}
			}
			leftId = rightId;
		}
	}

	/*for (auto& ptLefts : tandemLefts)
	{
		Logger::get().debug() << "Lefts of " << ptLefts.first;
		for (auto& left : ptLefts.second)
		{
			Logger::get().debug() << "\t" << left;
		}
		Logger::get().debug() << "Rights of " << ptLefts.first;
		for (auto& right : tandemRights[ptLefts.first])
		{
			Logger::get().debug() << "\t" << right;
		}
	}*/

	for (auto& seqPoints : _gluePoints)
	{
		std::vector<GluePoint> newPoints;
		size_t leftId = 0;
		size_t rightId = 0;
		while (rightId < seqPoints.second.size())
		{
			while (rightId < seqPoints.second.size() && 
				   seqPoints.second[leftId].pointId == 
				   		seqPoints.second[rightId].pointId)
			{
				++rightId;
			}

			//int32_t span = seqPoints.second[rightId - 1].position - 
			//			   seqPoints.second[leftId].position;
			//if (rightId - leftId == 1 || span > Parameters::get().minimumOverlap)
			if (rightId - leftId == 1)
			{
				for (size_t i = leftId; i < rightId; ++i)
				{
					newPoints.push_back(seqPoints.second[i]);
				}
			}
			else	//see if we can collapse this tandem repeat
			{
				size_t tandemId = seqPoints.second[leftId].pointId;
				bool leftDetemined = tandemLefts[tandemId].size() == 1;
				bool rightDetemined = tandemRights[tandemId].size() == 1;
				if (!leftDetemined && !rightDetemined)
				{
					for (size_t i = leftId; i < rightId; ++i)
					{
						newPoints.push_back(seqPoints.second[i]);
					}
				}
				else if(leftDetemined && !rightDetemined)
				{
					newPoints.push_back(seqPoints.second[rightId - 1]);
				}
				else if (!leftDetemined && rightDetemined)
				{
					newPoints.push_back(seqPoints.second[leftId]);
				}
				//if both determined, just don't add this guy
			}
			leftId = rightId;
		}
		seqPoints.second = newPoints;
	}
}

void RepeatGraph::initializeEdges(const OverlapContainer& asmOverlaps)
{
	Logger::get().debug() << "Initializing edges";

	typedef std::pair<GraphNode*, GraphNode*> NodePair;
	std::unordered_map<NodePair, std::vector<SequenceSegment>, 
					   pairhash> parallelSegments;
	std::unordered_map<NodePair, NodePair, pairhash> complEdges;

	std::unordered_map<size_t, GraphNode*> nodeIndex;
	auto idToNode = [&nodeIndex, this](size_t nodeId)
	{
		if (!nodeIndex.count(nodeId))
		{
			nodeIndex[nodeId] = this->addNode();
		}
		return nodeIndex[nodeId];
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

			GraphNode* leftNode = idToNode(gpLeft.pointId);
			GraphNode* rightNode = idToNode(gpRight.pointId);
			NodePair fwdPair = std::make_pair(leftNode, rightNode);

			GraphNode* complLeftNode = idToNode(complLeft.pointId);
			GraphNode* complRightNode = idToNode(complRight.pointId);
			NodePair revPair = std::make_pair(complLeftNode, complRightNode);

			int32_t seqLen = _asmSeqs.seqLen(gpLeft.seqId);
			parallelSegments[fwdPair].emplace_back(gpLeft.seqId, seqLen, 
												   gpLeft.position, 
							  					   gpRight.position);
			parallelSegments[revPair]
				.push_back(parallelSegments[fwdPair].back().complement());

			complEdges[fwdPair] = revPair;
			complEdges[revPair] = fwdPair;
		}
	}

	auto segIntersect = [] (const SequenceSegment& s, const OverlapRange& o)
	{
		return std::min(o.curEnd, s.end) - std::max(o.curBegin, s.start);
	};

	std::unordered_set<NodePair, pairhash> usedPairs;
	for (auto& nodePairSeqs : parallelSegments)
	{
		if (usedPairs.count(nodePairSeqs.first)) continue;
		usedPairs.insert(complEdges[nodePairSeqs.first]);

		//cluster segments based on their overlaps
		std::vector<SetNode<SequenceSegment*>*> segmentsClusters;
		for (auto& seg : nodePairSeqs.second) 
		{
			segmentsClusters.push_back(new SetNode<SequenceSegment*>(&seg));
		}
		for (auto& segOne : segmentsClusters)
		{
			for (auto& segTwo : segmentsClusters)
			{
				if (findSet(segOne) == findSet(segTwo)) continue;

				auto& overlaps = asmOverlaps.getOverlapIndex()
												.at(segOne->data->seqId);
				for (auto& ovlp : overlaps)
				{
					if (ovlp.extId != segTwo->data->seqId) continue;

					int32_t intersectOne = segIntersect(*segOne->data, ovlp);
					int32_t intersectTwo = segIntersect(*segTwo->data, 
														ovlp.reverse());
					float rateOne = (float)intersectOne / 
						(segOne->data->end - segOne->data->start);
					float rateTwo = (float)intersectTwo / 
						(segTwo->data->end - segTwo->data->start);

					if (rateOne > 0.5 && rateTwo > 0.5)
						//abs(intersectOne - intersectTwo) < _maxSeparation)
						//abs(intersectOne - intersectTwo) < 
						//	std::max(intersectOne, intersectTwo) / 2)
					{
						unionSet(segOne, segTwo);
						break;
					}
				}
			}
		}
		std::unordered_map<SetNode<SequenceSegment*>*, 
						   std::vector<SequenceSegment*>> edgeClusters;
		for (auto& setNode : segmentsClusters)
		{
			edgeClusters[findSet(setNode)].push_back(setNode->data);
		}
		//

		//add edge foe each cluster
		std::vector<SequenceSegment> usedSegments;
		for (auto& edgeClust : edgeClusters)
		{
			//in case we have complement edges within the node pair
			auto& anySegment = *edgeClust.second.front();
			if (std::find(usedSegments.begin(), usedSegments.end(), anySegment) 
						  != usedSegments.end()) continue;

			GraphNode* leftNode = nodePairSeqs.first.first;
			GraphNode* rightNode = nodePairSeqs.first.second;
			GraphEdge newEdge(leftNode, rightNode, FastaRecord::Id(_nextEdgeId));
			for (auto& seg : edgeClust.second)
			{
				newEdge.seqSegments.push_back(*seg);
				usedSegments.push_back(seg->complement());
			}

			//check if it's self-complmenet
			bool selfComplement = std::find(usedSegments.begin(), 
						usedSegments.end(), anySegment) != usedSegments.end();
			newEdge.selfComplement = selfComplement;

			this->addEdge(std::move(newEdge));
			if (!selfComplement)
			{
				leftNode = complEdges[nodePairSeqs.first].first;
				rightNode = complEdges[nodePairSeqs.first].second;
				GraphEdge* complEdge = this->addEdge(GraphEdge(leftNode, rightNode, 
												FastaRecord::Id(_nextEdgeId + 1)));
				for (auto& seg : edgeClust.second)
				{
					complEdge->seqSegments.push_back(seg->complement());
				}
			}

			_nextEdgeId += 2;
		}

		for (auto& s : segmentsClusters) delete s;
	}

	this->logEdges();
}


void RepeatGraph::logEdges()
{
	typedef std::pair<SequenceSegment*, GraphEdge*> SegEdgePair;
	std::unordered_map<FastaRecord::Id, 
					   std::vector<SegEdgePair>> sequenceEdges;
	for (auto& edge : this->iterEdges())
	{
		for (auto& segment : edge->seqSegments)
		{
			sequenceEdges[segment.seqId].push_back({&segment, edge});
		}
	}
	for (auto& seqEdgesPair : sequenceEdges)
	{
		std::sort(seqEdgesPair.second.begin(), seqEdgesPair.second.end(),
				  [](const SegEdgePair& s1, const SegEdgePair& s2)
				  	{return s1.first->start < s2.first->start;});
	}
	for (auto& seqEdgesPair : sequenceEdges)
	{
		if (!seqEdgesPair.first.strand()) continue;

		for (size_t i = 0; i < seqEdgesPair.second.size(); ++i)
		{
			SequenceSegment* segment = seqEdgesPair.second[i].first;
			GraphEdge* edge = seqEdgesPair.second[i].second;

			std::string unique = edge->seqSegments.size() == 1 ? "*" : " ";
			Logger::get().debug() << unique << "\t" 
								  << edge->edgeId.signedId() << "\t" 
								  << _asmSeqs.seqName(segment->seqId) << "\t"
								  << segment->start << "\t" 
								  << segment->end << "\t"
								  << segment->end - segment->start;
		}
	}
	Logger::get().debug() << "Total edges: " << _nextEdgeId / 2;
}


GraphPath RepeatGraph::complementPath(const GraphPath& path)
{
	GraphPath complEdges;
	for (auto itEdge = path.rbegin(); itEdge != path.rend(); ++itEdge)
	{
		complEdges.push_back(_idToEdge.at((*itEdge)->edgeId.rc()));
	}

	return complEdges;
}

GraphEdge* RepeatGraph::complementEdge(GraphEdge* edge)
{
	return _idToEdge.at(edge->edgeId.rc());
}

GraphNode* RepeatGraph::complementNode(GraphNode* node)
{
	if (!node->outEdges.empty())
	{
		return this->complementEdge(node->outEdges.front())->nodeRight;
	}
	else if(!node->inEdges.empty())
	{
		return this->complementEdge(node->inEdges.front())->nodeLeft;
	}
	return nullptr;
}
