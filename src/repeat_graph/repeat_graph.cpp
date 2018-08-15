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
	struct Point2d
	{
		Point2d(FastaRecord::Id curId = FastaRecord::ID_NONE, 
				int32_t curPos = 0, 
			    FastaRecord::Id extId = FastaRecord::ID_NONE, 
				int32_t extPos = 0):
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

bool GraphEdge::isTip() const
{
	return nodeLeft->inEdges.empty() || nodeRight->outEdges.empty();
}

std::unordered_set<GraphEdge*> GraphEdge::adjacentEdges()
{
	std::unordered_set<GraphEdge*> edges;
	for (auto& e: nodeLeft->inEdges) edges.insert(e);
	for (auto& e: nodeLeft->outEdges) edges.insert(e);
	for (auto& e: nodeRight->inEdges) edges.insert(e);
	for (auto& e: nodeRight->outEdges) edges.insert(e);
	edges.erase(this);
	return edges;
}

void RepeatGraph::build()
{
	//getting overlaps
	VertexIndex asmIndex(_asmSeqs, 
						 (int)Config::get("repeat_graph_kmer_sample"));
	asmIndex.countKmers(/*min freq*/ 1, /*genome size*/ 0);
	asmIndex.setRepeatCutoff(/*min freq*/ 1);
	asmIndex.buildIndex(/*min freq*/ 2);

	float badEndAdj = (float)Config::get("repeat_graph_ovlp_ident") / 2; //TODO: add to config
	OverlapDetector asmOverlapper(_asmSeqs, asmIndex, 
								  (int)Config::get("maximum_jump"), 
								  Parameters::get().minimumOverlap,
								  /*no overhang*/ 0, /*all overlaps*/ 0,
								  /*keep alignment*/ true, /*only max*/ false,
								  (float)Config::get("repeat_graph_ovlp_ident"),
								  badEndAdj);

	OverlapContainer asmOverlaps(asmOverlapper, _asmSeqs);
	asmOverlaps.findAllOverlaps();
	asmOverlaps.buildIntervalTree();
	asmOverlaps.overlapDivergenceStats();

	this->getGluepoints(asmOverlaps);
	this->collapseTandems();
	this->initializeEdges(asmOverlaps);
	//this->markChimericEdges();
}

void RepeatGraph::getGluepoints(OverlapContainer& asmOverlaps)
{
	//Process alignment ends to construct glueing points for the
	//repeat graph construction.
	//Note, that there are two kinds of clustering here:
	//based on X and Y coordinates (position + seqId) of the points,
	//that might often define different subsets
	
	Logger::get().debug() << "Computing gluepoints";
	typedef SetNode<Point2d> SetPoint2d;
	std::unordered_map<FastaRecord::Id, SetVec<Point2d>> endpoints;

	//first, extract endpoints from all overlaps.
	//each point has X and Y coordinates (curSeq and extSeq)
	for (auto& seq : _asmSeqs.iterSeqs())
	{
		for (auto& ovlp : asmOverlaps.lazySeqOverlaps(seq.id))
		{
			//if (_filteredSeqs.count(ovlp.curId) ||
			//	_filteredSeqs.count(ovlp.extId)) continue;

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
	std::vector<SetPoint2d*> flatEndpoints;
	for (auto& seqPoints : endpoints)
	{
		flatEndpoints.insert(flatEndpoints.end(), seqPoints.second.begin(), 
							 seqPoints.second.end());
	}
	auto clusters = groupBySet(flatEndpoints);

	typedef SetNode<Point1d> SetPoint1d;
	std::unordered_map<FastaRecord::Id, SetVec<Point1d>> tempGluepoints;
	std::unordered_map<SetPoint1d*, SetPoint1d*> complements;

	//we will now split each cluster based on it's Y coordinates
	//and project these subgroups to the corresponding sequences
	for (auto& clustEndpoints : clusters)
	{
		//first, simply add projections for each point from the cluster
		FastaRecord::Id clustSeq = clustEndpoints.second.front().curId;
		if (!clustSeq.strand()) continue;	//only for forward strands

		std::vector<int32_t> positions;
		for (auto& ep : clustEndpoints.second) 
		{
			positions.push_back(ep.curPos);
		}
		int32_t clusterXpos = median(positions);

		std::vector<Point1d> clusterPoints;
		clusterPoints.emplace_back(clustSeq, clusterXpos);

		SetVec<Point2d> extCoords;
		for (auto& ep : clustEndpoints.second)
		{
			extCoords.push_back(new SetPoint2d(ep));
		}
		
		//Important part: we need also add extra projections
		//for gluepoints that lie within other existing overlaps
		//(handles situations with 'repeat hierarchy', when some
		//repeats are parts of the other bigger repeats)
		for (auto& interval : asmOverlaps
				.getCoveringOverlaps(clustSeq, clusterXpos - 1, 
									 clusterXpos + 1))
		{
			if (interval.value->curEnd - clusterXpos > _maxSeparation &&
				clusterXpos - interval.value->curBegin > _maxSeparation)
			{
				int32_t projectedPos = interval.value->project(clusterXpos);
				extCoords.push_back(new SetPoint2d(Point2d(clustSeq, clusterXpos,
											   	   interval.value->extId, 
											   	   projectedPos)));
			}
		}

		//Finally, cluster the projected points based on Y coordinates
		std::sort(extCoords.begin(), extCoords.end(),
				  [](const SetPoint2d* p1, const SetPoint2d* p2)
				  {return (p1->data.extId != p2->data.extId) ? 
				  		   p1->data.extId < p2->data.extId :
						   p1->data.extPos < p2->data.extPos;});
		for (size_t i = 0; i < extCoords.size() - 1; ++i)
		{
			auto* p1 = extCoords[i];
			auto* p2 = extCoords[i + 1];
			if (p1->data.extId == p2->data.extId &&
				abs(p1->data.extPos - p2->data.extPos) < _maxSeparation)
			{
				unionSet(p1, p2);
			}

		}
		auto extClusters = groupBySet(extCoords);

		//now, get coordinates for each cluster
		for (auto& extClust : extClusters)
		{
			std::vector<int32_t> positions;
			for (auto& ep : extClust.second) 
			{
				positions.push_back(ep.extPos);
			}
			int32_t clusterYpos = median(positions);

			FastaRecord::Id extSeq = extClust.second.front().extId;
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

			//inserting into sorted vector
			auto cmp = [] (const SetPoint1d* gp, int32_t pos)
								{return gp->data.pos < pos;};
			size_t i = std::lower_bound(seqGluepoints.begin(), 
										seqGluepoints.end(),
										clustPt.pos, cmp) - seqGluepoints.begin();
			auto cmp2 = [] (int32_t pos, const SetPoint1d* gp)
								{return pos < gp->data.pos;};
			size_t ci = std::upper_bound(complGluepoints.begin(), 
										 complGluepoints.end(),
										 complPt.pos, cmp2) - complGluepoints.begin();
			if (!seqGluepoints.empty())
			{
				if (i > 0 && 
					clustPt.pos - seqGluepoints[i - 1]->data.pos < _maxSeparation)
				{
					toMerge.push_back(seqGluepoints[i - 1]);
				}
				if (i < seqGluepoints.size() && 
					seqGluepoints[i]->data.pos - clustPt.pos < _maxSeparation)
				{
					toMerge.push_back(seqGluepoints[i]);
				}
			}

			seqGluepoints.insert(seqGluepoints.begin() + i, 
								 new SetPoint1d(clustPt));
			auto fwdPtr = seqGluepoints[i];
			complGluepoints.insert(complGluepoints.begin() + ci,
								   new SetPoint1d(complPt));
			auto revPtr = complGluepoints[ci];

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
		int32_t clusterSize = group.back()->data.pos - 
							  group.front()->data.pos;

		//big cluster corresponding to a tandem repeat - 
		//split it into many short edges
		if (clusterSize > _maxSeparation)
		{
			_gluePoints[reprPoint->data.seqId]
				.emplace_back(setToId[findSet(reprPoint)],
							  reprPoint->data.seqId, group.front()->data.pos);

			int32_t repeats = std::floor(clusterSize / _maxSeparation);
			int32_t mode = clusterSize / repeats;
			for (int32_t i = 1; i < repeats; ++i)
			{
				int32_t pos = group.front()->data.pos + mode * i;
				_gluePoints[reprPoint->data.seqId]
					.emplace_back(setToId[findSet(reprPoint)],
								  reprPoint->data.seqId, pos);

			}

			_gluePoints[reprPoint->data.seqId]
				.emplace_back(setToId[findSet(reprPoint)],
							  reprPoint->data.seqId, group.back()->data.pos);
		}
		//"normal" endpoint - just take a consensus
		else
		{
			std::vector<int32_t> positions;
			for (auto& ep : group) 
			{
				positions.push_back(ep->data.pos);
			}
			int32_t clusterXpos = median(positions);

			_gluePoints[reprPoint->data.seqId]
				.emplace_back(setToId[findSet(reprPoint)],
							  reprPoint->data.seqId, clusterXpos);
		}

	};

	for (auto& seqGluepoints : tempGluepoints)
	{
		std::vector<SetPoint1d*> currentGroup;
		for (auto& gp : seqGluepoints.second)
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

	//add flanking points, if needed
	const int MAX_TIP = Parameters::get().minimumOverlap;
	for (auto& seq : _asmSeqs.iterSeqs())
	{
		if (!seq.id.strand()) continue;
		//if (_filteredSeqs.count(seq.id)) continue;
		auto& seqPoints = _gluePoints[seq.id];
		auto& complPoints = _gluePoints[seq.id.rc()];

		if (seqPoints.empty() || seqPoints.front().position > MAX_TIP)
		{
			seqPoints.emplace(seqPoints.begin(), pointId++, 
							  seq.id, 0);
			complPoints.emplace_back(pointId++, seq.id.rc(),
							  		 _asmSeqs.seqLen(seq.id) - 1);
		}
		if (seqPoints.size() == 1 || 
			_asmSeqs.seqLen(seq.id) - seqPoints.back().position > MAX_TIP)
		{
			seqPoints.emplace_back(pointId++, seq.id, 
								   _asmSeqs.seqLen(seq.id) - 1);
			complPoints.emplace(complPoints.begin(), pointId++, 
							  	seq.id.rc(), 0);
		}
	}

	//ensure coordinates are symmetric
	for (auto& seq : _asmSeqs.iterSeqs())
	{
		if (!seq.id.strand()) continue;
		//if (_filteredSeqs.count(seq.id)) continue;
		auto& seqPoints = _gluePoints[seq.id];
		auto& complPoints = _gluePoints[seq.id.rc()];

		int32_t seqLen = _asmSeqs.seqLen(seq.id);
		for (size_t i = 0; i < seqPoints.size(); ++i)
		{
			complPoints[seqPoints.size() - i - 1].position = 
				seqLen - seqPoints[i].position - 1;
		}
	}

	int numGluepoints = 0;
	for (auto& seqRec : _gluePoints) numGluepoints += seqRec.second.size();
	Logger::get().debug() << "Created " << numGluepoints << " gluepoints";
}

//Cleaning up some messy tandem repeats that are actually
//artifatcs of the alignment
void RepeatGraph::collapseTandems()
{
	std::unordered_map<size_t, std::unordered_set<size_t>> tandemLefts;
	std::unordered_map<size_t, std::unordered_set<size_t>> tandemRights;
	std::unordered_set<size_t> bigTandems;

	/*std::unordered_map<size_t, size_t> complPoints;
	for (auto& seqPoints : _gluePoints)
	{
		auto& complSeq = _gluePoints[seqPoints.first.rc()];
		if (seqPoints.second.size() != complSeq.size()) throw std::runtime_error("AAA");
		for (size_t i = 0; i < seqPoints.second.size(); ++i)
		{
			complPoints[seqPoints.second[i].pointId] =
				complSeq[seqPoints.second.size() - i - 1].pointId;
		}
	}*/

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
				if (seqPoints.second[rightId - 1].position - 
						seqPoints.second[leftId].position > 
						Parameters::get().minimumOverlap)
				{
					bigTandems.insert(tandemId);
					//bigTandems.insert(complPoints[tandemId]);
				}

				if (rightId < seqPoints.second.size())
				{
					tandemRights[tandemId]
						.insert(seqPoints.second[rightId].pointId);
				}
				else
				{
					tandemRights[tandemId].insert(-1);
				}
				if (leftId > 0)
				{
					tandemLefts[tandemId]
						.insert(seqPoints.second[leftId - 1].pointId);
				}
				else
				{
					tandemLefts[tandemId].insert(-1);
				}
			}
			leftId = rightId;
		}
	}

	int collapsedLeft = 0;
	int collapsedRight = 0;
	int collapsedBoth = 0;
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

			size_t tandemId = seqPoints.second[leftId].pointId;
			//size_t complId = complPoints[tandemId];
			if (rightId - leftId == 1 || bigTandems.count(tandemId))
			{
				for (size_t i = leftId; i < rightId; ++i)
				{
					newPoints.push_back(seqPoints.second[i]);
				}
			}
			else	//see if we can collapse this tandem repeat
			{
				//making sure graph remains symmetric
				bool leftDetemined = tandemLefts[tandemId].size() == 1;
									 //tandemRights[complId].size() == 1;
				bool rightDetemined = tandemRights[tandemId].size() == 1;
									  //tandemLefts[complId].size() == 1;
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
					++collapsedLeft;
				}
				else if (!leftDetemined && rightDetemined)
				{
					newPoints.push_back(seqPoints.second[leftId]);
					++collapsedRight;
				}
				else
				{
					int32_t newPos = (seqPoints.second[rightId - 1].position +
									  seqPoints.second[leftId].position) / 2;
					newPoints.push_back(seqPoints.second[leftId]);
					newPoints.back().position = newPos;
					++collapsedBoth;
				}
			}
			leftId = rightId;
		}
		seqPoints.second = newPoints;
	}

	Logger::get().debug() << "Tandems removed: " << collapsedLeft << " left, " 
		<< collapsedRight << " right, " << collapsedBoth << " both";
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
		if (seqEdgesPair.second.size() < 2) continue;

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

	auto segIntersect = [] (const SequenceSegment& s, int32_t intBegin, 
							int32_t intEnd)
	{
		return std::max(std::min(intEnd, s.end) - 
						std::max(intBegin, s.start), 0);
	};

	std::unordered_set<NodePair, pairhash> usedPairs;
	for (auto& nodePairSeqs : parallelSegments)
	{
		if (usedPairs.count(nodePairSeqs.first)) continue;
		usedPairs.insert(complEdges[nodePairSeqs.first]);

		//creating set and building index
		SetVec<SequenceSegment*> segmentSets;
		typedef SetNode<SequenceSegment*> SetSegment;
		std::unordered_map<FastaRecord::Id, 
						   std::vector<SetSegment*>> segmentIndex;
		for (auto& seg : nodePairSeqs.second) 
		{
			segmentSets.push_back(new SetNode<SequenceSegment*>(&seg));
			segmentIndex[seg.seqId].push_back(segmentSets.back());
		}
		for (auto& seqSegments : segmentIndex)
		{
			std::sort(seqSegments.second.begin(), seqSegments.second.end(),
					  [](const SetSegment* s1, const SetSegment* s2)
					  {return s1->data->start < s2->data->start;});
		}


		//cluster segments based on their overlaps
		for (auto& setOne : segmentSets)
		{
			for (auto& interval : asmOverlaps
						.getCoveringOverlaps(setOne->data->seqId, 
											 setOne->data->start,
											 setOne->data->end))
			{
				auto& ovlp = *interval.value;
				int32_t intersectOne = 
					segIntersect(*setOne->data, ovlp.curBegin, ovlp.curEnd);
				if (intersectOne < _maxSeparation) continue;

				auto& ss = segmentIndex[ovlp.extId];
				auto cmpBegin = [] (const SetSegment* s, int32_t pos)
								    {return s->data->start < pos;};
				auto cmpEnd = [] (const SetSegment* s, int32_t pos)
								    {return s->data->end < pos;};
				auto startRange = std::lower_bound(ss.begin(), ss.end(),
												   ovlp.extBegin, cmpEnd);
				auto endRange = std::lower_bound(ss.begin(), ss.end(),
												 ovlp.extEnd, cmpBegin);
				if (endRange != ss.end()) ++endRange;
				for (;startRange != endRange; ++startRange)
				{
					auto* setTwo = *startRange;
					if (findSet(setOne) == findSet(setTwo)) continue;

					//projecting the interval endpoints
					//(overlap might be covering the actual segment)
					int32_t projStart = ovlp.project(setOne->data->start);
					int32_t projEnd = ovlp.project(setOne->data->end);
					int32_t projIntersect =
						segIntersect(*setTwo->data, projStart, projEnd);

					if (projIntersect > setOne->data->length() / 2 && 
						projIntersect > setTwo->data->length() / 2)
					{
						unionSet(setOne, setTwo);
					}
				}
			}
		}
		auto edgeClusters = groupBySet(segmentSets);
		/*if (edgeClusters.size() > 5)
		{
			Logger::get().debug() << "Node with " << segmentSets.size() << " segments";
			Logger::get().debug() << "\tclusters: " << edgeClusters.size();
			for (auto& edgeClust : edgeClusters)
			{
				int sumLen = 0;
				for (auto s : edgeClust.second) sumLen += s->length();
				Logger::get().debug() << "\t\tcl: " << edgeClust.second.size()
					<< " " << sumLen / edgeClust.second.size();
				
				if (edgeClust.second.size() < 10)
				{
					for (auto s : edgeClust.second)
					{
						std::string ovlpCovered = " ";
						for (auto& interval : asmOverlaps
								.getCoveringOverlaps(s->seqId, 
											 s->start, s->end))
						{
							auto& ovlp = *interval.value;
							int32_t intersectOne = 
								segIntersect(*s, ovlp.curBegin, ovlp.curEnd);
							if (intersectOne > 0.75f * ovlp.curRange()) 
								ovlpCovered = "*" + std::to_string(ovlp.extId.rawId());
						}

						Logger::get().debug() << "\t\t\t" << s->seqId 
							<< " " << s->start << " " << s->length() <<
							" " << ovlpCovered;
					}
				}
				else
				{
					Logger::get().debug() << "\t\t\t...";
				}
			}
		}*/

		//add edge for each cluster
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

/*
void RepeatGraph::markChimericEdges()
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

	int extraSelfCompl = 0;
	for (auto& seqEdgesPair : sequenceEdges)
	{
		if (!seqEdgesPair.first.strand()) continue;

		for (size_t i = 1; i < seqEdgesPair.second.size() - 1; ++i)
		{
			GraphEdge* leftEdge = seqEdgesPair.second[i - 1].second;
			GraphEdge* midEdge = seqEdgesPair.second[i].second;
			GraphEdge* rightEdge = seqEdgesPair.second[i + 1].second;

			if (leftEdge->edgeId == rightEdge->edgeId.rc() &&
				leftEdge->edgeId != midEdge->edgeId &&
				rightEdge->edgeId != midEdge->edgeId &&
				midEdge->seqSegments.size() == 1 &&
				midEdge->length() < Parameters::get().minimumOverlap)
			{
				midEdge->selfComplement = true;
				this->complementEdge(midEdge)->selfComplement = true;
				++extraSelfCompl;
			}
		}
	}
	Logger::get().debug() << "Extra self-complements: " << extraSelfCompl;
}*/


GraphPath RepeatGraph::complementPath(const GraphPath& path) const
{
	GraphPath complEdges;
	for (auto itEdge = path.rbegin(); itEdge != path.rend(); ++itEdge)
	{
		complEdges.push_back(_idToEdge.at((*itEdge)->edgeId.rc()));
	}

	assert(!complEdges.empty());
	return complEdges;
}

GraphEdge* RepeatGraph::complementEdge(GraphEdge* edge) const
{
	return _idToEdge.at(edge->edgeId.rc());
}

GraphNode* RepeatGraph::complementNode(GraphNode* node) const
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
