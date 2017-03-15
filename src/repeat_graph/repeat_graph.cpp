//(c) 2016 by Authors
//This file is a part of ABruijn program.
//Released under the BSD license (see LICENSE file)


#include "../sequence/overlap.h"
#include "../sequence/vertex_index.h"
#include "../sequence/config.h"
#include "utils.h"

#include "repeat_graph.h"
#include "disjoint_set.h"

#include <deque>

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


GraphPath RepeatGraph::complementPath(const GraphPath& path)
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
