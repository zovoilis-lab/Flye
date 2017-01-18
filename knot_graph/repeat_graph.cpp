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
	this->buildGraph(asmOverlaps);
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
			unionSet(&overlapClusters[ovlp.curId].back(),
					 &overlapClusters[ovlp.extId].back());

			//reverse complement
			int32_t curLen = _asmSeqs.seqLen(ovlp.curId);
			int32_t extLen = _asmSeqs.seqLen(ovlp.extId);
			OverlapRange complOvlp = ovlp.complement(curLen, extLen);
			//

			overlapClusters[complOvlp.curId].emplace_back(complOvlp);
			overlapClusters[complOvlp.extId].emplace_back(complOvlp.reverse());
			unionSet(&overlapClusters[complOvlp.curId].back(),
					 &overlapClusters[complOvlp.extId].back());
			//
		}
	}
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
	}

	for (auto& seqClusterPair : overlapClusters)
	{
		//forming intersection-based clusters
		for (auto& ovlp : seqClusterPair.second)
		{
			ovlp.parent = &ovlp;
			ovlp.rank = 0;
		}

		for (auto& ovlp1 : seqClusterPair.second)
		{
			for (auto& ovlp2 : seqClusterPair.second)
			{
				if (ovlp1.data.curIntersect(ovlp2.data) > _maxSeparation)
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

			auto knotId = knotMappings[clustHash.second.front()];
			_repeatClusters[seqClusterPair.first]
					.emplace_back(seqClusterPair.first, knotId,
								  clustStart, clustEnd);
		}
	}
}

void RepeatGraph::buildGraph(const OverlapContainer& asmOverlaps)
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
}

void RepeatGraph::initializeEdges()
{
	size_t edgeId = 0;

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

			bool repetitive = false;
			size_t fwdClust = 0;
			size_t revClust = 0;

			float maxOvlp = 0;
			for (auto& cluster : _repeatClusters[gpLeft.seqId])
			{
				int32_t overlap = std::min(gpRight.position, cluster.end) -
								  std::max(gpLeft.position, cluster.start);
				float rate = float(overlap) / (gpRight.position - 
											   gpLeft.position);

				maxOvlp = std::max(maxOvlp, rate);
				if (rate > 0.5)
				{
					repetitive = true;
					fwdClust = cluster.clusterId;
					break;
				}
			}
			for (auto& cluster : _repeatClusters[complLeft.seqId])
			{
				int32_t overlap = std::min(complRight.position, cluster.end) -
								  std::max(complLeft.position, cluster.start);
				float rate = float(overlap) / (complRight.position - 
											   complLeft.position);
				if (rate > 0.5)
				{
					revClust = cluster.clusterId;
					break;
				}
			}

			_graphEdges[seqEdgesPair.first]
					.push_back(GraphEdge(gpLeft, gpRight, repetitive, 
										 FastaRecord::Id(edgeId), fwdClust));
			_graphEdges[complId]
					.push_back(GraphEdge(complLeft, complRight, repetitive, 
										 FastaRecord::Id(edgeId + 1), 
										 revClust));
			edgeId += 2;
		}
	}
	
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
	}
}


void RepeatGraph::chainReadAlignments(const SequenceContainer& edgeSeqs,
									  std::vector<EdgeAlignment> ovlps)
{
	const int32_t JUMP = 1500;

	std::sort(ovlps.begin(), ovlps.end(),
			  [](const EdgeAlignment& e1, const EdgeAlignment& e2)
			  	{return e1.overlap.curBegin < e2.overlap.curBegin;});

	typedef std::vector<EdgeAlignment*> Chain;

	std::list<Chain> activeChains;
	//Logger::get().debug() << ovlps.size();
	for (auto& edgeAlignment : ovlps)
	{
		std::list<Chain> newChains;
		for (auto& chain : activeChains)
		{
			//can extend?
			int32_t readDiff = edgeAlignment.overlap.curBegin - 
						   	   chain.back()->overlap.curEnd;
			int32_t graphDiff = edgeAlignment.overlap.extBegin +
								edgeSeqs.seqLen(chain.back()->overlap.extId) -
						   	    chain.back()->overlap.extEnd;

			if (JUMP > readDiff && readDiff > 0 &&
				JUMP > graphDiff && graphDiff > 0 &&
				chain.back()->edge->gpRight.pointId == 
				edgeAlignment.edge->gpLeft.pointId)
			{
				newChains.push_back(chain);
				chain.push_back(&edgeAlignment);
			}
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
		int32_t chainSpan = chain.back()->overlap.curEnd - 
				chain.front()->overlap.curBegin;

		int32_t overhang = _readSeqs.seqLen(chain.back()->overlap.curId) -
							chain.back()->overlap.curEnd;
		if (chainSpan > maxSpan && overhang < Constants::maximumOverhang)
		{
			maxSpan = chainSpan;
			maxChain = &chain;
		}
	}

	if (maxChain && maxChain->size() > 1)
	{
		Logger::get().debug() << _readSeqs.seqName(ovlps.front().overlap.curId);
		
		for (auto& edge : *maxChain)
		{
			Logger::get().debug() << edge->edge->gpLeft.seqId << "\t" 
								  << edge->edge->gpLeft.position << "\t"
								  << edge->edge->gpRight.position << "\t"
								  << edge->edge->repetitive << "\t"
								  << edge->edge->edgeId << "\t"
								  << edge->overlap.curBegin << "\t" 
								  << edge->overlap.curEnd;
		}
	}
}


void RepeatGraph::resolveRepeats()
{
	std::unordered_map<FastaRecord::Id, GraphEdge*> idToEdge;
	SequenceContainer pathsContainer;

	for (auto& seqEdgesPair : _graphEdges)
	{
		for (auto& edge : seqEdgesPair.second)
		{
			
			size_t len = edge.gpRight.position - edge.gpLeft.position;
			std::string sequence = _asmSeqs.getSeq(seqEdgesPair.first)
										.substr(edge.gpLeft.position, len);
			pathsContainer.addSequence(FastaRecord(sequence, "", edge.edgeId));
			idToEdge[edge.edgeId] = &edge;
		}
	}

	VertexIndex pathsIndex(pathsContainer);
	pathsIndex.countKmers(1);
	pathsIndex.buildIndex(1, 5000, 1);
	OverlapDetector readsOverlapper(pathsContainer, pathsIndex, 
									500, 1000, 0);
	OverlapContainer readsContainer(readsOverlapper, _readSeqs);
	readsContainer.findAllOverlaps();
	
	for (auto& seqOverlaps : readsContainer.getOverlapIndex())
	{
		std::vector<EdgeAlignment> alignments;
		for (auto& ovlp : filterOvlp(seqOverlaps.second))
		{
			alignments.push_back({ovlp, idToEdge[ovlp.extId]});
		}
		this->chainReadAlignments(pathsContainer, alignments);

		/*
		std::unordered_map<GraphEdge*, 
						   std::vector<const OverlapRange*>> byEdge;
		for (auto& ovlp : seqOverlaps.second)
		{
			byEdge[idToEdge[ovlp.extId]].push_back(&ovlp);
		}

		int uniqueMatches = 0;
		std::vector<OverlapRange> chosenOverlaps;
		for (auto& edgeAlignments : byEdge)
		{
			OverlapRange maxOverlap;
			for (auto& ovlp : edgeAlignments.second)
			{
				if (maxOverlap.curRange() < ovlp->curRange())
				{
					maxOverlap = *ovlp;
				}
			}
			if (!edgeAlignments.first->repetitive)
			{
				if (maxOverlap.extBegin > 500 && 
					pathsContainer.seqLen(maxOverlap.extId) - 
										maxOverlap.extEnd > 500)
				{
					continue;
				}
			}

			chosenOverlaps.push_back(maxOverlap);
			if (!edgeAlignments.first->repetitive)
			{
				++uniqueMatches;
			}
		}

		if (uniqueMatches < 2) continue;

		Logger::get().debug() << _readSeqs.seqName(seqOverlaps.first);
		for (auto& ovlp : chosenOverlaps)
		{
			GraphEdge& edge = *idToEdge[ovlp.extId];
			Logger::get().debug() << edge.gpLeft.seqId << "\t" 
								  << edge.gpLeft.position << "\t"
								  << edge.gpRight.position << "\t"
								  << edge.repetitive << "\t"
								  << edge.edgeId << "\t"
								  << ovlp.curBegin << "\t" << ovlp.curEnd;
		}*/
	}
}

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
}

void RepeatGraph::outputDot(const std::string& filename, bool collapseRepeats)
{
	const std::string COLORS[] = {"red", "darkgreen", "blue", "goldenrod", 
								  "cadetblue",
								  "darkorchid", "aquamarine1", "darkgoldenrod1",
								  "deepskyblue1", "darkolivegreen3"};

	std::ofstream fout(filename);
	fout << "digraph {\n";

	std::unordered_map<std::pair<size_t, size_t>, 
					   std::vector<int32_t>, pairhash> edgesLengths;
	std::unordered_map<std::pair<size_t, size_t>, size_t, pairhash> edgesIds;

	for (auto& seqEdgesPair : _graphEdges)
	{
		for (auto& edge : seqEdgesPair.second)
		{
			if (!edge.repetitive)
			{
				fout << "\"" << edge.gpLeft.pointId << "\" -> \"" 
					 << edge.gpRight.pointId
					 << "\" [label = \"" << edge.edgeId << " "
					 << edge.gpRight.position - edge.gpLeft.position
					 << "\", color = \"black\"] ;\n";
			}
			else
			{
				edgesLengths[std::make_pair(edge.gpLeft.pointId, 
											edge.gpRight.pointId)]
										 .push_back(edge.gpRight.position - 
										   			edge.gpLeft.position);
				edgesIds[std::make_pair(edge.gpLeft.pointId, 
										edge.gpRight.pointId)] = edge.clusterId;
			}
		}
	}

	for (auto& edgeClust : edgesLengths)
	{
		if (collapseRepeats)
		{
			int64_t sum = 0;
			for (auto& dist : edgeClust.second) sum += dist;
			int32_t avgLengths = sum / edgeClust.second.size();

			std::string color = COLORS[edgesIds[edgeClust.first] % 10];
			fout << "\"" << edgeClust.first.first << "\" -> \"" 
				 << edgeClust.first.second
				 << "\" [label = \"" << avgLengths << " ("
				 << edgeClust.second.size()
				 << ")\", color = \"" << color << "\", penwidth = 3] ;\n";
		}
	}

	fout << "}\n";
}
