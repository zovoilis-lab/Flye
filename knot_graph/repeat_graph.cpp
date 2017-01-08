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
	const int32_t THD = 1500;
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

	//get all repetitive regions
	for (auto& seqOvlps : asmOverlaps.getOverlapIndex())
	{
		typedef SetNode<const OverlapRange*> OvlpSet;
		std::list<OvlpSet> ovlpSet;
		for (auto& ovlp : seqOvlps.second) ovlpSet.emplace_back(&ovlp);

		for (auto& nodeOne : ovlpSet)
		{
			for (auto& nodeTwo : ovlpSet)
			{
				if (nodeOne.data->curIntersect(*nodeTwo.data) > 0)
				{
					unionSet(&nodeOne, &nodeTwo);
				}
			}
		}
		std::unordered_map<OvlpSet*, std::vector<OvlpSet*>> clusters;
		for (auto& node : ovlpSet)
		{
			clusters[findSet(&node)].push_back(&node);
		}

		for (auto& clustNodes : clusters)
		{
			int32_t clustStart = std::numeric_limits<int32_t>::max();
			int32_t clustEnd = 0;

			for (auto& node : clustNodes.second)
			{
				clustStart = std::min(clustStart, node->data->curBegin);
				clustEnd = std::max(clustEnd, node->data->curEnd);
			}

			_repetitiveRegions[seqOvlps.first]
					.emplace_back(clustStart, clustEnd);
		}
	}
	//
	
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
				abs(p1.data.curPos - p2.data.curPos) <= THD)
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
					abs(p1.data.extPos - p2.data.extPos) <= THD)
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
		//add all cluster-specific points into the shared set
		std::vector<GlueSet*> toMerge;
		for (auto& clustPt : clusterPoints)
		{
			bool used = false;
			for (auto& glueNode : glueSets)
			{
				if (glueNode.data.seqId == clustPt.seqId &&
					abs(glueNode.data.pos - clustPt.position) < THD)
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


	for (auto& seqPoints : _gluePoints)
	{
		for (auto& pt : seqPoints.second)
		{
			Logger::get().debug() << pt.seqId << "\t" 
								  << pt.position << "\t" << pt.pointId;
		}
	}
}

void RepeatGraph::outputDot(const std::string& filename)
{
	std::ofstream fout(filename);
	fout << "digraph {\n";

	for (auto& seqEdgesPair : _gluePoints)
	{
		for (size_t i = 0; i < seqEdgesPair.second.size() - 1; ++i)
		{
			GluePoint gpLeft = seqEdgesPair.second[i];
			GluePoint gpRight = seqEdgesPair.second[i + 1];

			bool repetitive = false;
			for (auto& interval : _repetitiveRegions[gpLeft.seqId])
			{
				if (gpLeft.position >= interval.first - 1500 && 
					gpRight.position <= interval.second + 1500)
				{
					repetitive = true;
					break;
				}
			}
			std::string color = repetitive ? "red" : "black";

			//fout << "\"" << gpLeft.pointId << "\" -> \"" << gpRight.pointId
			//	 << "\" [label = \"" << _asmSeqs.seqName(gpLeft.seqId) 
			//	 << "[" << gpLeft.position << ":" 
			//	 << gpRight.position << "]\"] ;\n";
			fout << "\"" << gpLeft.pointId << "\" -> \"" << gpRight.pointId
				 << "\" [label = \"" << gpRight.position - gpLeft.position
				 << "\", color = \"" << color << "\"] ;\n";
		}
	}
	fout << "}\n";
}

