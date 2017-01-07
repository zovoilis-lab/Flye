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

struct GluePoint
{
	GluePoint(size_t id = 0, FastaRecord::Id seqId = FastaRecord::ID_NONE,
			  int32_t position = 0):
		pointId(id), seqId(seqId), position(position) {}

	size_t 	pointId;
	FastaRecord::Id seqId;
	int32_t	position;
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
									Constants::maximumOverhang);
	OverlapContainer asmOverlaps(asmOverlapper, _asmSeqs);
	asmOverlaps.findAllOverlaps();
	
	//cluster interval endpoints
	typedef SetNode<Point> Endpoint;
	std::vector<Endpoint> endpoints;
	for (auto& seqOvlps : asmOverlaps.getOverlapIndex())
	{
		for (auto& ovlp : seqOvlps.second)
		{
			endpoints.emplace_back(Point(ovlp.curId, ovlp.curBegin, 
										 ovlp.extId, ovlp.extBegin));
			endpoints.emplace_back(Point(ovlp.curId, ovlp.curEnd, 
										 ovlp.extId, ovlp.extEnd));
		}
	}

	for (auto& p1 : endpoints)
	{
		for (auto& p2 : endpoints)
		{
			if (p1.data.curId == p2.data.curId &&
				abs(p1.data.curPos - p2.data.curPos) <= 300)
			{
				unionSet(&p1, &p2);
			}
		}
	}

	//for each cluster, get a consensus position and 
	//the corresponding Y coordinates
	std::unordered_map<Endpoint*, std::vector<Endpoint*>> clusters;
	for (auto& endpoint : endpoints)
	{
		clusters[findSet(&endpoint)].push_back(&endpoint);
	}

	std::vector<GluePoint> gluePoints;

	//for each cluster:
	for (auto& clustEndpoints : clusters)
	{
		int32_t sum = 0;
		for (auto& ep : clustEndpoints.second)
		{
			sum += ep->data.curPos;
		}
		int32_t clusterXpos = sum / clustEndpoints.second.size();
		FastaRecord::Id clustSeq = clustEndpoints.second.front()->data.curId;

		bool processed = true;
		for (auto& gp : gluePoints)
		{
			if (clustSeq == gp.seqId && abs(clusterXpos - gp.position) < 300)
			{
				processed = true;
				break;
			}
		}
		if (processed) break;

		gluePoints.emplace_back(gluePoints.size(), clustSeq, clusterXpos);

		std::vector<Endpoint> extCoords;
		for (auto& ep : clustEndpoints.second)
		{
			extCoords.push_back(*ep);
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
					abs(p1.data.extPos - p2.data.extPos) <= 300)
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
			int32_t sum = 0;
			for (auto& ep : extClust.second)
			{
				sum += ep->data.extPos;
			}
			int32_t clusterYpos = sum / extClust.second.size();
			FastaRecord::Id extSeq = extClust.second.front()->data.extId;

			bool processed = true;
			for (auto& gp : gluePoints)
			{
				if (extSeq == gp.seqId && abs(clusterYpos - gp.position) < 300)
				{
					processed = true;
					break;
				}
			}
			if (processed) break;

			gluePoints.emplace_back(gluePoints.size(), extSeq, clusterYpos);
		}
	}
}
