
#include "assembly_graph.h"
#include "config.h"
#include "logger.h"
#include "disjoint_set.h"
#include <fstream>

void AssemblyGraph::construct(OverlapDetector& ovlpDetector)
{
	auto ovlpIndex = ovlpDetector.getOverlapIndex();

	//forming overlap-based clusters
	typedef SetNode<OverlapRange> DjsOverlap;
	std::unordered_map<FastaRecord::Id, 
					   std::list<DjsOverlap>> overlapClusters;
	//std::list<OverlapRange> reverseOverlaps;
	for (auto& ovlpHash : ovlpIndex)
	{
		for (auto& ovlp : ovlpHash.second)
		{
			overlapClusters[ovlp.curId].emplace_back(ovlp);
			overlapClusters[ovlp.extId].emplace_back(ovlp.reverse());
			unionSet(&overlapClusters[ovlp.curId].back(), 
					 &overlapClusters[ovlp.extId].back());
		}
	}/*
	for (auto& ovlpHash : overlapClusters)
	{
		std::sort(ovlpHash.second.begin(), ovlpHash.second.end(),
		 		  [](const SetNode<OverlapRange*>& s1, 
				  	 const SetNode<OverlapRange*>& s2)
					 	{return s1.data->curBegin < s2.data->curBegin;});
	}*/
	for (auto& ovlpHash : overlapClusters)
	{
		for (auto& ovlp1 : ovlpHash.second)
		{
			for (auto& ovlp2 : ovlpHash.second)
			{
				if (ovlp1.data.curIntersect(ovlp2.data) > 0)
				{
					unionSet(&ovlp1, &ovlp2);
				}
			}
		}
	}

	//getting cluster assignments
	std::unordered_map<DjsOverlap*, Knot::Id> knotMappings;
	std::unordered_map<DjsOverlap*, Knot::Id> reprToKnot;
	size_t nextKnotId = 0;

	for (auto& ovlpHash : overlapClusters)
	{
		for (auto& ovlp : ovlpHash.second)
		{
			if (!reprToKnot.count(findSet(&ovlp)))
			{
				_knots.push_back(Knot(nextKnotId));
				reprToKnot[findSet(&ovlp)] = nextKnotId++;
			}
			knotMappings[&ovlp] = reprToKnot[findSet(&ovlp)];
		}
	}

	Logger::get().debug() << "Number of knots: " << _knots.size();

	//forming intersection-based clusters
	for (auto& ovlpHash : overlapClusters)
	{
		for (auto& ovlp : ovlpHash.second)
		{
			ovlp.parent = &ovlp;
			ovlp.rank = 0;
		}

		for (auto& ovlp1 : ovlpHash.second)
		{
			for (auto& ovlp2 : ovlpHash.second)
			{
				if (ovlp1.data.curIntersect(ovlp2.data) > 0)
				{
					unionSet(&ovlp1, &ovlp2);
				}
			}
		}

		std::unordered_map<DjsOverlap*, 
						   std::vector<DjsOverlap*>> intersectClusters;
		for (auto& ovlp : ovlpHash.second)
		{
			intersectClusters[findSet(&ovlp)].push_back(&ovlp);
		}

		for (auto clustHash : intersectClusters)
		{
			int32_t clustStart = std::numeric_limits<int32_t>::max();
			int32_t clustEnd = std::numeric_limits<int32_t>::min();
			for (auto& ovlp : clustHash.second)
			{
				clustStart = std::min(clustStart, ovlp->data.curBegin);
				clustEnd = std::max(clustStart, ovlp->data.curEnd);
			}

			auto knotId = knotMappings[clustHash.second.front()];
			_edges[ovlpHash.first].push_back(Edge(ovlpHash.first, clustStart, 
												  clustEnd, knotId));
		}
	}
	
	for (auto& seqEdgesPair : _edges)
	{
		std::sort(seqEdgesPair.second.begin(), seqEdgesPair.second.end(), 
				  [](Edge& e1, Edge& e2){return e1.seqBegin < e2.seqBegin;});
		std::cout << seqEdgesPair.first << " " 
				  << seqEdgesPair.second.size() << std::endl;
	}
}

void AssemblyGraph::outputDot(const std::string& filename)
{
}
