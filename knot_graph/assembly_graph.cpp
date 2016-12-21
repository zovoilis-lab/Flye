
#include "assembly_graph.h"
#include "config.h"
#include "logger.h"
#include "disjoint_set.h"
#include <fstream>

namespace
{
	struct SeqKnot
	{
		SeqKnot(int32_t start, int32_t end, Knot::Id knot):
			start(start), end(end), knot(knot) {}
		int32_t start;
		int32_t end;
		Knot::Id knot;
	};
}

void AssemblyGraph::construct(const OverlapContainer& ovlpContainer)
{
	const int OVLP_THR = 0;

	//forming overlap-based clusters
	typedef SetNode<OverlapRange> DjsOverlap;
	std::unordered_map<FastaRecord::Id, 
					   std::list<DjsOverlap>> overlapClusters;
	for (auto& ovlpHash : ovlpContainer.getOverlapIndex())
	{
		for (auto& ovlp : ovlpHash.second)
		{
			overlapClusters[ovlp.curId].emplace_back(ovlp);
			overlapClusters[ovlp.extId].emplace_back(ovlp.reverse());
			unionSet(&overlapClusters[ovlp.curId].back(), 
					 &overlapClusters[ovlp.extId].back());
		}
	}
	for (auto& ovlpHash : overlapClusters)
	{
		for (auto& ovlp1 : ovlpHash.second)
		{
			for (auto& ovlp2 : ovlpHash.second)
			{
				if (ovlp1.data.curIntersect(ovlp2.data) > OVLP_THR)
				{
					unionSet(&ovlp1, &ovlp2);
				}
			}
		}
	}

	//getting cluster assignments
	std::unordered_map<DjsOverlap*, Knot::Id> knotMappings;
	std::unordered_map<DjsOverlap*, Knot::Id> reprToKnot;
	size_t nextKnotId = SEQ_END + 1;

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
				if (ovlp1.data.curIntersect(ovlp2.data) > OVLP_THR)
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

		std::vector<SeqKnot> seqKnots;
		for (auto clustHash : intersectClusters)
		{
			int32_t clustStart = std::numeric_limits<int32_t>::max();
			int32_t clustEnd = std::numeric_limits<int32_t>::min();
			for (auto& ovlp : clustHash.second)
			{
				clustStart = std::min(clustStart, ovlp->data.curBegin);
				clustEnd = std::max(clustEnd, ovlp->data.curEnd);
			}

			auto knotId = knotMappings[clustHash.second.front()];
			seqKnots.emplace_back(clustStart, clustEnd, knotId);
		}
		///

		//constructing the graph
		std::sort(seqKnots.begin(), seqKnots.end(), 
				  [](const SeqKnot& k1, const SeqKnot& k2) 
				  		{return k1.start < k2.start;});
		auto& seqEdges = _edges[seqClusterPair.first];

		seqEdges.push_back(Edge(seqClusterPair.first, 
							    0, seqKnots.front().start,
							    SEQ_BEGIN, seqKnots.front().knot));
		_knots[SEQ_BEGIN].outEdges.push_back(&seqEdges.back());
		_knots[seqKnots.front().knot].inEdges.push_back(&seqEdges.back());

		for (size_t i = 0; i < seqKnots.size() - 1; ++i)
		{
			seqEdges.push_back(Edge(seqClusterPair.first, 
							   		seqKnots[i].end, seqKnots[i + 1].start,
							   		seqKnots[i].knot, seqKnots[i + 1].knot));
			_knots[seqKnots[i].knot].outEdges.push_back(&seqEdges.back());
			_knots[seqKnots[i + 1].knot].inEdges.push_back(&seqEdges.back());
		}

		auto seqLength = _seqAssembly.seqLen(seqClusterPair.first);
		seqEdges.push_back(Edge(seqClusterPair.first, 
						   		seqKnots.back().end, seqLength,
						   		seqKnots.back().knot, SEQ_END));
		_knots[seqKnots.back().knot].outEdges.push_back(&seqEdges.back());
		_knots[SEQ_END].inEdges.push_back(&seqEdges.back());
		////
	}
}

void AssemblyGraph::outputDot(const std::string& filename)
{
	std::ofstream fout(filename);
	fout << "digraph {\n";
	for (Knot& knot : _knots)
	{
		if (knot.knotId > SEQ_END)
		{
			fout << "\"" << knot.knotId << "\"[color = red, label = repeat];\n";
		}
	}

	size_t nextNewNode = _knots.size();
	for (auto& seqEdgesPair : _edges)
	{
		for (auto& edge : seqEdgesPair.second)
		{
			size_t firstNode = edge.knotBegin;
			if (firstNode == SEQ_BEGIN) firstNode = nextNewNode++;
			size_t lastNode = edge.knotEnd;
			if (lastNode == SEQ_END) lastNode = nextNewNode++;

			fout << "\"" << firstNode << "\" -> \"" << lastNode
				 << "\" [label = \"" << _seqAssembly.seqName(edge.seqId) 
				 << "[" << edge.seqBegin << ":" << edge.seqEnd << "]\"] ;\n";
		}
	}
	fout << "}\n";
}

namespace
{
	struct EdgePos
	{
		EdgePos(Edge* edge, int32_t pos): edge(edge), position(pos) {}
		EdgePos(): edge(nullptr), position(0) {}

		Edge* edge;
		int32_t position;
	};
}

void AssemblyGraph::untangle(const OverlapContainer& ovlpContainer)
{
	std::vector<Edge*> inEdges;
	std::vector<Edge*> outEdges;
	const int OVERHANG = 500;

	for (auto& knot : _knots)
	{
		if (knot.knotId <= SEQ_END) continue;

		for (Edge* edge : knot.inEdges)
		{
			inEdges.push_back(edge);
		}
		for (Edge* edge : knot.outEdges)
		{
			outEdges.push_back(edge);
		}
	}

	int numConnections = 0;
	for (auto seqOvelaps : ovlpContainer.getOverlapIndex())
	{
		std::unordered_set<Edge*> inSupport;
		std::unordered_set<Edge*> outSupport;
		for (auto& ovlp : seqOvelaps.second)
		{
			for (Edge* edge : inEdges)
			{
				if (ovlp.extId != edge->seqId) continue;
				
				int32_t intersect = std::min(ovlp.extEnd, edge->seqEnd) -
									std::max(ovlp.extBegin, edge->seqBegin);

				//TODO: overhangs wrt to unique edges
				//TODO: agreement with the graph structure
				if (intersect < 0 ||
					edge->seqEnd - ovlp.extEnd > OVERHANG ||
					edge->seqEnd - ovlp.extBegin < OVERHANG ||
					std::min(ovlp.curBegin, edge->seqBegin) > OVERHANG) 
				{
					continue;
				}

				inSupport.insert(edge);
				//std::cout << "IN: " << _seqAssembly.seqName(edge->seqId)
				//		  << " " << edge->seqEnd
				//		  << " " << edge->seqEnd - ovlp.extBegin << std::endl;
			}

			for (Edge* edge : outEdges)
			{
				if (ovlp.extId != edge->seqId) continue;
				
				int32_t intersect = std::min(ovlp.extEnd, edge->seqEnd) -
									std::max(ovlp.extBegin, edge->seqBegin);

				int32_t readOvhg = _seqReads.seqLen(ovlp.curId) - ovlp.curEnd;
				int32_t asmOvhg = _seqAssembly.seqLen(ovlp.extId) - ovlp.extEnd;
				if (intersect < 0 ||
					ovlp.extBegin - edge->seqBegin > OVERHANG ||
					ovlp.extEnd - edge->seqBegin < OVERHANG ||
					std::min(readOvhg, asmOvhg) > OVERHANG)
				{
					continue;
				}

				outSupport.insert(edge);
				//std::cout << "OUT: " << _seqAssembly.seqName(edge->seqId)
				//		  << " " << edge->seqEnd
				//		  << " " << ovlp.extEnd - edge->seqBegin << std::endl;
			}
		}

		if (inSupport.size() == 1 && outSupport.size() == 1)
		{
			Edge* leftEdge = *inSupport.begin();
			Edge* rightEdge = *outSupport.begin();

			if (leftEdge->knotEnd != rightEdge->knotBegin)
			{
				++numConnections;

				std::cout << _seqReads.seqName(seqOvelaps.first) << std::endl;
				std::cout << "Connection: " << _seqAssembly.seqName(leftEdge->seqId)
						  << " " << leftEdge->seqEnd << " "
						  << _seqAssembly.seqName(rightEdge->seqId) << " " 
						  << rightEdge->seqBegin << std::endl;
			}
		}
	}
	std::cout << numConnections << std::endl;
}
