
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
			DjsOverlap* ptrOne = &overlapClusters[ovlp.curId].back(); 
			overlapClusters[ovlp.extId].emplace_back(ovlp.reverse());
			DjsOverlap* ptrTwo = &overlapClusters[ovlp.extId].back(); 
			unionSet(ptrOne, ptrTwo);
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
	size_t nextKnotId = _knots.size();

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

	Logger::get().debug() << "Number of knots: " << _knots.size() - 2;

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
		//std::cout << seqClusterPair.first << " " << intersectClusters.size() << std::endl;

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
		if (knot.knotId > SEQ_END && 
			(!knot.inEdges.empty() || !knot.outEdges.empty()))
		{
			//fout << "\"" << knot.knotId << "\"[color = red, label = repeat];\n";
			fout << "\"" << knot.knotId << "\" -> \"" << -knot.knotId
				 << "\" [label = \"repeat\", color = \"red\"] ;\n";
		}
	}

	size_t nextNewNode = _knots.size();
	for (auto& seqEdgesPair : _edges)
	{
		for (auto& edge : seqEdgesPair.second)
		{
			Knot::Id firstNode = -edge.knotBegin;
			if (firstNode == SEQ_BEGIN) firstNode = nextNewNode++;
			Knot::Id lastNode = edge.knotEnd;
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
	
	template<class T>
	void vecRemove(std::vector<T>& v, T val)
	{
		v.erase(std::remove(v.begin(), v.end(), val), v.end()); 
	}
}

void AssemblyGraph::untangle(const OverlapContainer& ovlpContainer)
{
	auto connections = this->getConnections(ovlpContainer);

	std::unordered_map<Edge*, std::unordered_map<Edge*, int>> votes;
	for (auto& conn : connections)
	{
		++votes[conn.first][conn.second];
	}

	for (auto edgeConnections : votes)
	{
		Edge* pairedEdge = nullptr;
		if (edgeConnections.second.size() == 1)
		{
			pairedEdge = edgeConnections.second.begin()->first;
		}
		else
		{
			int maxVotes = 0;
			Edge* maxEdge = nullptr;
			for (auto connCount : edgeConnections.second) 
			{
				if (maxVotes < connCount.second)
				{
					maxVotes = connCount.second;
					maxEdge = connCount.first;
				}
			}
			bool reliable = true;
			for (auto connCount : edgeConnections.second) 
			{
				if (connCount.first != maxEdge && 
					maxVotes / connCount.second < 2) reliable = false;
			}
			if (reliable) pairedEdge = maxEdge;
		}

		if (!pairedEdge)
		{
			Logger::get().debug() << "Ambigious "
					<< edgeConnections.first->seqId
				    << "\t" << edgeConnections.first->seqEnd;

			continue;
		}

		Logger::get().debug() << "Connection " 
				<< edgeConnections.first->seqId
				<< "\t" << edgeConnections.first->seqEnd << "\t"
				<< pairedEdge->seqId << "\t" << pairedEdge->seqBegin;

		//modify the graph
		Knot::Id oldKnot = edgeConnections.first->knotEnd;
		Knot::Id newKnot = _knots.size();
		_knots.push_back(Knot(newKnot));
		vecRemove(_knots[oldKnot].inEdges, edgeConnections.first);
		vecRemove(_knots[oldKnot].outEdges, pairedEdge);
		_knots[newKnot].inEdges.push_back(edgeConnections.first);
		_knots[newKnot].outEdges.push_back(pairedEdge);
		edgeConnections.first->knotEnd = newKnot;
		pairedEdge->knotBegin = newKnot;

		//std::cout << "Votes for " << edgeVotes.first->seqId
		//		  << " " << edgeVotes.first->seqEnd << std::endl;

		/*for (auto votesIter : edgeConnections.second)
		{
			std::cout << "\t" << votesIter.first->seqId << " "
					  << votesIter.first->seqBegin << " " 
					  << votesIter.second << std::endl;
		}*/
	}
}

std::vector<std::pair<Edge*, Edge*>>
	AssemblyGraph::getConnections(const OverlapContainer& ovlpContainer)
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

	std::vector<std::pair<Edge*, Edge*>> connections;
	for (auto seqOvelaps : ovlpContainer.getOverlapIndex())
	{
		std::unordered_map<Edge*, int32_t> inSupport;
		std::unordered_map<Edge*, int32_t> outSupport;
		for (auto& ovlp : seqOvelaps.second)
		{
			for (Edge* edge : inEdges)
			{
				if (ovlp.extId != edge->seqId) continue;
				
				int32_t intersect = std::min(ovlp.extEnd, edge->seqEnd) -
									std::max(ovlp.extBegin, edge->seqBegin);

				//TODO: overhangs wrt to unique edges
				
				if (intersect < 0 ||
					edge->seqEnd - ovlp.extEnd > OVERHANG ||
					edge->seqEnd - ovlp.extBegin < OVERHANG ||
					std::min(ovlp.curBegin, edge->seqBegin) > OVERHANG) 
				{
					continue;
				}

				int32_t readIntersection = ovlp.curBegin + edge->seqEnd - 
										   ovlp.extBegin;
				inSupport[edge] = readIntersection;
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

				int32_t readIntersection = ovlp.curEnd - 
										   (ovlp.extEnd - edge->seqBegin);
				outSupport[edge] = readIntersection;
			}
		}

		if (inSupport.size() == 1 && outSupport.size() == 1)
		{
			auto leftEdge = *inSupport.begin();
			auto rightEdge = *outSupport.begin();

			if (leftEdge.first->knotEnd == rightEdge.first->knotBegin &&
				rightEdge.second - leftEdge.second > Parameters::minimumOverlap)
			{
				connections.push_back({leftEdge.first, rightEdge.first});
				/*
				std::cout << _seqReads.seqName(seqOvelaps.first) << std::endl;
				std::cout << "Connection: " << leftEdge.first->seqId
						  << "\t" << leftEdge.first->seqEnd << "\t"
						  << rightEdge.first->seqId << "\t" 
						  << rightEdge.first->seqBegin << "\t" << leftEdge.second 
						  << "\t" << rightEdge.second << std::endl;*/
			}
		}
	}
	return connections;
}
