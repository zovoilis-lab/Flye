
#include "assembly_graph.h"
#include "config.h"
#include "logger.h"
#include "disjoint_set.h"
#include "overlap.h"

#include <fstream>
#include <map>

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

			_asmOverlaps.push_back(ovlp);
			_asmOverlaps.push_back(ovlp.reverse());
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
	}

	this->addEdgesPointers();
}


void AssemblyGraph::addEdgesPointers()
{
	for (auto& seqEdgesPair : _edges)
	{
		//prev and next
		auto prevEdge = seqEdgesPair.second.begin();
		auto nextEdge = prevEdge;
		++nextEdge;
		while (nextEdge != seqEdgesPair.second.end())
		{
			prevEdge->nextEdge = &(*nextEdge);
			nextEdge->prevEdge = &(*prevEdge);
			++prevEdge;
			++nextEdge;
		}

		//complementary edges
		if (!seqEdgesPair.first.strand()) continue;

		auto& complEdges = _edges[seqEdgesPair.first.rc()];
		if (seqEdgesPair.second.size() != complEdges.size())
		{
			throw std::runtime_error("Number of edges not equal for "
									 "complementary sequences");
		}

		auto fwdIt = seqEdgesPair.second.begin();
		auto cmpIt = complEdges.rbegin();
		while (fwdIt != seqEdgesPair.second.end())
		{
			fwdIt->complementEdge = &(*cmpIt);
			cmpIt->complementEdge = &(*fwdIt);
			++fwdIt;
			++cmpIt;
		}
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

/*
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
}*/

void AssemblyGraph::generatePathCandidates()
{
	const int32_t FLANK = 10000;

	for (auto& seqEdges : _edges)
	{
		for (auto& edge : seqEdges.second)
		{
			if (!edge.nextEdge) continue;

			//add database entry
			auto& asmSeq = _seqAssembly.getSeq(edge.seqId);
			int32_t leftFlank = std::max(edge.seqEnd - FLANK, 0);
			int32_t rightFlank = std::min(edge.nextEdge->seqBegin + FLANK, 
										  (int32_t)asmSeq.size());
			std::string sequence = asmSeq.substr(leftFlank, 
												 rightFlank - leftFlank);
			int32_t repeatStart = std::min(FLANK, edge.seqEnd);
			int32_t repeatEnd = repeatStart + edge.nextEdge->seqBegin - 
								edge.seqEnd;
			_pathCandidates.emplace_back(sequence, repeatStart, repeatEnd,
									     &edge, edge.nextEdge);
			//
			
			for (auto& outEdge : _knots[edge.knotEnd].outEdges)
			{
				if (outEdge == edge.nextEdge) continue;

				//find corresponding overlap
				OverlapRange maxOverlap;
				bool found = false;
				for (auto& ovlp : _asmOverlaps)
				{
					if (ovlp.curId == edge.seqId &&
						ovlp.curBegin >= edge.seqEnd &&
						ovlp.curEnd <= edge.nextEdge->seqBegin &&
						ovlp.extId == outEdge->seqId &&
						ovlp.extBegin >= outEdge->prevEdge->seqEnd &&
						ovlp.extEnd <= outEdge->seqBegin)
					{
						found = true;
						if (maxOverlap.curRange() < ovlp.curRange())
						{
							maxOverlap = ovlp;
						}
					}
				}
				if (!found) continue;
				//
				
				//create database entry
				auto& asmLeft = _seqAssembly.getSeq(edge.seqId);
				auto& asmRight = _seqAssembly.getSeq(outEdge->seqId);
				int32_t leftFlank = std::max(edge.seqEnd - FLANK, 0);
				int32_t rightFlank = std::min(outEdge->seqBegin + FLANK, 
											  (int32_t)asmRight.size());
				std::string sequence = 
					asmLeft.substr(leftFlank, maxOverlap.curEnd - leftFlank) +
					asmRight.substr(maxOverlap.extEnd, 
									rightFlank - maxOverlap.extEnd);
				int32_t repeatStart = std::min(FLANK, edge.seqEnd);
				int32_t repeatEnd = repeatStart + maxOverlap.curRange() +
									outEdge->seqBegin - maxOverlap.extEnd;
				_pathCandidates.emplace_back(sequence, repeatStart, repeatEnd,
											 &edge, outEdge);
				//
			}
		}
	}
}

void AssemblyGraph::outputFasta(const std::string& filename)
{
	std::ofstream fout(filename);
	size_t nextNewNode = _knots.size();
	for (auto& seqEdgesPair : _edges)
	{
		for (auto& edge : seqEdgesPair.second)
		{
			Knot::Id firstNode = -edge.knotBegin;
			if (firstNode == SEQ_BEGIN) firstNode = nextNewNode++;
			Knot::Id lastNode = edge.knotEnd;
			if (lastNode == SEQ_END) lastNode = nextNewNode++;

			fout << ">" << _seqAssembly.seqName(edge.seqId)
				 << "[" << edge.seqBegin << ":" << edge.seqEnd << "] "
				 << firstNode << " " << lastNode << std::endl;
			fout << _seqAssembly.getSeq(edge.seqId)
						.substr(edge.seqBegin, edge.seqEnd - edge.seqBegin) 
				 << std::endl;
		}
	}
	for (Knot& knot : _knots)
	{
		if (knot.knotId > SEQ_END && 
			(!knot.inEdges.empty() || !knot.outEdges.empty()))
		{
			fout << ">repeat " << knot.knotId << " " 
				 << -knot.knotId << std::endl;

			//a hack for complete genomes
			bool found = false;
			for (Edge* e1 : knot.inEdges)
			{
				for (Edge* e2 : knot.outEdges)
				{
					if (e1->seqId == e2->seqId)
					{
						found = true;

						fout << _seqAssembly.getSeq(e1->seqId)
								.substr(e1->seqEnd, e2->seqBegin - e1->seqEnd) 
				 			 << std::endl;
						break;
					}
				}
				if (found) break;
			}
		}
	}
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

std::vector<Connection> AssemblyGraph::getConnections()
{
	std::vector<Connection> connections;
	std::unordered_map<FastaRecord::Id, PathCandidate*> idToPath;

	SequenceContainer pathsContainer;
	for (size_t i = 0; i < _pathCandidates.size(); ++i)
	{
		pathsContainer.addSequence(FastaRecord(_pathCandidates[i].sequence, "", 
											   FastaRecord::Id(i)));
		idToPath[FastaRecord::Id(i)] = &_pathCandidates[i];
	}
	VertexIndex pathsIndex(pathsContainer);
	pathsIndex.countKmers(1);
	pathsIndex.buildIndex(1, 500);

	OverlapDetector readsOverlapper(pathsContainer, pathsIndex, 
									500, 
									Parameters::minimumOverlap,
									Constants::maximumOverhang);
	OverlapContainer readsContainer(readsOverlapper, _seqReads);
	readsContainer.findAllOverlaps();
	
	for (auto& seqOvelaps : readsContainer.getOverlapIndex())
	{
		//auto overlaps = filterOvlp(seqOvelaps.second);
		std::unordered_set<FastaRecord::Id> supported;
		OverlapRange maxOverlap;
		for (auto& ovlp : seqOvelaps.second)
		{
			if (ovlp.curRange() > maxOverlap.curRange())
			{
				maxOverlap = ovlp;
			}
		}

		PathCandidate* path = idToPath[maxOverlap.extId];
		if (path->repeatStart - maxOverlap.extBegin > 1000 &&
			maxOverlap.extEnd - path->repeatEnd > 1000)
		{
			std::cout << seqOvelaps.first << "\t";
			std::cout << path->inEdge->seqEnd << "\t" 
					  << path->outEdge->seqBegin << "\t"
					  << maxOverlap.curBegin << "\t" << maxOverlap.curEnd << "\t"
					  << maxOverlap.extBegin << "\t" << maxOverlap.extEnd
					  << std::endl;
			supported.insert(maxOverlap.extId);
		}

		if (supported.size() == 1)
		{
			PathCandidate* path = idToPath[*supported.begin()];
			connections.push_back({path->inEdge, path->outEdge});
			connections.push_back({path->outEdge->complementEdge,
								   path->inEdge->complementEdge});
		}
	}

	return connections;
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

void AssemblyGraph::untangle()
{
	this->generatePathCandidates();
	auto connections = this->getConnections();

	typedef std::pair<Edge*, Edge*> EdgePair;
	std::unordered_map<Knot::Id, 
					   std::unordered_map<EdgePair, int, pairhash>> edges;

	for (auto& conn : connections)
	{
		Knot::Id knot = conn.inEdge->knotEnd;
		++edges[knot][EdgePair(conn.inEdge, conn.outEdge)];
	}

	for (auto& knotEdges : edges)
	{
		Logger::get().debug() << "Knot" << knotEdges.first;
		std::vector<EdgePair> sortedEdges;
		for (auto& edgeCount : knotEdges.second)
		{
			sortedEdges.push_back(edgeCount.first);
		}
		std::sort(sortedEdges.begin(), sortedEdges.end(),
				  [&knotEdges](const EdgePair& e1, const EdgePair& e2)
				  {return knotEdges.second[e1] > knotEdges.second[e2];});

		std::unordered_set<Edge*> usedInEdges;
		std::unordered_set<Edge*> usedOutEdges;

		for (auto& edgePair : sortedEdges)
		{
			if (!usedInEdges.count(edgePair.first) &&
				!usedOutEdges.count(edgePair.second))
			{

				Logger::get().debug() << "\tConnection " 
					<< edgePair.first->seqId
					<< "\t" << edgePair.first->seqEnd << "\t"
					<< edgePair.second->seqId << "\t" << edgePair.second->seqBegin
					<< "\t" << knotEdges.second[edgePair];

				usedInEdges.insert(edgePair.first);
				usedOutEdges.insert(edgePair.second);

				//modify the graph
				Knot::Id oldKnot = edgePair.first->knotEnd;
				Knot::Id newKnot = _knots.size();
				_knots.push_back(Knot(newKnot));
				vecRemove(_knots[oldKnot].inEdges, edgePair.first);
				vecRemove(_knots[oldKnot].outEdges, edgePair.second);
				_knots[newKnot].inEdges.push_back(edgePair.first);
				_knots[newKnot].outEdges.push_back(edgePair.second);
				edgePair.first->knotEnd = newKnot;
				edgePair.second->knotBegin = newKnot;
				//
			}
		}
	}
}
