//(c) 2016 by Authors
//This file is a part of ABruijn program.
//Released under the BSD license (see LICENSE file)

#include <cmath>


#include "repeat_resolver.h"
#include "../common/config.h"
#include "../common/utils.h"
#include "../common/parallel.h"
#include "bipartie_mincost.h"


void RepeatResolver::separatePath(const GraphPath& graphPath, 
								  SequenceSegment readSegment, 
								  FastaRecord::Id newId)
{
	//first edge
	GraphNode* leftNode = _graph.addNode();
	vecRemove(graphPath.front()->nodeRight->inEdges, graphPath.front());
	graphPath.front()->nodeRight = leftNode;
	leftNode->inEdges.push_back(graphPath.front());
	int32_t pathCoverage = (graphPath.front()->meanCoverage +
						    graphPath.back()->meanCoverage) / 2;

	//repetitive edges in the middle
	for (size_t i = 1; i < graphPath.size() - 1; ++i)
	{
		graphPath[i]->resolved = true;
		graphPath[i]->meanCoverage = 
			std::max(graphPath[i]->meanCoverage - pathCoverage, 0);
	}

	GraphNode* rightNode = leftNode;
	if (graphPath.size() > 2)
	{
		rightNode = _graph.addNode();
		GraphEdge* newEdge = _graph.addEdge(GraphEdge(leftNode, rightNode,
													  newId));
		newEdge->seqSegments.push_back(readSegment);
		newEdge->meanCoverage = _multInf.getMeanCoverage();
	}

	//last edge
	vecRemove(graphPath.back()->nodeLeft->outEdges, graphPath.back());
	graphPath.back()->nodeLeft = rightNode;
	rightNode->outEdges.push_back(graphPath.back());
}

int RepeatResolver::resolveConnections(const std::vector<Connection>& connections)
{
	///////////
	/*
	std::unordered_map<GraphEdge*, std::unordered_map<GraphEdge*, int>> stats;
	for (auto& conn : connections)
	{
		++stats[conn.path.front()][conn.path.back()];
	}
	
	for (auto& leftEdge : stats)
	{
		Logger::get().debug() << "For " << leftEdge.first->edgeId.signedId() << " "
			<< leftEdge.first->seqSegments.front().seqId << " "
			<< leftEdge.first->seqSegments.front().end;
		
		for (auto& rightEdge : leftEdge.second)
		{
			Logger::get().debug() << "\t" << rightEdge.first->edgeId.signedId() << " "
				<< rightEdge.first->seqSegments.front().seqId << " "
				<< rightEdge.first->seqSegments.front().start << " " 
				<< rightEdge.second;
		}
		Logger::get().debug() << "";
	}*/
	///////////
	std::unordered_map<GraphEdge*, int> leftCoverage;
	std::unordered_map<GraphEdge*, int> rightCoverage;
	
	//create bipartie graph matrix
	std::unordered_map<GraphEdge*, size_t> leftEdgesId;
	std::unordered_map<size_t, GraphEdge*> leftIdToEdge;
	size_t nextLeftId = 0;
	std::unordered_map<GraphEdge*, size_t> rightEdgesId;
	std::unordered_map<size_t, GraphEdge*> rightIdToEdge;
	size_t nextRightId = 0;

	for (auto& conn : connections)
	{
		GraphEdge* leftEdge = conn.path.front();
		GraphEdge* rightEdge = conn.path.back();
		++leftCoverage[leftEdge];
		++rightCoverage[rightEdge];

		if (!leftEdgesId.count(leftEdge))
		{
			leftEdgesId[leftEdge] = nextLeftId;
			leftIdToEdge[nextLeftId++] = leftEdge;
		}
		if (!rightEdgesId.count(rightEdge))
		{
			rightEdgesId[rightEdge] = nextRightId;
			rightIdToEdge[nextRightId++] = rightEdge;
		}
	}

	size_t numNodes = std::max(leftEdgesId.size(), rightEdgesId.size());
	BipartieTable table;
	table.assign(numNodes, std::vector<double>(numNodes, 0));
	for (auto& conn : connections)
	{
		GraphEdge* leftEdge = conn.path.front();
		GraphEdge* rightEdge = conn.path.back();
		if (leftEdge->edgeId == rightEdge->edgeId ||
			leftEdge->edgeId == rightEdge->edgeId.rc()) continue;

		//solving min cost mathcing
		--table[leftEdgesId[leftEdge]][rightEdgesId[rightEdge]];
	}
	auto edges = bipartieMincost(table);
	typedef std::pair<size_t, size_t> MatchPair;
	std::vector<MatchPair> matchingPairs;
	for (size_t i = 0; i < edges.size(); ++i)
	{
		matchingPairs.emplace_back(i, edges[i]);
	}

	std::unordered_set<FastaRecord::Id> usedEdges;
	std::vector<Connection> uniqueConnections;
	int totalLinks = 0;
	int unresolvedLinks = 0;
	for (auto match : matchingPairs)
	{
		GraphEdge* leftEdge = leftIdToEdge[match.first];
		GraphEdge* rightEdge = rightIdToEdge[match.second];

		int support = -table[match.first][match.second];
		float confidence = 2.0f * support / (leftCoverage[leftEdge] + 
											 rightCoverage[rightEdge]);
		if (!support) continue;
		if (usedEdges.count(leftEdge->edgeId)) continue;
		usedEdges.insert(rightEdge->edgeId.rc());

		Logger::get().debug() << "\tConnection " 
			<< leftEdge->edgeId.signedId()
			<< "\t" << rightEdge->edgeId.signedId()
			<< "\t" << support / 2 << "\t" << confidence;

		if (confidence < Constants::minRepeatResSupport) 
		{
			++unresolvedLinks;
			continue;
		}
		//if (support < 4) continue;

		totalLinks += 2;
		for (auto& conn : connections)
		{
			//TODO: choose representetive read more carefully
			if (conn.path.front() == leftEdge && 
				conn.path.back() == rightEdge)
			{
				uniqueConnections.push_back(conn);
				break;
			}
		}
	}

	for (auto& conn : uniqueConnections)
	{
		GraphPath complPath = _graph.complementPath(conn.path);
		SequenceSegment complSegment = conn.readSequence.complement();

		FastaRecord::Id edgeId = _graph.newEdgeId();
		this->separatePath(conn.path, conn.readSequence, edgeId);
		this->separatePath(complPath, complSegment, edgeId.rc());
	}

	Logger::get().debug() << "Resolved: " << totalLinks / 2 << " links: "
						  << connections.size() / 2;
	Logger::get().debug() << "Unresolved: " << unresolvedLinks / 2;

	return totalLinks / 2;
}

//TODO: operate on unbranching paths rather than single edges
void RepeatResolver::findRepeats()
{
	std::unordered_map<GraphEdge*, 
					   std::unordered_map<GraphEdge*, int>> outConnections;

	//all good at the beginning
	for (auto& edge : _graph.iterEdges())
	{
		edge->repetitive = false;
	}

	//mark edges with high coverage as repetitive
	for (auto& edge : _graph.iterEdges())
	{
		if (!edge->edgeId.strand()) continue;

		if (edge->meanCoverage > _multInf.getUniqueCovThreshold() * 2)
		{
			edge->repetitive = true;
			_graph.complementEdge(edge)->repetitive = true;
		}
		//plus tanem repeat
		if (edge->isLooped())
		{
			std::unordered_set<FastaRecord::Id> seen;
			for (auto& seg : edge->seqSegments)
			{
				if (seen.count(seg.seqId))
				{
					edge->repetitive = true;
					_graph.complementEdge(edge)->repetitive = true;
				}
				seen.insert(seg.seqId);
			}
		}
	}

	//Now, using read alignments
	//extract read alignments
	for (auto& readPath : _aligner.getAlignments())
	{
		if (readPath.size() < 2) continue;
		int overhang = std::max(readPath.front().overlap.curBegin,
								readPath.back().overlap.curLen - 
									readPath.back().overlap.curEnd);
		if (overhang > Constants::maximumOverhang) continue;

		for (size_t i = 0; i < readPath.size() - 1; ++i)
		{
			if (readPath[i].edge == readPath[i + 1].edge &&
				readPath[i].edge->isLooped()) continue;
			if (readPath[i].edge->edgeId == 
				readPath[i + 1].edge->edgeId.rc()) continue;

			GraphEdge* complLeft = _graph.complementEdge(readPath[i].edge);
			GraphEdge* complRight = _graph.complementEdge(readPath[i + 1].edge);
			++outConnections[readPath[i].edge][readPath[i + 1].edge];
			++outConnections[complRight][complLeft];
		}
	}
	//summarizes multiplicity
	auto edgeMultiplicity = [this, &outConnections] (GraphEdge* edge)
	{
		int maxSupport = 0;
		for (auto& outConn : outConnections[edge])
		{
			if (maxSupport < outConn.second)
			{
				maxSupport =  outConn.second;
			}
		}

		int repeatMult = 0;
		int uniqueMult = 0;
		int minSupport = maxSupport / Constants::outPathsRatio;
		int minCoverage = _multInf.getMeanCoverage() / Constants::readCovRate;
		for (auto& outConn : outConnections[edge]) 
		{
			if (outConn.second > minSupport && 
				outConn.first->meanCoverage > minCoverage)
			{
				outConn.first->repetitive ? ++repeatMult : ++uniqueMult;
			}
		}
		return std::min(repeatMult, 1) + uniqueMult;
	};

	//order might be important, process short edges first
	std::vector<GraphEdge*> sortedEdges;
	for (auto& edge : _graph.iterEdges()) sortedEdges.push_back(edge);
	std::sort(sortedEdges.begin(), sortedEdges.end(),
			  [](const GraphEdge* e1, const GraphEdge* e2) 
			  {return e1->length() < e2->length();});
	for (auto& edge : sortedEdges)
	{
		if (!edge->edgeId.strand()) continue;

		GraphEdge* complEdge = _graph.complementEdge(edge);
		int rightMult = edgeMultiplicity(edge);
		int leftMult = edgeMultiplicity(complEdge);
		int mult = std::max(leftMult, rightMult);
		if (mult > 1) 
		{
			edge->repetitive = true;
			complEdge->repetitive = true;
		}
	}
	
	//propagate within unbranching paths, jumping over loops
	std::unordered_set<GraphEdge*> visited;
	for (GraphEdge* edge : _graph.iterEdges())
	{
		if (!edge->isRepetitive() || edge->isLooped() ||
			visited.count(edge)) continue;

		GraphEdge* curEdge = edge;
		visited.insert(edge);
		while(true)
		{
			if (curEdge->nodeRight->inEdges.size() != 1 ||
				curEdge->nodeRight->outEdges.size() != 1) break;

			GraphEdge* nextEdge = curEdge->nodeRight->outEdges.front();
			if (visited.count(nextEdge)) break;
			visited.insert(curEdge);
			curEdge = nextEdge;
			curEdge->repetitive = true;
			_graph.complementEdge(curEdge)->repetitive = true;
		}
	}

	//////////////////
	///////some logging
	for (auto& edge : _graph.iterEdges())
	{
		if (!edge->edgeId.strand()) continue;
		GraphEdge* complEdge = _graph.complementEdge(edge);
		int rightMult = edgeMultiplicity(edge);
		int leftMult = edgeMultiplicity(complEdge);
		int mult = std::max(leftMult, rightMult);
		bool match = (edge->multiplicity > 1) == (edge->repetitive);
		if (!match && edge->multiplicity != 0 && !edge->resolved)
		{
			std::string star = edge->repetitive ? "R" : " ";
			std::string loop = edge->isLooped() ? "L" : " ";
			Logger::get().debug() << star << " " << loop << " " 
				<< edge->edgeId.signedId()
				<< "\t" << edge->multiplicity << " -> " << mult << " ("
				<< leftMult << "," << rightMult << ") " << edge->length() << "\t"
				<< edge->meanCoverage;
		}
	}
	for (auto& edgeList : outConnections)
	{
		bool match = (edgeList.first->multiplicity > 1) == 
					  (edgeList.first->repetitive);
		if (!match && edgeList.first->multiplicity != 0 && 
			!edgeList.first->resolved)
		{
			Logger::get().debug() << "Outputs: " << edgeList.first->edgeId.signedId()
				<< " " << edgeList.first->multiplicity;
			for (auto& outEdgeCount : edgeList.second)
			{
				std::string star = outEdgeCount.first->repetitive ? "R" : " ";
				std::string loop = outEdgeCount.first->isLooped() ? "L" : " ";
				Logger::get().debug() << star << " " << loop << " " 
					<< outEdgeCount.first->edgeId.signedId() << "\t" << outEdgeCount.second;
			}
			Logger::get().debug() << "";
		}
	}
}

void RepeatResolver::resolveRepeats()
{
	while (true)
	{
		auto connections = this->getConnections();
		int resolvedConnections = this->resolveConnections(connections);
		if (!resolvedConnections) break;

		_aligner.updateAlignments();
		this->findRepeats();
	}

	this->removeUnsupportedEdges();
	this->clearResolvedRepeats();
}


std::vector<RepeatResolver::Connection> 
	RepeatResolver::getConnections()
{
	
	auto safeEdge = [this](GraphEdge* edge)
	{
		return !edge->isRepetitive() &&
			edge->meanCoverage > _multInf.getMeanCoverage() / 
									Constants::readCovRate;
	};

	int totalSafe = 0;
	for (GraphEdge* edge : _graph.iterEdges())
	{
		if (edge->edgeId.strand() && safeEdge(edge)) ++totalSafe;
	}
	Logger::get().debug() << "Total unique edges: " << totalSafe;

	std::vector<Connection> readConnections;
	for (auto& readPath : _aligner.getAlignments())
	{
		GraphPath currentPath;
		int32_t readStart = 0;
		for (auto& aln : readPath)
		{
			if (currentPath.empty()) 
			{
				if (!safeEdge(aln.edge)) continue;
				readStart = aln.overlap.curEnd + aln.overlap.extLen - 
							aln.overlap.extEnd;
			}

			currentPath.push_back(aln.edge);
			if (safeEdge(aln.edge) && currentPath.size() > 1)
			{
				if (!currentPath.back()->nodeLeft->isBifurcation() &&
					!currentPath.front()->nodeRight->isBifurcation()) continue;

				GraphPath complPath = _graph.complementPath(currentPath);

				int32_t readEnd = aln.overlap.curBegin - aln.overlap.extBegin;
				readEnd = std::max(readStart + 100, readEnd);
				SequenceSegment segment(aln.overlap.curId, aln.overlap.curLen, 
										readStart, readEnd);
				segment.readSequence = true;
				SequenceSegment complSegment = segment.complement();

				readConnections.push_back({currentPath, segment});
				readConnections.push_back({complPath, complSegment});

				currentPath.clear();
				currentPath.push_back(aln.edge);
				readStart = aln.overlap.curEnd + aln.overlap.extLen - 
							aln.overlap.extEnd;
			}
		}
	}

	return readConnections;
}

void RepeatResolver::removeUnsupportedEdges()
{
	int coverageThreshold = _multInf.getMeanCoverage() / Constants::readCovRate;
	Logger::get().debug() << "Read coverage cutoff: " << coverageThreshold;

	std::unordered_set<GraphEdge*> edgesRemove;
	for (auto& edge : _graph.iterEdges())
	{
		GraphEdge* complEdge = _graph.complementEdge(edge);
		if (edge->meanCoverage <= coverageThreshold &&
			!edge->resolved)
		{
			edgesRemove.insert(edge);
			edgesRemove.insert(complEdge);
		}
	}
	for (auto& edge : edgesRemove) _graph.removeEdge(edge);
	Logger::get().debug() << "Removed " << edgesRemove.size() 
		<< " unsupported edges";
}

void RepeatResolver::clearResolvedRepeats()
{
	const int MIN_LOOP = Parameters::get().minimumOverlap;
	auto nextEdge = [](GraphNode* node)
	{
		for (auto edge : node->outEdges)
		{
			if (!edge->isLooped()) return edge;
		}
		return (GraphEdge*)nullptr;
	};

	auto shouldRemove = [](GraphEdge* edge)
	{
		//return edge->isRepetitive() && edge->resolved;
		return edge->resolved;
	};

	std::unordered_set<GraphNode*> toRemove;

	for (auto& node : _graph.iterNodes())
	{
		//separated nodes
		if (node->neighbors().size() == 0)
		{
			bool resolved = true;
			for (auto& edge : node->outEdges) 
			{
				if (!shouldRemove(edge) &&
					edge->length() > MIN_LOOP) resolved = false;
			}

			if (resolved) toRemove.insert(node);
		}

		//other nodes
		if (!node->isEnd()) continue;

		GraphEdge* direction = nextEdge(node);
		if (!direction) continue;

		GraphPath traversed;
		traversed.push_back(direction);
		GraphNode* curNode = direction->nodeRight;
		while (curNode->isResolved())
		{
			traversed.push_back(nextEdge(curNode));
			curNode = traversed.back()->nodeRight;
		}
		if (traversed.empty()) continue;

		bool removeLast = curNode->isEnd();
		bool resolvedRepeat = true;
		for (auto& edge : traversed) 
		{
			if (!shouldRemove(edge)) resolvedRepeat = false;
		}

		GraphPath complPath = _graph.complementPath(traversed);
		if (resolvedRepeat)
		{
			//first-last
			toRemove.insert(traversed.front()->nodeLeft);
			if (removeLast) toRemove.insert(complPath.front()->nodeLeft);

			//middle nodes
			for (size_t i = 0; i < traversed.size() - 1; ++i)
			{
				toRemove.insert(traversed[i]->nodeRight);
				toRemove.insert(complPath[i]->nodeRight);
			}

			//last-first
			if (removeLast) toRemove.insert(traversed.back()->nodeRight);
			toRemove.insert(complPath.back()->nodeRight);
		}
	}

	for (auto node : toRemove) _graph.removeNode(node);
}
