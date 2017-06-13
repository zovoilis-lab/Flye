//(c) 2016 by Authors
//This file is a part of ABruijn program.
//Released under the BSD license (see LICENSE file)

#include <cmath>


#include "repeat_resolver.h"
#include "../common/config.h"
#include "../common/utils.h"
#include "bipartie_mincost.h"


GraphAlignment
	RepeatResolver::chainReadAlignments(const SequenceContainer& edgeSeqs,
								 	    std::vector<EdgeAlignment> ovlps)
{
	std::sort(ovlps.begin(), ovlps.end(),
			  [](const EdgeAlignment& e1, const EdgeAlignment& e2)
			  	{return e1.overlap.curBegin < e2.overlap.curBegin;});

	typedef std::vector<EdgeAlignment*> Chain;

	std::list<Chain> activeChains;
	for (auto& edgeAlignment : ovlps)
	{
		std::list<Chain> newChains;
		int32_t maxSpan = 0;
		Chain* maxChain = nullptr;
		for (auto& chain : activeChains)
		{
			OverlapRange& nextOvlp = edgeAlignment.overlap;
			OverlapRange& prevOvlp = chain.back()->overlap;

			int32_t readDiff = nextOvlp.curBegin - prevOvlp.curEnd;
			int32_t graphDiff = nextOvlp.extBegin +
								prevOvlp.extLen - prevOvlp.extEnd;

			if (_readJump > readDiff && readDiff > 0 &&
				_readJump > graphDiff && graphDiff > 0 &&
				abs(readDiff - graphDiff) < _readJump / Constants::farJumpRate &&
				chain.back()->edge->nodeRight == edgeAlignment.edge->nodeLeft)
			{
				int32_t readSpan = nextOvlp.curEnd -
								   chain.front()->overlap.curBegin;
				if (readSpan > maxSpan)
				{
					maxSpan = readSpan;
					maxChain = &chain;
				}
			}
		}
		
		if (maxChain)
		{
			newChains.push_back(*maxChain);
			maxChain->push_back(&edgeAlignment);
		}

		activeChains.splice(activeChains.end(), newChains);
		activeChains.push_back({&edgeAlignment});
	}

	int32_t maxSpan = 0;
	Chain* maxChain = nullptr;
	for (auto& chain : activeChains)
	{
		int32_t readSpan = chain.back()->overlap.curEnd - 
						   chain.front()->overlap.curBegin;
		//if (readSpan > Parameters::get().minimumOverlap && readSpan > maxSpan)
		if (readSpan > maxSpan)
		{
			maxSpan = readSpan;
			maxChain = &chain;
		}
	}

	GraphAlignment result;
	if (maxChain)
	{
		//check length consistency
		int32_t readSpan = maxChain->back()->overlap.curEnd - 
						   maxChain->front()->overlap.curBegin;
		int32_t graphSpan = maxChain->front()->overlap.extRange();
		for (size_t i = 1; i < maxChain->size(); ++i)
		{
			graphSpan += (*maxChain)[i]->overlap.extEnd +
						 (*maxChain)[i - 1]->overlap.extLen - 
						 (*maxChain)[i - 1]->overlap.extEnd;	
		}
		float lengthDiff = abs(readSpan - graphSpan);
		float meanLength = (readSpan + graphSpan) / 2.0f;
		if (lengthDiff > meanLength / Constants::overlapDivergenceRate)
		{
			return {};
		}

		for (auto& aln : *maxChain) result.push_back(*aln);
	}

	return result;
}

void RepeatResolver::separatePath(const GraphPath& graphPath, 
								  SequenceSegment readSegment, size_t newId)
{
	//first edge
	GraphNode* leftNode = _graph.addNode();
	vecRemove(graphPath.front()->nodeRight->inEdges, graphPath.front());
	graphPath.front()->nodeRight = leftNode;
	leftNode->inEdges.push_back(graphPath.front());

	//repetitive edges in the middle
	for (size_t i = 1; i < graphPath.size() - 1; ++i)
	{
		--graphPath[i]->multiplicity;
	}

	GraphNode* rightNode = leftNode;
	if (graphPath.size() > 2)
	{
		rightNode = _graph.addNode();
		GraphEdge* newEdge = _graph.addEdge(GraphEdge(leftNode, rightNode,
													  FastaRecord::Id(newId)));
		newEdge->multiplicity = 1;
		newEdge->seqSegments.push_back(readSegment);
		//newEdge->readSequence = true;
	}

	//last edge
	vecRemove(graphPath.back()->nodeLeft->outEdges, graphPath.back());
	graphPath.back()->nodeLeft = rightNode;
	rightNode->outEdges.push_back(graphPath.back());
}

void RepeatResolver::resolveConnections(const std::vector<Connection>& connections)
{
	//Logger::get().debug() << "Resolving repeats";

	///////////
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
	}
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

	const float MIN_SUPPORT = 0.5f;
	std::unordered_set<FastaRecord::Id> usedEdges;
	std::vector<Connection> uniqueConnections;
	int totalLinks = 0;
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
			<< "\t" << support << "\t" << confidence;

		if (support < MIN_SUPPORT) continue;

		totalLinks += 2;
		for (auto& conn : connections)
		{
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

		this->separatePath(conn.path, conn.readSequence, _graph._nextEdgeId);
		this->separatePath(complPath, complSegment, _graph._nextEdgeId + 1);
		_graph._nextEdgeId += 2;
	}

	Logger::get().debug() << "Edges: " << totalLinks << " links: "
						  << connections.size();
}


void RepeatResolver::resolveRepeats()
{
	int PATHS_TO_SKIP = 5;
	std::unordered_map<GraphEdge*, int> coveredEdges;
	for (auto& readPath : _readAlignments)
	{
		std::vector<GraphEdge*> safeEdges;
		for (auto& aln : readPath)
		{
			if (!aln.edge->isRepetitive() &&
				aln.edge->length() > Constants::maximumJump) 
			{
				safeEdges.push_back(aln.edge);
			}
		}
		if (safeEdges.size() > 2)
		{
			for (size_t i = 1; i < safeEdges.size() - 1; ++i)
			{
				coveredEdges[safeEdges[i]] += 1;
			}
		}
	}
	std::unordered_set<GraphEdge*> skipEdges;
	for (auto& edge : coveredEdges)
	{
		if (edge.second >= PATHS_TO_SKIP)
		{
			skipEdges.insert(edge.first);
			Logger::get().debug() << "Skip: " << edge.first->edgeId.signedId() 
				<< " " << edge.second;
		}
	}

	auto connections = this->getConnections(skipEdges);
	this->resolveConnections(connections);

	//one more time
	connections = this->getConnections(skipEdges);
	this->resolveConnections(connections);

	this->clearResolvedRepeats();
}

std::vector<RepeatResolver::Connection> 
	RepeatResolver::getConnections(const std::unordered_set<GraphEdge*> skipEdges)
{
	auto safeEdge = [&skipEdges](GraphEdge* edge)
	{
		return !edge->isRepetitive() && 
			   edge->length() > Constants::maximumJump &&
			   !skipEdges.count(edge);
	};

	std::vector<Connection> readConnections;
	for (auto& readPath : _readAlignments)
	{
		/*if (readPath.size() > 1)
		{
			Logger::get().debug() 
				<< _readSeqs.seqName(readPath.front().overlap.curId);
		}*/

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
			
			/*if (readPath.size() > 1)
			{
				Logger::get().debug() << aln.edge->edgeId.signedId() << "\t" 
									  << aln.overlap.curBegin << "\t"
									  << aln.overlap.curEnd << "\t"
									  << aln.overlap.curRange();
			}*/

			currentPath.push_back(aln.edge);
			if (safeEdge(aln.edge) && currentPath.size() > 1)
			{
				if (!currentPath.back()->nodeLeft->isBifurcation() &&
					!currentPath.front()->nodeRight->isBifurcation()) continue;

				//check that the path is still viable
				//in case of 2nd RR iteration
				bool inconsistent = false;
				for (size_t i = 0; i < currentPath.size() - 1; ++i)
				{
					if (currentPath[i]->nodeRight != 
						currentPath[i + 1]->nodeLeft) inconsistent = true;
				}
				if (inconsistent) continue;

				GraphPath complPath = _graph.complementPath(currentPath);

				int32_t readEnd = aln.overlap.curBegin - aln.overlap.extBegin;
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

	std::unordered_set<GraphNode*> toRemove;

	for (auto& node : _graph.iterNodes())
	{
		//separated nodes
		if (node->neighbors().size() == 0)
		{
			bool resolved = true;
			for (auto& edge : node->outEdges) 
			{
				if (edge->multiplicity > 0 &&
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
			if (edge->multiplicity > 0) resolvedRepeat = false;
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



void RepeatResolver::alignReads()
{
	const int MAX_KMER_COUNT = 1000;

	//create database
	std::unordered_map<FastaRecord::Id, 
					   std::pair<GraphEdge*, SequenceSegment*>> idToSegment;
	SequenceContainer pathsContainer;

	for (auto& edge : _graph.iterEdges())
	{
		for (auto& segment : edge->seqSegments)
		{
			size_t len = segment.end - segment.start;
			FastaRecord::DnaRepr sequence = _asmSeqs.getSeq(segment.seqId)
												.substr(segment.start, len);
			auto& newRec = pathsContainer.addSequence(sequence, "");
			idToSegment[newRec.id] = {edge, &segment};
			//idToSegment[newRec.id.rc()] = {edge, &segment};
		}
	}

	//index it and align reads
	VertexIndex pathsIndex(pathsContainer);
	pathsIndex.countKmers(1);
	pathsIndex.buildIndex(1, MAX_KMER_COUNT, 1);
	//OverlapDetector readsOverlapper(pathsContainer, pathsIndex, _readJump,
	//								_maxSeparation - _readOverhang, 0);
	OverlapDetector readsOverlapper(pathsContainer, pathsIndex, _readJump,
									_maxSeparation, 0);
	OverlapContainer readsOverlaps(readsOverlapper, _readSeqs, false);
	readsOverlaps.findAllOverlaps();

	Logger::get().debug() << "Threading reads through the graph";
	//get connections
	int numAligned = 0;
	for (auto& readId : _readSeqs.getIndex())
	{
		auto& overlaps = readsOverlaps.getOverlapIndex().at(readId.first);
		std::vector<EdgeAlignment> alignments;

		for (auto& ovlp : overlaps)
		{
			if (idToSegment.count(ovlp.extId))
			{
				alignments.push_back({ovlp, idToSegment[ovlp.extId].first,
									  idToSegment[ovlp.extId].second});
			}
		}

		_readAlignments.push_back(this->chainReadAlignments(pathsContainer, 
															alignments));
		if (!_readAlignments.back().empty()) ++numAligned;
	}

	Logger::get().debug() << "Aligned " << numAligned << " / " 
		<< _readSeqs.getIndex().size();
}
