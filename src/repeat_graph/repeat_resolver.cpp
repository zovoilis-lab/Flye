//(c) 2016 by Authors
//This file is a part of ABruijn program.
//Released under the BSD license (see LICENSE file)

#include <cmath>


#include "repeat_resolver.h"
#include "../common/config.h"
#include "../common/utils.h"
#include "bipartie_mincost.h"


namespace
{
	struct Chain
	{
		Chain(): score(0) {}
		std::vector<EdgeAlignment*> aln;
		int32_t score;
	};
}

GraphAlignment
	RepeatResolver::chainReadAlignments(const SequenceContainer& edgeSeqs,
								 	    std::vector<EdgeAlignment> ovlps)
{
	std::sort(ovlps.begin(), ovlps.end(),
			  [](const EdgeAlignment& e1, const EdgeAlignment& e2)
			  	{return e1.overlap.curBegin < e2.overlap.curBegin;});

	std::list<Chain> activeChains;
	for (auto& edgeAlignment : ovlps)
	{
		std::list<Chain> newChains;
		int32_t maxScore = 0;
		Chain* maxChain = nullptr;
		for (auto& chain : activeChains)
		{
			OverlapRange& nextOvlp = edgeAlignment.overlap;
			OverlapRange& prevOvlp = chain.aln.back()->overlap;

			int32_t readDiff = nextOvlp.curBegin - prevOvlp.curEnd;
			int32_t graphDiff = nextOvlp.extBegin +
								prevOvlp.extLen - prevOvlp.extEnd;
			int32_t maxDiscordance = std::max(_readJump / Constants::farJumpRate,
											  _maxSeparation);

			if (_readJump > readDiff && readDiff > -500 &&
				_readJump > graphDiff && graphDiff > -500 &&
				abs(readDiff - graphDiff) < maxDiscordance &&
				chain.aln.back()->edge->nodeRight == edgeAlignment.edge->nodeLeft)
			{
				int32_t score = chain.score + nextOvlp.score;
				if (score > maxScore)
				{
					maxScore = score;
					maxChain = &chain;
				}
			}
		}
		
		if (maxChain)
		{
			newChains.push_back(*maxChain);
			maxChain->aln.push_back(&edgeAlignment);
			maxChain->score = maxScore;
		}

		activeChains.splice(activeChains.end(), newChains);
		activeChains.push_back(Chain());
		activeChains.back().aln.push_back(&edgeAlignment);
		activeChains.back().score = edgeAlignment.overlap.score;
	}

	int32_t maxScore = 0;
	Chain* maxChain = nullptr;
	for (auto& chain : activeChains)
	{
		if (chain.score > maxScore)
		{
			maxScore = chain.score;
			maxChain = &chain;
		}
	}

	GraphAlignment result;
	if (maxChain)
	{
		//check length consistency
		/*int32_t readSpan = maxChain->aln.back()->overlap.curEnd - 
						   maxChain->aln.front()->overlap.curBegin;
		int32_t graphSpan = maxChain->aln.front()->overlap.extRange();
		for (size_t i = 1; i < maxChain->aln.size(); ++i)
		{
			graphSpan += maxChain->aln[i]->overlap.extEnd +
						 maxChain->aln[i - 1]->overlap.extLen - 
						 maxChain->aln[i - 1]->overlap.extEnd;	
		}
		float lengthDiff = abs(readSpan - graphSpan);
		float meanLength = (readSpan + graphSpan) / 2.0f;
		if (lengthDiff > meanLength / Constants::overlapDivergenceRate)
		{
			return {};
		}*/

		for (auto& aln : maxChain->aln) result.push_back(*aln);
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
			<< "\t" << support / 2 << "\t" << confidence;

		if (support < Constants::minRepeatResSupport) continue;

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

	Logger::get().debug() << "Resolved: " << totalLinks / 2 << " links: "
						  << connections.size() / 2;
}


void RepeatResolver::findRepeats(int uniqueCovThreshold)
{
	std::unordered_map<GraphEdge*, 
					   std::unordered_map<GraphEdge*, int>> outConnections;

	//first, by coverage
	for (auto& edge : _graph.iterEdges())
	{
		auto complEdge = _graph.complementPath({edge}).front();
				if (edge->meanCoverage > uniqueCovThreshold * 2) 
		{
			edge->repetitive = true;
			complEdge->repetitive = true;
		}
	}

	for (auto& readPath : _readAlignments)
	{
		if (readPath.size() < 2) continue;

		if (readPath.front().overlap.curBegin > 
			Constants::maximumOverhang) continue;
		if (readPath.back().overlap.curLen - readPath.back().overlap.curEnd > 
			Constants::maximumOverhang) continue;

		for (size_t i = 0; i < readPath.size() - 1; ++i)
		{
			if (readPath[i].edge == readPath[i + 1].edge &&
				readPath[i].edge->isLooped()) continue;

			bool goodLeft = readPath[i].overlap.extRange() > 
								_readJump ||
							(0 < i && i < readPath.size());
			bool goodRight = readPath[i + 1].overlap.extRange() > 
								_readJump ||
							(0 < i + 1 && i + 1 < readPath.size());
			if (goodLeft && goodRight)
			{
				++outConnections[readPath[i].edge][readPath[i + 1].edge];
			}
		}
		/*std::stringstream pathStr;
		for (auto& aln : readPath)
		{
			pathStr << "(" << aln.edge->edgeId.signedId() << " " <<
				aln.overlap.extLen << " " << aln.overlap.extBegin 
				<< " " << aln.overlap.extEnd << " " << aln.overlap.curBegin
				<< " " << aln.overlap.curEnd << ") -> ";
		}
		Logger::get().debug() << _readSeqs.seqName(readPath.front().overlap.curId);
		Logger::get().debug() << "Path: " << pathStr.str();*/
	}
	
	/*for (auto& edgeList : outConnections)
	{
		Logger::get().debug() << "Outputs: " << edgeList.first->edgeId.signedId()
			<< " " << edgeList.first->multiplicity;
		for (auto& outEdgeCount : edgeList.second)
		{
			Logger::get().debug() << "\t" << outEdgeCount.first->edgeId.signedId()
				<< " " << outEdgeCount.second << " " << outEdgeCount.first->isLooped();
		}
		Logger::get().debug() << "";
	}*/

	const int RATIO_THLD = 10;
	for (auto& edge : outConnections)
	{
		if (!edge.first->edgeId.strand()) continue;

		int rightCoverage = 0;
		int rightMult = 0;
		for (auto& outConn : edge.second)
		{
			rightCoverage = std::max(rightCoverage, outConn.second);
		}
		for (auto& outConn : edge.second) 
		{
			if (outConn.second > rightCoverage / RATIO_THLD) ++rightMult;
		}

		auto complEdge = _graph.complementPath({edge.first}).front();
		int leftCoverage = 0;
		int leftMult = 0;
		if (outConnections.count(complEdge))
		{
			for (auto& outConn : outConnections.at(complEdge))
			{
				leftCoverage = std::max(leftCoverage, outConn.second);
			}
			for (auto& outConn : outConnections.at(complEdge)) 
			{
				if (outConn.second > leftCoverage / RATIO_THLD) ++leftMult;
			}
		}

		int mult = std::max(leftMult, rightMult);
		if (mult > 1) 
		{
			edge.first->repetitive = true;
			complEdge->repetitive = true;
		}

		////////
		std::string match = (edge.first->multiplicity != 1) != 
							(edge.first->repetitive) ? "*" : " ";

		Logger::get().debug() << match << " " << edge.first->edgeId.signedId()
			<< " " << edge.first->multiplicity << " -> " << mult << " ("
			<< leftMult << "," << rightMult << ") " << edge.first->length() << "\t"
			<< edge.first->meanCoverage;
	}
}


void RepeatResolver::resolveRepeats()
{
	std::unordered_set<GraphEdge*> skipEdges;
	/*
	int PATHS_TO_SKIP = 5;
	std::unordered_map<GraphEdge*, int> coveredEdges;
	for (auto& readPath : _readAlignments)
	{
		std::vector<GraphEdge*> safeEdges;
		for (auto& aln : readPath)
		{
			if (!aln.edge->isRepetitive() && aln.edge->length() > _readJump) 
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
	for (auto& edge : coveredEdges)
	{
		if (edge.second >= PATHS_TO_SKIP)
		{
			skipEdges.insert(edge.first);
			Logger::get().debug() << "Skip: " << edge.first->edgeId.signedId() 
				<< " " << edge.second;
		}
	}*/

	auto connections = this->getConnections(skipEdges);
	this->resolveConnections(connections);

	//one more time
	//connections = this->getConnections(skipEdges);
	//this->resolveConnections(connections);

	this->clearResolvedRepeats();
}

std::vector<RepeatResolver::Connection> 
	RepeatResolver::getConnections(const std::unordered_set<GraphEdge*> skipEdges)
{
	
	auto safeEdge = [&skipEdges, this](GraphEdge* edge)
	{
		return !edge->isRepetitive() && 
			   edge->length() > _readJump &&
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
	//std::ofstream alnDump("../alignment_dump.txt");

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
	//								_maxSeparation - 2 * _readOverhang, 0);
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
		if (!_readAlignments.back().empty())
		{
			++numAligned;
			/*alnDump << _readSeqs.seqName(readId.first)
				<< "\t" << _readSeqs.seqLen(readId.first) << std::endl;
			for (auto& aln : _readAlignments.back())
			{
				alnDump << "\t" << aln.edge->edgeId.signedId()
					<< "\t" << aln.overlap.extBegin << "\t" << aln.overlap.extEnd
					<< "\t" << aln.overlap.extRange() << "\t|\t"
					<< aln.overlap.curBegin << "\t" << aln.overlap.curEnd
					<< "\t" << aln.overlap.curRange() << std::endl;
			}*/
		}
	}

	Logger::get().debug() << "Aligned " << numAligned << " / " 
		<< _readSeqs.getIndex().size();
}
