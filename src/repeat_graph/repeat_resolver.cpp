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

	//repetitive edges in the middle
	int64_t covSum = 0;
	for (size_t i = 1; i < graphPath.size() - 1; ++i)
	{
		//--graphPath[i]->multiplicity;
		graphPath[i]->resolved = true;
		covSum += graphPath[i]->meanCoverage;
	}
	int32_t numEdges = graphPath.size() - 2;
	int32_t meanCov = numEdges ? covSum / numEdges : 0;

	GraphNode* rightNode = leftNode;
	if (graphPath.size() > 2)
	{
		rightNode = _graph.addNode();
		GraphEdge* newEdge = _graph.addEdge(GraphEdge(leftNode, rightNode,
													  newId));
		newEdge->seqSegments.push_back(readSegment);
		newEdge->meanCoverage = meanCov;
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


void RepeatResolver::findRepeats()
{
	std::unordered_map<GraphEdge*, 
					   std::unordered_map<GraphEdge*, int>> outConnections;

	//all good at the beginning
	for (auto& edge : _graph.iterEdges())
	{
		edge->repetitive = false;
	}

	//first, by coverage
	for (auto& edge : _graph.iterEdges())
	{
		if (!edge->edgeId.strand()) continue;

		GraphEdge* complEdge = _graph.complementEdge(edge);
		edge->repetitive = false;
		complEdge->repetitive = false;

		if (edge->meanCoverage > _multInf.getUniqueCovThreshold() * 2 ||
		   (edge->isLooped() && edge->length() < Parameters::get().minimumOverlap))
		{
			edge->repetitive = true;
			complEdge->repetitive = true;
		}
	}

	//then, by read alignments
	for (auto& readPath : _readAlignments)
	{
		GraphAlignment filteredPath = readPath;
		if (filteredPath.size() < 2) continue;

		for (size_t i = 0; i < filteredPath.size() - 1; ++i)
		{
			if (filteredPath[i].edge == filteredPath[i + 1].edge &&
				filteredPath[i].edge->isLooped()) continue;

			GraphEdge* complLeft = _graph.complementEdge(filteredPath[i].edge);
			GraphEdge* complRight = _graph.complementEdge(filteredPath[i + 1].edge);
			++outConnections[filteredPath[i].edge][filteredPath[i + 1].edge];
			++outConnections[complRight][complLeft];
		}
	}
	
	for (auto& edgeList : outConnections)
	{
		if (edgeList.first-> multiplicity == 1 &&
			edgeList.second.size() == 1) continue;

		Logger::get().debug() << "Outputs: " << edgeList.first->edgeId.signedId()
			<< " " << edgeList.first->multiplicity;
		for (auto& outEdgeCount : edgeList.second)
		{
			Logger::get().debug() << "\t" << outEdgeCount.first->edgeId.signedId()
				<< " " << outEdgeCount.second << " " << outEdgeCount.first->isLooped();
		}
		Logger::get().debug() << "";
	}

	auto edgeMultiplicity = [this, &outConnections] (GraphEdge* edge)
	{
		int maxSupport = 0;
		int maxCoverage = 0;
		for (auto& outConn : outConnections[edge])
		{
			if (maxSupport < outConn.second)
			{
				maxSupport =  outConn.second;
				maxCoverage = outConn.first->meanCoverage;
			}
		}

		int multiplicity = 0;
		int minSupport = maxSupport / Constants::outPathsRatio;
		int minCoverage = maxCoverage / Constants::outPathsRatio;
		for (auto& outConn : outConnections[edge]) 
		{
			if (outConn.second > minSupport && 
				outConn.first->meanCoverage > minCoverage) ++multiplicity;
		}
		return multiplicity;
	};

	for (auto& edge : _graph.iterEdges())
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

		///////
		bool match = (edge->multiplicity > 1) == (edge->repetitive);
		if (!match && edge->multiplicity != 0)
		{
			std::string star = edge->repetitive ? "R" : " ";
			std::string loop = edge->isLooped() ? "L" : " ";
			Logger::get().debug() << star << " " << loop << " " 
				<< edge->edgeId.signedId()
				<< " " << edge->multiplicity << " -> " << mult << " ("
				<< leftMult << "," << rightMult << ") " << edge->length() << "\t"
				<< edge->meanCoverage;
		}
		///////
	}

	//mark all unprocessed edges as repetitive
	for (auto& edge : _graph.iterEdges())
	{
		GraphEdge* complEdge = _graph.complementEdge(edge);
		if (!outConnections.count(edge) && !outConnections.count(complEdge)
			&& edge->length() < Parameters::get().minimumOverlap)
		{
			Logger::get().debug() << "Not updated: " << edge->edgeId.signedId();
			edge->repetitive = true;
			complEdge->repetitive = true;
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

		this->updateAlignments();
		this->findRepeats();
	}

	this->removeUnsupportedEdges();
	this->clearResolvedRepeats();
}

void RepeatResolver::updateAlignments()
{
	//removes alignments that are no longer supported by the graph
	std::vector<GraphAlignment> newAlignments;
	int split = 0;
	for (auto& aln : _readAlignments)
	{
		GraphAlignment curAlignment;
		for (size_t i = 0; i < aln.size() - 1; ++i)
		{
			curAlignment.push_back(aln[i]);

			if (aln[i].edge->nodeRight != aln[i + 1].edge->nodeLeft)
			{
				++split;
				newAlignments.push_back(curAlignment);
				curAlignment.clear();
			}
		}

		curAlignment.push_back(aln.back());
		newAlignments.push_back(curAlignment);
	}
	Logger::get().debug() << "Split " << split << " alignments";

	_readAlignments = newAlignments;

	//mark resolved repeats
	/*
	int determinedRepeats = 0;
	for (auto& edge : _graph.iterEdges())
	{
		if (!edge->isRepetitive()) continue;

		if ((edge->nodeRight->outEdges.size() == 1 &&
			!edge->nodeRight->outEdges.front()->repetitive) ||
			(edge->nodeLeft->inEdges.size() == 1 &&
			!edge->nodeLeft->inEdges.front()->repetitive))
		{
			++determinedRepeats;
			edge->repetitive = false;
		}
	}

	Logger::get().debug() << "Determined " << determinedRepeats << " repeats";
	return determinedRepeats;*/
}

std::vector<RepeatResolver::Connection> 
	RepeatResolver::getConnections()
{
	
	auto safeEdge = [this](GraphEdge* edge)
	{
		return !edge->isRepetitive() && 
			   edge->length() > Constants::maxSeparation;
	};

	int totalSafe = 0;
	for (GraphEdge* edge : _graph.iterEdges())
	{
		if (edge->edgeId.strand() && safeEdge(edge)) ++totalSafe;
	}
	Logger::get().debug() << "Total unique edges: " << totalSafe;

	std::vector<Connection> readConnections;
	for (auto& readPath : _readAlignments)
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
		if (edge->meanCoverage <= coverageThreshold)
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
		return edge->isRepetitive() && edge->resolved;
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

void RepeatResolver::alignReads()
{
	//std::ofstream alnDump("../alignment_dump.txt");

	//create database
	std::unordered_map<FastaRecord::Id, 
					   std::pair<GraphEdge*, SequenceSegment>> idToSegment;
	SequenceContainer pathsContainer;

	for (auto& edge : _graph.iterEdges())
	{
		if (!edge->edgeId.strand()) continue;

		for (auto& segment : edge->seqSegments)
		{
			size_t len = segment.end - segment.start;
			auto sequence = _asmSeqs.getSeq(segment.seqId)
										.substr(segment.start, len);
			auto& newRec = pathsContainer.addSequence(sequence, "");

			idToSegment[newRec.id] = {edge, segment};
			idToSegment[newRec.id.rc()] = {_graph.complementEdge(edge), 
										   segment.complement()};
		}
	}

	//index it and align reads
	VertexIndex pathsIndex(pathsContainer);
	pathsIndex.countKmers(1);
	pathsIndex.buildIndex(1, Constants::readAlignMaxKmer, 
						  Constants::readAlignKmerSample);
	OverlapDetector readsOverlapper(pathsContainer, pathsIndex, 
									Constants::maximumJump,
									Constants::maxSeparation, /*overhang*/ 0);
	OverlapContainer readsOverlaps(readsOverlapper, _readSeqs, 
								   /*onlyMax*/ false);

	std::vector<FastaRecord::Id> allQueries;
	int64_t totalLength = 0;
	for (auto& readId : _readSeqs.getIndex())
	{
		if (_readSeqs.seqLen(readId.first) > Constants::maxSeparation)
		{
			totalLength += _readSeqs.seqLen(readId.first);
			allQueries.push_back(readId.first);
		}
	}
	std::mutex indexMutex;
	int numAligned = 0;
	int64_t alignedLength = 0;
	std::function<void(const FastaRecord::Id&)> alignRead = 
	[this, &indexMutex, &numAligned, &readsOverlaps, 
		&idToSegment, &pathsContainer, &alignedLength] 
		(const FastaRecord::Id& seqId)
	{
		auto overlaps = readsOverlaps.seqOverlaps(seqId);
		std::vector<EdgeAlignment> alignments;
		for (auto& ovlp : overlaps)
		{
			alignments.push_back({ovlp, idToSegment[ovlp.extId].first,
								  idToSegment[ovlp.extId].second});
		}
		std::sort(alignments.begin(), alignments.end(),
		  [](const EdgeAlignment& e1, const EdgeAlignment& e2)
			{return e1.overlap.curBegin < e2.overlap.curBegin;});
		auto readChain = this->chainReadAlignments(pathsContainer, alignments);

		if (readChain.empty()) return;
		indexMutex.lock();
		_readAlignments.push_back(readChain);
		++numAligned;
		alignedLength += readChain.back().overlap.curEnd - 
						 readChain.front().overlap.curBegin;
		indexMutex.unlock();
	};

	processInParallel(allQueries, alignRead, 
					  Parameters::get().numThreads, true);

	Logger::get().debug() << "Aligned " << numAligned << " / " 
		<< allQueries.size();
	Logger::get().debug() << "Aligned length " << alignedLength << " / " 
		<< totalLength << " " << (float)alignedLength / totalLength;
}

namespace
{
	struct Chain
	{
		Chain(): score(0) {}
		std::vector<const EdgeAlignment*> aln;
		int32_t score;
	};
}

GraphAlignment
	RepeatResolver::chainReadAlignments(const SequenceContainer& edgeSeqs,
								 	    const std::vector<EdgeAlignment>& ovlps) const
{
	std::list<Chain> activeChains;
	for (auto& edgeAlignment : ovlps)
	{
		std::list<Chain> newChains;
		int32_t maxScore = 0;
		Chain* maxChain = nullptr;
		for (auto& chain : activeChains)
		{
			const OverlapRange& nextOvlp = edgeAlignment.overlap;
			const OverlapRange& prevOvlp = chain.aln.back()->overlap;

			int32_t readDiff = nextOvlp.curBegin - prevOvlp.curEnd;
			int32_t graphDiff = nextOvlp.extBegin +
								prevOvlp.extLen - prevOvlp.extEnd;
			int32_t maxDiscordance = std::max(Constants::maximumJump / 
													Constants::farJumpRate,
											  Constants::maxSeparation);

			if (Constants::maximumJump > readDiff && 
					readDiff > Constants::alnOverlap &&
				Constants::maximumJump > graphDiff && 
					graphDiff > Constants::alnOverlap  &&
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
		for (auto& aln : maxChain->aln) result.push_back(*aln);
	}

	return result;
}

