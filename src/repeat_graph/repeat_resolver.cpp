//(c) 2016 by Authors
//This file is a part of ABruijn program.
//Released under the BSD license (see LICENSE file)

#include <cmath>

#include "repeat_resolver.h"
#include "config.h"
#include "utils.h"
#include "bipartie_mincost.h"

#include <simplex.h>
#include <variable.h>


RepeatResolver::GraphAlignment
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
		if (edgeAlignment.overlap.curBegin < _readOverhang)
		{
			activeChains.push_back({&edgeAlignment});
		}
	}

	int32_t maxSpan = 0;
	Chain* maxChain = nullptr;
	std::unordered_set<FastaRecord::Id> inEdges;
	std::unordered_set<FastaRecord::Id> outEdges;
	for (auto& chain : activeChains)
	{
		//check right overhang
		int32_t overhang = chain.back()->overlap.curLen - 
						   chain.back()->overlap.curEnd;
		if (overhang > _readOverhang) continue;

		int32_t readSpan = chain.back()->overlap.curEnd - 
						   chain.front()->overlap.curBegin;
		if (readSpan > maxSpan)
		{
			maxSpan = readSpan;
			maxChain = &chain;
		}

		inEdges.insert(chain.front()->edge->edgeId);
		outEdges.insert(chain.back()->edge->edgeId);
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

		if (inEdges.size() != 1 || outEdges.size() != 1) return {};

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
		newEdge->readSequence = true;
	}

	//last edge
	vecRemove(graphPath.back()->nodeLeft->outEdges, graphPath.back());
	graphPath.back()->nodeLeft = rightNode;
	rightNode->outEdges.push_back(graphPath.back());
}

void RepeatResolver::resolveConnections(const std::vector<Connection>& connections)
{
	Logger::get().debug() << "Resolving repeats";

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
			<< leftEdge->seqSegments.front().seqId
			<< "\t" << leftEdge->seqSegments.front().end << "\t"
			<< rightEdge->seqSegments.front().seqId
			<< "\t" << rightEdge->seqSegments.front().start
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
		int32_t complStart = _readSeqs.seqLen(conn.readSequence.seqId) - 
							 conn.readSequence.end - 1;
		int32_t complEnd = _readSeqs.seqLen(conn.readSequence.seqId) - 
							 conn.readSequence.start - 1;
		SequenceSegment complSegment(conn.readSequence.seqId.rc(), complStart, 
									 complEnd);

		this->separatePath(conn.path, conn.readSequence, _graph._nextEdgeId);
		this->separatePath(complPath, complSegment, _graph._nextEdgeId + 1);
		_graph._nextEdgeId += 2;
	}

	Logger::get().debug() << "Edges: " << totalLinks << " links: "
						  << connections.size();
}

void RepeatResolver::resolveRepeats()
{
	std::vector<Connection> readConnections;
	for (auto& readPath : _readAlignments)
	{
		/*
		if (!readPath.empty())
		{
			Logger::get().debug() << _readSeqs.seqName(readId.first);
		}*/

		GraphPath currentPath;
		int32_t readStart = 0;
		for (auto& aln : readPath)
		{
			if (currentPath.empty()) 
			{
				if (aln.edge->isRepetitive() || 
					aln.edge->length() < Constants::maximumJump) continue;
				readStart = aln.overlap.curEnd + aln.overlap.extLen - 
							aln.overlap.extEnd;
			}
			
			/*if (!aln.edge->isRepetitive())
			{
				Logger::get().debug() << aln.edge->edgeId << "\t" 
									  << aln.edge->seqSegments.front().seqId << "\t"
									  << aln.edge->seqSegments.front().start << "\t"
									  << aln.edge->seqSegments.front().end << "\t"
									  << aln.overlap.curRange();
			}*/

			currentPath.push_back(aln.edge);
			if (!aln.edge->isRepetitive() && currentPath.size() > 1 &&
				aln.edge->length() >= Constants::maximumJump)
			{
				GraphPath complPath = _graph.complementPath(currentPath);

				int32_t readEnd = aln.overlap.curBegin - aln.overlap.extBegin;
				int32_t complStart = aln.overlap.curLen - readEnd - 1;
				int32_t complEnd = aln.overlap.curLen - readStart - 1;
				SequenceSegment segment(aln.overlap.curId, readStart, readEnd);
				SequenceSegment complSegment(aln.overlap.curId.rc(), complStart,
											 complEnd);

				readConnections.push_back({currentPath, segment});
				readConnections.push_back({complPath, complSegment});

				currentPath.clear();
				currentPath.push_back(aln.edge);
				readStart = aln.overlap.curEnd + aln.overlap.extLen - 
							aln.overlap.extEnd;
			}
		}
	}

	this->resolveConnections(readConnections);
	this->clearResolvedRepeats();
}


void RepeatResolver::correctEdgesMultiplicity()
{
	std::unordered_map<GraphEdge*, int64_t> edgesCoverage;
	for (auto& path : _readAlignments)
	{
		for (size_t i = 0; i < path.size(); ++i)
		{
			if (0 < i && i < path.size() - 1)
			{
				edgesCoverage[path[i].edge] += path[i].edge->length();
			}
			else
			{
				edgesCoverage[path[i].edge] += path[i].overlap.extRange();
			}
		}
	}

	int64_t sumCov = 0;
	int64_t sumLength = 0;
	for (auto edgeCov : edgesCoverage)
	{
		sumCov += edgeCov.second;
		sumLength += edgeCov.first->length();
	}
	int meanCoverage = sumCov / sumLength;
	Logger::get().debug() << "Mean edge coverage: " << meanCoverage;

	for (auto edgeCov : edgesCoverage)
	{
		if (!edgeCov.first->edgeId.strand()) continue;
		if (edgeCov.first->isLooped() && 
			edgeCov.first->length() < Constants::maximumJump) continue;

		GraphEdge* complEdge = _graph.complementPath({edgeCov.first}).front();
		int normCov = (edgeCov.second + edgesCoverage[complEdge]) / 
									(2 * edgeCov.first->length());
		edgeCov.first->coverage = (float)normCov;
		complEdge->coverage = (float)normCov;

		int estMult = std::max(1.0f, roundf((float)normCov / meanCoverage));
		edgeCov.first->multiplicity = estMult;
		complEdge->multiplicity = estMult;
	}

	this->correctWeights();
}

void RepeatResolver::clearResolvedRepeats()
{
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
		//single looped nodes
		if (node->neighbors().size() == 0)
		{
			bool resolvedRepeat = true;
			for (auto& edge : node->outEdges) 
			{
				if (edge->multiplicity > 0) resolvedRepeat = false;
			}

			bool smallLoop = node->outEdges.size() == 1 && 
				node->outEdges.front()->isRepetitive() &&
				node->outEdges.front()->length() < Constants::maximumJump;

			if (resolvedRepeat || smallLoop) toRemove.insert(node);
		}

		//other nodes
		if (node->neighbors().size() != 1) continue;

		GraphEdge* direction = nextEdge(node);
		if (!direction) continue;

		GraphPath traversed;
		traversed.push_back(direction);
		GraphNode* curNode = direction->nodeRight;
		while (curNode->neighbors().size() == 2)
		{
			traversed.push_back(nextEdge(curNode));
			curNode = traversed.back()->nodeRight;
		}
		if (traversed.empty()) continue;

		bool removeLast = curNode->neighbors().size() == 1;
		bool resolvedRepeat = true;
		for (auto& edge : traversed) 
		{
			if (edge->multiplicity != 0) resolvedRepeat = false;
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

void RepeatResolver::correctWeights()
{
	using namespace optimization;

	std::unordered_map<FastaRecord::Id, size_t> edgeToId;
	size_t nextId = 0;

	auto edgeNumber = [&edgeToId, &nextId](GraphEdge* edge)
	{
		if (!edgeToId.count(edge->edgeId))
		{
			edgeToId[edge->edgeId] = nextId;
			edgeToId[edge->edgeId.rc()] = nextId;
			++nextId;
		}
		return edgeToId[edge->edgeId];
	};

	//enumerating
	for (auto edge : _graph.iterEdges())
	{
		if (!edge->edgeId.strand() || edge->isLooped()) continue;
		edgeNumber(edge);
	}

	//formulate linear programming
	Simplex simplex("");
	for (auto edge : _graph.iterEdges())
	{
		if (!edge->edgeId.strand() || edge->isLooped()) continue;

		size_t varNumber = edgeNumber(edge);

		Matrix eye(1, nextId, 0);
		eye(varNumber) = 1;

		std::string varName = std::to_string(edge->edgeId.signedId());
		simplex.add_variable(new Variable(&simplex, varName.c_str()));
	
		simplex.add_constraint(Constraint(eye, CT_MORE_EQUAL, 
										  (float)edge->multiplicity));    
		simplex.add_constraint(Constraint(eye, CT_LESS_EQUAL, 999.0f));
	}

	std::unordered_set<GraphNode*> processedNodes;
	for (auto node : _graph.iterNodes())
	{
		if (node->inEdges.empty() || node->outEdges.empty()) continue;
		if (node->neighbors().size() < 2) continue;
		if (processedNodes.count(node)) continue;

		GraphNode* complNode = _graph.complementNode(node);
		processedNodes.insert(complNode);

		std::vector<int> coefficients(nextId, 0);
		for (auto edge : node->inEdges) 
		{
			if (!edge->isLooped()) coefficients[edgeNumber(edge)] += 1;
		}
		for (auto edge : node->outEdges) 
		{
			if (!edge->isLooped()) coefficients[edgeNumber(edge)] -= 1;
		}

		//for (int x : coefficients) std::cout << x << "\t";
		//std::cout << "= 0\n";

        Matrix coefMatrix(1, nextId, 0);
		bool nonTrivial = false;
		for (size_t i = 0; i < nextId; ++i) 
		{
			coefMatrix(i) = coefficients[i];
			if (coefficients[i]) nonTrivial = true;
		}
		if (nonTrivial)
		{
        	simplex.add_constraint(Constraint(coefMatrix, CT_EQUAL, 0.0f));
		}
	}

    Matrix costs(1, nextId, 1.0f);
    simplex.set_objective_function(ObjectiveFunction(OFT_MINIMIZE, costs));   

	simplex.solve();
	if (!simplex.has_solutions() || simplex.must_be_fixed() || 
		simplex.is_unlimited()) throw std::runtime_error("Error while solving LP");

	//simplex.print_solution();
	for (auto edge : _graph.iterEdges())
	{
		if (!edge->isLooped())
		{
			edge->multiplicity = simplex.get_solution()(edgeNumber(edge));
		}
	}
}


void RepeatResolver::alignReads()
{
	//create database
	std::unordered_map<FastaRecord::Id, 
					   std::pair<GraphEdge*, SequenceSegment*>> idToSegment;
	SequenceContainer pathsContainer;

	for (auto& edge : _graph.iterEdges())
	{
		for (auto& segment : edge->seqSegments)
		{
			size_t len = segment.end - segment.start;
			std::string sequence = _asmSeqs.getSeq(segment.seqId)
												   .substr(segment.start, len);
			auto& newRec = pathsContainer.addSequence(sequence, "");
			idToSegment[newRec.id] = {edge, &segment};
			//idToSegment[newRec.id.rc()] = {edge, &segment};
		}
	}

	//index it and align reads
	VertexIndex pathsIndex(pathsContainer);
	pathsIndex.countKmers(1);
	pathsIndex.buildIndex(1, 5000, 1);
	OverlapDetector readsOverlapper(pathsContainer, pathsIndex, 
									_readJump, _maxSeparation, 0);
	OverlapContainer readsOverlaps(readsOverlapper, _readSeqs, false);
	readsOverlaps.findAllOverlaps();

	Logger::get().debug() << "Threading reads through the graph";
	//get connections
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
	}
}
