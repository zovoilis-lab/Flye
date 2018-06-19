//(c) 2016 by Authors
//This file is a part of ABruijn program.
//Released under the BSD license (see LICENSE file)

#include <cmath>


#include "repeat_resolver.h"
#include "graph_processing.h"
#include "../common/config.h"
#include "../common/utils.h"
#include "../common/parallel.h"

#include <lemon/list_graph.h>
#include <lemon/matching.h>

//Given the path in the graph with a resolved repeat inside,
//separates in into a single unbranching path. The first
//and the last edges of the graphPath parameter
//should correspond to the flanking unique edges
void RepeatResolver::separatePath(const GraphPath& graphPath, 
								  SequenceSegment readSegment, 
								  FastaRecord::Id newId)
{
	//first edge
	GraphNode* leftNode = _graph.addNode();
	vecRemove(graphPath.front()->nodeRight->inEdges, graphPath.front());
	graphPath.front()->nodeRight = leftNode;
	leftNode->inEdges.push_back(graphPath.front());
	//int32_t pathCoverage = (graphPath.front()->meanCoverage +
	//					    graphPath.back()->meanCoverage) / 2;

	//repetitive edges in the middle
	for (size_t i = 1; i < graphPath.size() - 1; ++i)
	{
		graphPath[i]->resolved = true;
		//--graphPath[i]->multiplicity;
		//graphPath[i]->meanCoverage = 
		//	std::max(graphPath[i]->meanCoverage - pathCoverage, 0);
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

//Resolves all repeats simulateously through the graph mathcing optimization,
//Given the reads connecting unique edges (or pairs of edges in the transitions graph)
int RepeatResolver::resolveConnections(const std::vector<Connection>& connections, 
									   float minSupport)
{
	//Constructs transitions graph using the lemon library
	std::unordered_map<FastaRecord::Id, int> leftCoverage;
	std::unordered_map<FastaRecord::Id, int> rightCoverage;

	std::unordered_map<FastaRecord::Id, int> asmToLemon;
	std::unordered_map<int, FastaRecord::Id> lemonToAsm;
	lemon::ListGraph graph;
	lemon::ListGraph::EdgeMap<int> edgeWeights(graph);

	auto getEdge = [&graph](lemon::ListGraph::Node n1, lemon::ListGraph::Node n2)
	{
		for (lemon::ListGraph::IncEdgeIt edgeIt(graph, n1); 
			 edgeIt != lemon::INVALID; ++edgeIt) 
		{
			if (graph.oppositeNode(n1, edgeIt) == n2) return edgeIt;
		}
		return lemon::ListGraph::IncEdgeIt(lemon::INVALID);
	};

	for (auto& conn : connections)
	{
		GraphEdge* leftEdge = conn.path.front();
		GraphEdge* rightEdge = conn.path.back();

		if (leftEdge->edgeId == rightEdge->edgeId ||
			leftEdge->edgeId == rightEdge->edgeId.rc()) continue;

		++leftCoverage[leftEdge->edgeId];
		++rightCoverage[rightEdge->edgeId.rc()];

		if (!asmToLemon.count(leftEdge->edgeId))
		{
			auto newNode = graph.addNode();
			asmToLemon[leftEdge->edgeId] = graph.id(newNode);
			lemonToAsm[graph.id(newNode)] = leftEdge->edgeId;
		}
		if (!asmToLemon.count(rightEdge->edgeId.rc()))
		{
			auto newNode = graph.addNode();
			asmToLemon[rightEdge->edgeId.rc()] = graph.id(newNode);
			lemonToAsm[graph.id(newNode)] = rightEdge->edgeId.rc();
		}

		auto leftLemonNode = graph.nodeFromId(asmToLemon[leftEdge->edgeId]);
		auto rightLemonNode = graph.nodeFromId(asmToLemon[rightEdge->edgeId.rc()]);
		if (!graph.valid(getEdge(leftLemonNode, rightLemonNode)))
		{
			auto edge = graph.addEdge(leftLemonNode, rightLemonNode);
			edgeWeights[edge] = 0;
		}
		auto edge = getEdge(leftLemonNode, rightLemonNode);
		++edgeWeights[edge];
	}

	//copmutes maximum weight matching on this graph
	lemon::MaxWeightedMatching<lemon::ListGraph> matcher(graph, edgeWeights);
	matcher.run();

	//converting matching to the resolved paths on the graph
	std::unordered_set<FastaRecord::Id> usedEdges;
	std::vector<Connection> uniqueConnections;
	int unresolvedLinks = 0;
	for (auto lemonAsm : lemonToAsm)
	{
		auto mateNode = matcher.mate(graph.nodeFromId(lemonAsm.first));
		if (mateNode == lemon::INVALID) continue;

		FastaRecord::Id leftId = lemonAsm.second;
		FastaRecord::Id rightId = lemonToAsm[graph.id(mateNode)];
		int support = edgeWeights[getEdge(graph.nodeFromId(lemonAsm.first), 
										  mateNode)];

		if (usedEdges.count(leftId)) continue;
		usedEdges.insert(rightId);

		float confidence = (float)support / (leftCoverage[leftId] + 
									  		 rightCoverage[rightId]);

		Logger::get().debug() << "\tConnection " 
			<< leftId.signedId() << "\t" << rightId.rc().signedId()
			<< "\t" << support / 4 << "\t" << confidence;

		if (confidence < minSupport)
		{
			++unresolvedLinks;
			continue;
		}

		std::vector<Connection> spanningConnections;
		for (auto& conn : connections)
		{
			if ((conn.path.front()->edgeId == leftId && 
				 	conn.path.back()->edgeId == rightId.rc()) ||
				(conn.path.front()->edgeId == rightId && 
				 	conn.path.back()->edgeId == leftId.rc()))
			{
				spanningConnections.push_back(conn);
			}
		}
		std::sort(spanningConnections.begin(), spanningConnections.end(),
				  [](const Connection c1, const Connection c2)
				{return c1.readSequence.length() < c2.readSequence.length();});
		uniqueConnections
			.push_back(spanningConnections[spanningConnections.size() / 2]);
	}

	//separates the resolved paths in the graph
	for (auto& conn : uniqueConnections)
	{
		GraphPath complPath = _graph.complementPath(conn.path);
		SequenceSegment complSegment = conn.readSequence.complement();

		FastaRecord::Id edgeId = _graph.newEdgeId();
		this->separatePath(conn.path, conn.readSequence, edgeId);
		this->separatePath(complPath, complSegment, edgeId.rc());
	}

	Logger::get().debug() << "Resolved: " << uniqueConnections.size() << " links: "
						  << connections.size() / 2;
	Logger::get().debug() << "Unresolved: " << unresolvedLinks;

	return uniqueConnections.size();
}

//Classifies all edges into unique and repetitive based on the coverage + 
//alignment information - one of the key steps here.
void RepeatResolver::findRepeats()
{
	Logger::get().debug() << "Finding repeats";

	//all edges are unique at the beginning
	for (auto& edge : _graph.iterEdges())
	{
		edge->repetitive = false;
	}

	//Will operate on unbranching paths rather than single edges
	GraphProcessor proc(_graph, _asmSeqs, _readSeqs);
	auto unbranchingPaths = proc.getUnbranchingPaths();
	std::unordered_map<FastaRecord::Id, UnbranchingPath*> idToPath;
	for (auto& path : unbranchingPaths) idToPath[path.id] = &path;
	auto complPath = [&idToPath](UnbranchingPath* path)
	{
		if (idToPath.count(path->id.rc()))
		{
			return idToPath[path->id.rc()];
		}
		return path;	//self-complement
	};
	auto markRepetitive = [](UnbranchingPath* path)
	{
		for (auto& edge : path->path) edge->repetitive = true;
	};

	for (auto& path : unbranchingPaths)
	{
		if (!path.id.strand()) continue;

		//mark paths with high coverage as repetitive
		if (path.meanCoverage > _multInf.getUniqueCovThreshold() &&
			path.length < (int)Config::get("unique_edge_length"))
		{
			markRepetitive(&path);
			markRepetitive(complPath(&path));
			Logger::get().debug() << "Cov: " 
				<< path.edgesStr() << "\t" << path.length << "\t" 
				<< path.meanCoverage;
		}

		//self-complements are repetitive
		for (auto& edge : path.path)
		{
			if (edge->selfComplement)
			{
				Logger::get().debug() << "Self-compl: " << path.edgesStr();
				markRepetitive(&path);
				markRepetitive(complPath(&path));
				break;
			}
		}

		//tandem repeats
		if (path.path.size() == 1 && path.path.front()->isLooped())
		{
			std::unordered_set<FastaRecord::Id> seen;
			for (auto& seg : path.path.front()->seqSegments)
			{
				if (seen.count(seg.seqId))
				{
					markRepetitive(&path);
					markRepetitive(complPath(&path));
					Logger::get().debug() << "Tan: " 
						<< path.edgesStr() << "\t" << path.length << "\t" 
						<< path.meanCoverage;
					break;
				}
				seen.insert(seg.seqId);
			}
		}
	}

	//Finally, using the read alignments
	std::unordered_map<GraphEdge*, 
					   std::unordered_map<GraphEdge*, int>> outConnections;
	for (auto& readPath : _aligner.getAlignments())
	{
		if (readPath.size() < 2) continue;
		//int overhang = std::max(readPath.front().overlap.curBegin,
		//						readPath.back().overlap.curLen - 
		//							readPath.back().overlap.curEnd);
		//if (overhang > (int)Config::get("maximum_overhang")) continue;

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
	//computes the number of successors of the given edge in read alignments
	auto rightMultiplicity = [this, &outConnections] (GraphEdge* edge)
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
		int minSupport = maxSupport / (int)Config::get("out_paths_ratio");
		for (auto& outConn : outConnections[edge]) 
		{
			if (outConn.second > minSupport)
			{
				//all repeatitive successors are count as one
				outConn.first->repetitive ? ++repeatMult : ++uniqueMult;
			}
		}
		return std::min(repeatMult, 1) + uniqueMult;
	};

	//order might be important, process short edges first
	std::vector<UnbranchingPath*> sortedPaths;
	for (auto& path : unbranchingPaths) sortedPaths.push_back(&path);
	std::sort(sortedPaths.begin(), sortedPaths.end(),
			  [](const UnbranchingPath* p1, const UnbranchingPath* p2) 
			  {return p1->length < p2->length;});
	for (auto& path : sortedPaths)
	{
		if (!path->id.strand()) continue;
		if (path->path.front()->repetitive) continue;

		int rightMult = rightMultiplicity(path->path.back());
		int leftMult = rightMultiplicity(complPath(path)->path.back());
		int mult = std::max(leftMult, rightMult);
		if (mult > 1) 
		{
			markRepetitive(path);
			markRepetitive(complPath(path));

			Logger::get().debug() << "Mult: " 
				<< path->edgesStr() << "\t" << path->length << "\t" 
				<< path->meanCoverage << "\t" << mult << " ("
				<< leftMult << "," << rightMult << ")";
			
			for (auto& outEdgeCount : outConnections[path->path.back()])
			{
				std::string star = outEdgeCount.first->repetitive ? "R" : " ";
				std::string loop = outEdgeCount.first->isLooped() ? "L" : " ";
				std::string tip = outEdgeCount.first->isTip() ? "T" : " ";
				Logger::get().debug() << "+\t" << star << " " << loop << " " << tip << " "
					<< outEdgeCount.first->edgeId.signedId() << "\t" << outEdgeCount.second;
			}
			for (auto& outEdgeCount : outConnections[complPath(path)->path.back()])
			{
				std::string star = outEdgeCount.first->repetitive ? "R" : " ";
				std::string loop = outEdgeCount.first->isLooped() ? "L" : " ";
				std::string tip = outEdgeCount.first->isTip() ? "T" : " ";
				Logger::get().debug() << "-\t" << star << " " << loop << " " << tip << " "
					<< outEdgeCount.first->edgeId.signedId() << "\t" << outEdgeCount.second;
			}
		}
	}

	//now, check for this structure >-<, in case read alignments were not enough
	for (auto& path : sortedPaths)
	{
		if (path->path.front()->repetitive || path->isLooped()) continue;
		if (path->path.front()->nodeLeft->outEdges.size() > 1 ||
			path->path.back()->nodeRight->inEdges.size() > 1) continue;

		int numIn = 0;
		int numOut = 0;
		for (auto& edge: path->path.front()->nodeLeft->inEdges)
		{
			if (!edge->repetitive) ++numIn;
		}
		for (auto& edge: path->path.back()->nodeRight->outEdges)
		{
			if (!edge->repetitive) ++numOut;
		}
		if (numIn > 1 && numOut > 1)
		{
			Logger::get().debug() << "Structure: " << path->edgesStr();
			markRepetitive(path);
			markRepetitive(complPath(path));
		}
	}
}

void RepeatResolver::fixLongEdges()
{
	GraphProcessor proc(_graph, _asmSeqs, _readSeqs);
	auto unbranchingPaths = proc.getUnbranchingPaths();
	for (auto& path : unbranchingPaths)
	{
		if (!path.id.strand()) continue;

		if (!path.path.front()->selfComplement &&
			path.path.front()->repetitive &&
			path.length > (int)Config::get("unique_edge_length") &&
			(float)path.meanCoverage < _multInf.getUniqueCovThreshold())
		{
			for (auto& edge : path.path)
			{
				edge->repetitive = false;
				_graph.complementEdge(edge)->repetitive = false;
			}

			Logger::get().debug() << "Fixed: " 
				<< path.edgesStr() << "\t" << path.length << "\t" 
				<< path.meanCoverage;
		}
	}
}

//Iterates repeat detection and resolution until
//no new repeats are resolved
void RepeatResolver::resolveRepeats()
{
	bool perfectIter = true;
	while (true)
	{
		auto connections = this->getConnections();

		float minSupport = perfectIter ? 1.0f : 
						   Config::get("min_repeat_res_support");
		int resolvedConnections = 
			this->resolveConnections(connections, minSupport);
		perfectIter = !perfectIter;

		this->clearResolvedRepeats();
		_aligner.updateAlignments();

		if (!resolvedConnections && !perfectIter) break;
		this->findRepeats();
	}

	GraphProcessor proc(_graph, _asmSeqs, _readSeqs);
	proc.fixChimericJunctions();
	_aligner.updateAlignments();
}


//extracts conenctions between pairs of unique edges from
//read alignments
std::vector<RepeatResolver::Connection> 
	RepeatResolver::getConnections()
{
	
	auto safeEdge = [this](GraphEdge* edge)
	{
		return !edge->isRepetitive();
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
		GraphAlignment currentAln;
		int32_t readStart = 0;
		for (auto& aln : readPath)
		{
			if (currentAln.empty()) 
			{
				if (!safeEdge(aln.edge)) continue;
				readStart = aln.overlap.curEnd + aln.overlap.extLen - 
							aln.overlap.extEnd;
			}

			currentAln.push_back(aln);
			if (safeEdge(aln.edge))
			{
				if (currentAln.back().edge->nodeLeft->isBifurcation() ||
					currentAln.front().edge->nodeRight->isBifurcation()) 
				{

					int32_t flankScore = std::min(currentAln.front().overlap.curRange(),
												  currentAln.back().overlap.curRange());
					GraphPath currentPath;
					for (auto& aln : currentAln) currentPath.push_back(aln.edge);
					GraphPath complPath = _graph.complementPath(currentPath);

					int32_t readEnd = aln.overlap.curBegin - aln.overlap.extBegin;
					readEnd = std::max(readStart + 100, readEnd);	//TODO: less ad-hoc fix
					SequenceSegment segment(aln.overlap.curId, aln.overlap.curLen, 
											readStart, readEnd);
					segment.segType = SequenceSegment::Read;
					SequenceSegment complSegment = segment.complement();

					readConnections.push_back({currentPath, segment, flankScore});
					readConnections.push_back({complPath, complSegment, flankScore});
				}

				currentAln.clear();
				currentAln.push_back(aln);
				readStart = aln.overlap.curEnd + aln.overlap.extLen - 
							aln.overlap.extEnd;
			}
		}
	}

	return readConnections;
}

//cleans up the graph after repeat resolution
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
