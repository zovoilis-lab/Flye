//(c) 2016 by Authors
//This file is a part of ABruijn program.
//Released under the BSD license (see LICENSE file)

#include <deque>

#include "graph_processing.h"
#include "logger.h"

void GraphProcessor::simplify()
{
	this->trimTips();

	this->unrollLoops();
	this->condenceEdges();
	this->unrollLoops();
	this->condenceEdges();

	//this->updateEdgesMultiplicity();
}

void GraphProcessor::unrollLoops()
{
	auto unrollEdge = [](GraphEdge& loopEdge)
	{
		GraphEdge* prevEdge = nullptr;
		for (auto& inEdge : loopEdge.nodeLeft->inEdges)
		{
			if (inEdge != &loopEdge) prevEdge = inEdge;
		}
		GraphEdge* nextEdge = nullptr;
		for (auto& outEdge : loopEdge.nodeLeft->outEdges)
		{
			if (outEdge != &loopEdge) nextEdge = outEdge;
		}

		/*
		Logger::get().debug() << "In";
		for (auto seg : prevEdge->seqSegments)
		{
			Logger::get().debug() << seg.seqId << " " << seg.start << " " << seg.end;
		}
		Logger::get().debug() << "Out";
		for (auto seg : nextEdge->seqSegments)
		{
			Logger::get().debug() << seg.seqId << " " << seg.start << " " << seg.end;
		}*/

		auto growingSeqs = prevEdge->seqSegments;
		std::vector<SequenceSegment> updatedSeqs;
		for (auto& seq : growingSeqs)
		{
			while (true)
			{
				bool updated = false;
				for (auto& otherSeq : loopEdge.seqSegments)
				{
					if (seq.seqId != otherSeq.seqId ||
						seq.end != otherSeq.start) continue;

					updated = true;
					seq.end = otherSeq.end;
					break;
				}

				bool complete = false;
				for (auto& otherSeq : nextEdge->seqSegments)
				{
					if (seq.seqId != otherSeq.seqId ||
						seq.end != otherSeq.start) continue;

					complete = true;
					break;
				}

				if (complete)
				{
					updatedSeqs.push_back(seq);
					break;
				}
				if (!updated) break;
			}
		}
		if (!updatedSeqs.empty())
		{
			prevEdge->seqSegments = updatedSeqs;
			return true;
		}
		else
		{
			Logger::get().debug() << "Can't unroll";
			return false;
		}
	};

	std::unordered_set<GraphEdge*> toRemove;
	for (GraphEdge* edge : _graph.iterEdges())
	{
		if (edge->isLooped() &&
			edge->nodeLeft->inEdges.size() == 2 &&
			edge->nodeLeft->outEdges.size() == 2 &&
			edge->nodeLeft->neighbors().size() == 2)
		{
			if (unrollEdge(*edge))
			{
				Logger::get().debug() << "Unroll " << edge->edgeId.signedId();
				toRemove.insert(edge);
			}
		}
	}

	for (auto edge : toRemove)	{_graph.removeEdge(edge);};
}

void GraphProcessor::trimTips()
{
	std::unordered_set<GraphEdge*> toRemove;
	std::unordered_set<GraphEdge*> unknownEdges;
	for (GraphEdge* tipEdge : _graph.iterEdges())
	{
		int leftDegree = tipEdge->nodeLeft->inEdges.size();
		int rightDegree = tipEdge->nodeRight->outEdges.size();
		if (tipEdge->length() < _tipThreshold && 
			(leftDegree == 0 || rightDegree == 0))
		{
			toRemove.insert(tipEdge);
		}
	}

	Logger::get().debug() << toRemove.size() << " tips removed";
	for (auto edge : toRemove)	{_graph.removeEdge(edge);};
}

void GraphProcessor::condenceEdges()
{
	int edgesRemoved = 0;
	int edgesAdded = 0;

	auto collapseEdges = [] (const GraphPath& edges)
	{
		std::vector<GraphEdge> newEdges;
		std::list<SequenceSegment> growingSeqs(edges.front()->seqSegments.begin(),
											   edges.front()->seqSegments.end());
		assert(edges.size() > 1);
		size_t prevStart = 0;
		for (size_t i = 1; i < edges.size(); ++i)
		{
			auto prevSeqs = growingSeqs;
			for (auto prevSeg = growingSeqs.begin(); 
				 prevSeg != growingSeqs.end(); )
			{
				bool continued = false;
				for (auto& nextSeg : edges[i]->seqSegments)
				{
					if (prevSeg->seqId == nextSeg.seqId &&
						prevSeg->end == nextSeg.start)
					{
						continued = true;
						prevSeg->end = nextSeg.end;
					}
				}
				if (!continued)
				{
					prevSeg = growingSeqs.erase(prevSeg);
				}
				else
				{
					++prevSeg;
				}
			}

			if (growingSeqs.empty())
			{
				newEdges.emplace_back(edges[prevStart]->nodeLeft, 
									  edges[i - 1]->nodeRight);
				std::copy(prevSeqs.begin(), prevSeqs.end(),
				  		  std::back_inserter(newEdges.back().seqSegments));
				newEdges.back().multiplicity = edges[prevStart]->multiplicity;

				std::copy(edges[i]->seqSegments.begin(), 
						  edges[i]->seqSegments.end(), 
						  std::back_inserter(growingSeqs));
				prevStart = i;
			}
		}

		newEdges.emplace_back(edges[prevStart]->nodeLeft, 
							  edges.back()->nodeRight);
		std::copy(growingSeqs.begin(), growingSeqs.end(),
				  std::back_inserter(newEdges.back().seqSegments));
		newEdges.back().multiplicity = edges[prevStart]->multiplicity;

		return newEdges;
	};

	std::vector<GraphPath> toCollapse;
	std::unordered_set<FastaRecord::Id> usedDirections;
	for (auto& node : _graph.iterNodes())
	{
		if (!node->isBifurcation()) continue;

		for (auto& direction : node->outEdges)
		{
			if (usedDirections.count(direction->edgeId)) continue;
			usedDirections.insert(direction->edgeId);

			GraphNode* curNode = direction->nodeRight;
			GraphPath traversed;
			traversed.push_back(direction);
			while (!curNode->isBifurcation() &&
				   !curNode->outEdges.empty())
			{
				traversed.push_back(curNode->outEdges.front());
				curNode = curNode->outEdges.front()->nodeRight;
			}
			usedDirections.insert(traversed.back()->edgeId.rc());
			
			if (traversed.size() > 1)
			{
				toCollapse.emplace_back(std::move(traversed));
			}
		}
	}

	for (auto& unbranchingPath : toCollapse)
	{
		GraphPath complPath = _graph.complementPath(unbranchingPath);
		auto newEdges = collapseEdges(unbranchingPath);
		for (auto& edge : newEdges)
		{
			GraphEdge addFwd = edge;
			addFwd.edgeId = FastaRecord::Id(_graph._nextEdgeId);

			GraphEdge addRev(_graph.complementNode(edge.nodeRight),
							 _graph.complementNode(edge.nodeLeft),
							 FastaRecord::Id(_graph._nextEdgeId + 1));
			addRev.multiplicity = addFwd.multiplicity;

			//complementary segments
			for (auto seqSeg : addFwd.seqSegments)
			{
				int32_t seqLen = _asmSeqs.seqLen(seqSeg.seqId);
				addRev.seqSegments.emplace_back(seqSeg.seqId.rc(),
												seqLen - seqSeg.end - 1,
												seqLen - seqSeg.start - 1);
			}

			_graph.addEdge(std::move(addFwd));
			_graph.addEdge(std::move(addRev));
			_graph._nextEdgeId += 2;
		}

		for (auto edge : unbranchingPath) _graph.removeEdge(edge);
		for (auto edge : complPath) _graph.removeEdge(edge);

		edgesRemoved += unbranchingPath.size();
		edgesAdded += newEdges.size();
	}

	Logger::get().debug() << "Removed " << edgesRemoved << " edges";
	Logger::get().debug() << "Added " << edgesAdded << " edges";
}

void GraphProcessor::generateContigs()
{
	std::unordered_map<FastaRecord::Id, size_t> edgeIds;
	size_t nextEdgeId = 0;
	auto pathToId = [&edgeIds, &nextEdgeId](GraphPath path)
	{
		if (!edgeIds.count(path.front()->edgeId))
		{
			for (auto edge : path)
			{
				edgeIds[edge->edgeId] = nextEdgeId;
				edgeIds[edge->edgeId.rc()] = nextEdgeId + 1;
			}
			nextEdgeId += 2;
		}
		return FastaRecord::Id(edgeIds[path.front()->edgeId]);
	};
	
	std::unordered_set<GraphEdge*> visitedEdges;
	for (auto edge : _graph.iterEdges())
	{
		if (visitedEdges.count(edge)) continue;
		visitedEdges.insert(edge);

		GraphPath traversed;
		GraphNode* curNode = edge->nodeLeft;
		while (!curNode->isBifurcation() &&
			   !curNode->inEdges.empty() &&
			   !visitedEdges.count(curNode->inEdges.front()))
		{
			traversed.push_back(curNode->inEdges.front());
			visitedEdges.insert(traversed.back());
			curNode = curNode->inEdges.front()->nodeLeft;
		}

		std::reverse(traversed.begin(), traversed.end());
		traversed.push_back(edge);

		curNode = edge->nodeRight;
		while (!curNode->isBifurcation() &&
			   !curNode->outEdges.empty() &&
			   !visitedEdges.count(curNode->outEdges.front()))
		{
			traversed.push_back(curNode->outEdges.front());
			visitedEdges.insert(traversed.back());
			curNode = curNode->outEdges.front()->nodeRight;
		}

		FastaRecord::Id edgeId = pathToId(traversed);
		_contigs.emplace_back(traversed, edgeId);

		/*Logger::get().debug() << "Contig " << edgeId.signedId();
		for (auto& edge : traversed)
		{
			Logger::get().debug() << "\t" << edge->edgeId.signedId();
		}*/
	}
	Logger::get().info() << "Generated " << _contigs.size() / 2 << " contigs";
}

void GraphProcessor::outputContigsFasta(const std::string& filename)
{
	static const size_t FASTA_SLICE = 80;

	std::ofstream fout(filename);
	if (!fout.is_open()) throw std::runtime_error("Can't open " + filename);
	
	for (auto& contig : _contigs)
	{
		if (!contig.id.strand()) continue;

		std::string contigSequence;
		for (auto& edge : contig.path) 
		{
			if (!edge->seqSegments.empty())
			{
				const SequenceSegment& seg = edge->seqSegments.front();
				const std::string& edgeSeq = !edge->readSequence ? 
											 _asmSeqs.getSeq(seg.seqId) : 
											 _readSeqs.getSeq(seg.seqId);
				contigSequence += edgeSeq.substr(seg.start, seg.end - seg.start);
			}
			else
			{
				Logger::get().warning() << "Edge without sequence!";
			}
		}

		fout << ">edge_" << contig.id.signedId() << std::endl;
		for (size_t c = 0; c < contigSequence.length(); c += FASTA_SLICE)
		{
			fout << contigSequence.substr(c, FASTA_SLICE) << std::endl;
		}
	}
}

void GraphProcessor::outputContigsGraph(const std::string& filename)
{
	std::ofstream fout(filename);
	if (!fout.is_open()) throw std::runtime_error("Can't open " + filename);

	fout << "digraph {\n";
	
	///re-enumerating helper functions
	std::unordered_map<GraphNode*, int> nodeIds;
	int nextNodeId = 0;
	auto nodeToId = [&nodeIds, &nextNodeId](GraphNode* node)
	{
		if (!nodeIds.count(node))
		{
			nodeIds[node] = nextNodeId++;
		}
		return nodeIds[node];
	};

	const std::string COLORS[] = {"red", "darkgreen", "blue", "goldenrod", 
								  "cadetblue", "darkorchid", "aquamarine1", 
								  "darkgoldenrod1", "deepskyblue1", 
								  "darkolivegreen3"};
	std::unordered_map<FastaRecord::Id, size_t> colorIds;
	size_t nextColorId = 0;
	auto idToColor = [&colorIds, &nextColorId, &COLORS](FastaRecord::Id id)
	{
		if (!id.strand()) id = id.rc();
		if (!colorIds.count(id))
		{
			colorIds[id] = nextColorId;
			nextColorId = (nextColorId + 1) % 10;
		}
		return COLORS[colorIds[id]];
	};
	/////////////

	for (auto& contig : _contigs)
	{
		int32_t contigLength = 0;
		for (auto& edge : contig.path) contigLength += edge->length();

		if (contig.path.front()->isRepetitive())
		{
			std::string color = idToColor(contig.id);

			fout << "\"" << nodeToId(contig.path.front()->nodeLeft) 
				 << "\" -> \"" << nodeToId(contig.path.back()->nodeRight)
				 << "\" [label = \"" << contig.id.signedId() << 
				 " " << contigLength << " (" 
				 << contig.path.front()->multiplicity << ")\", color = \"" 
				 << color << "\" " << " penwidth = 3] ;\n";
		}
		else
		{
			fout << "\"" << nodeToId(contig.path.front()->nodeLeft) 
				 << "\" -> \"" << nodeToId(contig.path.back()->nodeRight)
				 << "\" [label = \"" << contig.id.signedId()
				 << " " << contigLength << "\", color = \"black\"] ;\n";
		}
	}

	fout << "}\n";

}
