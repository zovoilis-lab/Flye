//(c) 2016 by Authors
//This file is a part of ABruijn program.
//Released under the BSD license (see LICENSE file)

#include <deque>
#include <iomanip>

#include "graph_processing.h"
#include "../sequence/contig_generator.h"
#include "../common/logger.h"
#include "../common/config.h"

void GraphProcessor::condence()
{
	this->trimTips();
	//this->trimFakeLoops();

	//this->unrollLoops();
	this->condenceEdges();
	//this->unrollLoops();
	//this->condenceEdges();

	this->fixChimericJunctions();
	this->trimTips();
}

/*
void GraphProcessor::trimFakeLoops()
{
	std::unordered_set<GraphEdge*> toRemove;
	for (auto& edge : _graph.iterEdges())
	{
		if (edge->isLooped() && 
			edge->length() < Parameters::get().minimumOverlap)
		{
			std::unordered_set<FastaRecord::Id> seenSeqs;
			bool isFake = true;
			for (auto& seg : edge->seqSegments)
			{
				if (seenSeqs.count(seg.seqId))
				{
					isFake = false;
					break;
				}
				seenSeqs.insert(seg.seqId);
			}

			if (isFake)
			{
				toRemove.insert(edge);
				toRemove.insert(_graph.complementEdge(edge));
			}
		}
	}

	for (auto& edge : toRemove)	_graph.removeEdge(edge);
	Logger::get().debug() << "Removed " << toRemove.size() / 2 
		<< " fake loops";
}*/

void GraphProcessor::fixChimericJunctions()
{
	//a very specific case: 1 in - 1 out
	std::unordered_set<GraphNode*> simpleCases;
	for (auto& node : _graph.iterNodes())
	{
		if (!node->isBifurcation() &&
			(node->inEdges.front()->edgeId ==
			 node->outEdges.front()->edgeId.rc()))
		{
			simpleCases.insert(node);
		}
	}
	for (auto& node : simpleCases)
	{
		GraphNode* newNode = _graph.addNode();
		GraphEdge* cutEdge = node->outEdges.front();
		newNode->outEdges.push_back(cutEdge);
		cutEdge->nodeLeft = newNode;
		node->outEdges.clear();
	}

	//more common case: 2 in - 2 out
	std::unordered_set<GraphNode*> complexCases;
	for (auto& node : _graph.iterNodes())
	{
		if (node->inEdges.size() != 2 ||
			node->outEdges.size() != 2) continue;
		auto& inEdges = node->inEdges;
		auto& outEdges = node->outEdges;
		if (inEdges[0]->edgeId.rc() != outEdges[0]->edgeId)
		{
			//match INs with OUTs
			std::swap(inEdges[0], inEdges[1]);
		}

		if (inEdges[0]->edgeId.rc() == outEdges[0]->edgeId &&
			inEdges[1]->edgeId.rc() == outEdges[1]->edgeId)
		{
			complexCases.insert(node);
		}
	}
	for(auto& node : complexCases)
	{
		GraphNode* newNode = _graph.addNode();
		node->inEdges[1]->nodeRight = newNode;
		node->outEdges[0]->nodeLeft = newNode;
		newNode->inEdges.push_back(node->inEdges[1]);
		newNode->outEdges.push_back(node->outEdges[0]);

		node->inEdges.pop_back();
		node->outEdges.erase(node->outEdges.begin());
	}

	Logger::get().debug() << "Removed " 
		<< simpleCases.size() + complexCases.size()
		<< " chimeric junctions";
}

/*
void GraphProcessor::unrollLoops()
{
	auto unrollEdge = [this](GraphEdge& loopEdge)
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
		if (updatedSeqs.size() == prevEdge->seqSegments.size())
		{
			prevEdge->seqSegments = updatedSeqs;
			//Logger::get().debug() << "Unroll " << loopEdge.edgeId.signedId();
			return true;

		}
		//Logger::get().debug() << "Can't unroll " << loopEdge.edgeId.signedId();
		return false;
	};

	int unrollLoops = 0;

	std::unordered_set<GraphEdge*> toRemove;
	for (GraphEdge* edge : _graph.iterEdges())
	{
		if (!edge->isLooped()) continue;
		if (edge->nodeLeft->inEdges.size() == 2 &&
			edge->nodeLeft->outEdges.size() == 2 &&
			edge->nodeLeft->neighbors().size() == 2)
		{
			if (unrollEdge(*edge))
			{
				++unrollLoops;
				//Logger::get().debug() << "Unroll " << edge->edgeId.signedId();
			}
			toRemove.insert(edge);
		}
	}
	for (auto& edge : toRemove)	_graph.removeEdge(edge);

	Logger::get().debug() << "Unrolled " << unrollLoops 
		<< ", removed " << toRemove.size();
}*/

void GraphProcessor::trimTips()
{
	std::unordered_set<GraphNode*> toRemove;
	for (GraphEdge* tipEdge : _graph.iterEdges())
	{
		if (tipEdge->length() < _tipThreshold && tipEdge->isTip())
		{
			int leftDegree = tipEdge->nodeLeft->inEdges.size();
			toRemove.insert(leftDegree == 0 ? tipEdge->nodeLeft : 
							tipEdge->nodeRight);
		}
	}

	Logger::get().debug() << toRemove.size() << " tips removed";
	for (auto& edge : toRemove)	_graph.removeNode(edge);
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

			std::unordered_set<FastaRecord::Id> complementEdges;
			complementEdges.insert(direction->edgeId.rc());

			while (!curNode->isBifurcation() &&
				   !curNode->outEdges.empty() &&
				   !complementEdges.count(curNode->outEdges.front()->edgeId))
			{
				traversed.push_back(curNode->outEdges.front());
				complementEdges.insert(traversed.back()->edgeId.rc());
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
			addFwd.edgeId = _graph.newEdgeId();

			GraphEdge addRev(_graph.complementNode(edge.nodeRight),
							 _graph.complementNode(edge.nodeLeft),
							 addFwd.edgeId.rc());

			//complementary segments
			for (auto seqSeg : addFwd.seqSegments)
			{
				addRev.seqSegments.push_back(seqSeg.complement());
			}

			_graph.addEdge(std::move(addFwd));
			_graph.addEdge(std::move(addRev));
		}

		std::unordered_set<GraphEdge*> toRemove;
		for (auto& edge : unbranchingPath) toRemove.insert(edge);
		for (auto& edge : complPath) toRemove.insert(edge);
		for (auto& edge : toRemove) _graph.removeEdge(edge);

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
		bool circular = (traversed.front()->nodeLeft == 
						traversed.back()->nodeRight) &&
						traversed.front()->nodeLeft->outEdges.size() == 1;

		bool repetitive = traversed.front()->isRepetitive() || 
						  traversed.back()->isRepetitive();

		int contigLength = 0;
		int64_t sumCov = 0;
		std::string contentsStr;
		for (auto& edge : traversed) 
		{
			contentsStr += std::to_string(edge->edgeId.signedId()) + " -> ";
			contigLength += edge->length();
			sumCov += edge->meanCoverage * edge->length();
		}
		int meanCoverage = contigLength ? sumCov / contigLength : 0;

		_contigs.emplace_back(traversed, edgeId, circular, 
							  contigLength, meanCoverage);
		_contigs.back().repetitive = repetitive;

		if (edgeId.strand())
		{
			Logger::get().debug() << "Contig " << edgeId.signedId() 
							<< ": " << contentsStr;
		}
	}
	this->generateContigSequences(_contigs);
	Logger::get().info() << "Generated " << _contigs.size() / 2 << " contigs";
}

void GraphProcessor::generateContigSequences(std::vector<Contig>& contigs) const
{
	ContigGenerator gen;
	for (auto& contig : contigs)
	{
		std::unordered_map<FastaRecord::Id, int> seqIdFreq;
		for (auto& edge : contig.path) 
		{
			std::unordered_set<FastaRecord::Id> edgeSeqIds;
			for (auto& seg: edge->seqSegments) 
			{
				edgeSeqIds.insert(seg.seqId);
			}
			for (auto& seqId : edgeSeqIds)
			{
				seqIdFreq[seqId] += 1;
			}
		}

		ContigPath contigPath;
		int32_t prevFlank = 0;
		int32_t prevLength = 0;

		//Logger::get().debug() << "Contig: " << contig.id.signedId();
		//int32_t contigStart = 0;
		//int32_t theoreticalLen = 0;
		for (size_t i = 0; i < contig.path.size(); ++i) 
		{
			if (contig.path[i]->seqSegments.empty()) 
			{
				throw std::runtime_error("Edge without sequence");
			}

			//get the sequence with maximum frequency
			SequenceSegment* bestSegment = nullptr;
			for (auto& seg : contig.path[i]->seqSegments)
			{
				if (!bestSegment || 
					seqIdFreq[seg.seqId] > seqIdFreq[bestSegment->seqId])
				{
					bestSegment = &seg;
				}

				//auto name = (!seg.readSequence) ? 
				//			_asmSeqs.seqName(seg.seqId) :
				//			_readSeqs.seqName(seg.seqId);
				//Logger::get().debug() << "\t\t" << name << "\t"
				//	<< seg.start << "\t" << seg.end;
			}

			//auto name = (!bestSegment->readSequence) ? 
			//				_asmSeqs.seqName(bestSegment->seqId) :
			//				_readSeqs.seqName(bestSegment->seqId);
			//Logger::get().debug() << "\tChosen: " << name << "\t" 
			//	<< contigStart << "\t" << bestSegment->start << "\t" << bestSegment->end;
			//contigStart += bestSegment->end - bestSegment->start;

			//theoreticalLen += bestSegment->end - bestSegment->start;

			auto& sequence = (!bestSegment->readSequence) ? 
							  _asmSeqs.getSeq(bestSegment->seqId) :
							  _readSeqs.getSeq(bestSegment->seqId);

			int32_t leftFlank = std::min(Constants::maxSeparation,
										 bestSegment->start);
			if (i == 0) leftFlank = 0;
			int32_t rightFlank = std::min(Constants::maxSeparation,
											(int32_t)sequence.length() - 
											bestSegment->end);
			if (i == contig.path.size() - 1) rightFlank = 0;

			//Logger::get().debug() << "\tLeft Flank " << leftFlank
			//	<< " Right flank: " << rightFlank;

			contigPath.sequences
				.push_back(sequence.substr(bestSegment->start - leftFlank,
										   bestSegment->end - bestSegment->start
										   	+ leftFlank + rightFlank));

			if (i != 0)
			{
				int32_t overlapLen = prevFlank + leftFlank;
				OverlapRange ovlp;
				ovlp.curBegin = prevLength - overlapLen;
				ovlp.curEnd = prevLength;
				ovlp.curLen = prevLength;
				ovlp.extBegin = 0;
				ovlp.extEnd = overlapLen;
				ovlp.extLen = contigPath.sequences.back().length();
				contigPath.overlaps.push_back(ovlp);
			}
			prevFlank = rightFlank;
			prevLength = contigPath.sequences.back().length();
		}
		auto fastaRec = gen.generateLinear(contigPath);
		//Logger::get().debug() << "Expected: " << theoreticalLen 
		//	<< " generated " << fastaRec.sequence.length();
		contig.sequence = fastaRec.sequence.str();
	}
}

void GraphProcessor::dumpRepeats(const std::vector<GraphAlignment>& readAlignments,
								 const std::string& outFile)
{
	std::ofstream fout(outFile);
	if (!fout.is_open()) throw std::runtime_error("Can't open " + outFile);

	for (auto& contig : _contigs)
	{
		if (!contig.path.front()->isRepetitive()) continue;

		bool isSimple = true;
		std::unordered_set<GraphEdge*> inputs;
		for (auto& edge : contig.path.front()->nodeLeft->inEdges)
		{
			inputs.insert(edge);
			if (edge->isRepetitive()) isSimple = false;
		}
		std::unordered_set<GraphEdge*> outputs;
		for (auto& edge : contig.path.back()->nodeRight->outEdges)
		{
			outputs.insert(edge);
			if (edge->isRepetitive()) isSimple = false;
		}
		if (!isSimple || inputs.size() != outputs.size()) continue;

		std::unordered_set<GraphEdge*> innerEdges(contig.path.begin(), 
												  contig.path.end());

		std::unordered_set<FastaRecord::Id> allReads;
		std::unordered_map<GraphEdge*, 
						   std::unordered_set<FastaRecord::Id>> inputEdges;
		std::unordered_map<GraphEdge*, 
						   std::unordered_set<FastaRecord::Id>> outputEdges;


		fout << "#Repeat " << contig.id.signedId() << "\t"
			<< inputs.size() << std::endl;

		//classifying reads into inner, input, output
		for (auto& readAln : readAlignments)
		{
			bool repeatRead = false;
			for (auto& alnEdge : readAln)
			{
				if (innerEdges.count(alnEdge.edge)) repeatRead = true;
			}

			if (!repeatRead) continue;

			allReads.insert(readAln.front().overlap.curId);
			for (auto& alnEdge : readAln)
			{
				for (auto& inputEdge : inputs)
				{
					if (alnEdge.edge == inputEdge) 
					{
						inputEdges[inputEdge]
							.insert(readAln.front().overlap.curId);
					}
				}
				for (auto& outputEdge : outputs)
				{
					if (alnEdge.edge == outputEdge) 
					{
						outputEdges[outputEdge]
							.insert(readAln.front().overlap.curId);
					}
				}
			}
		}

		//dump as text
		fout << "\n#All reads\t" << allReads.size() << std::endl;
		for (auto& readId : allReads)
		{
			fout << _readSeqs.seqName(readId) << std::endl;
		}
		fout << std::endl;

		for (auto& inputEdge : inputs)
		{
			//TODO: more accurate version!
			int ctgId = 0;
			for (auto& ctg : _contigs)
			{
				for (auto& edge : ctg.path) 
				{
					if (edge == inputEdge) ctgId = ctg.id.signedId();
				}
			}

			fout << "#Input " << ctgId << "\t" 
				<< inputEdges[inputEdge].size() << std::endl;

			for (auto& readId : inputEdges[inputEdge])
			{
				fout << _readSeqs.seqName(readId) << std::endl;
			}
			fout << std::endl;
		}

		for (auto& outputEdge : outputs)
		{
			int ctgId = 0;
			for (auto& ctg : _contigs)
			{
				for (auto& edge : ctg.path) 
				{
					if (edge == outputEdge) ctgId = ctg.id.signedId();
				}
			}

			fout << "#Output " << ctgId << "\t" 
				<< outputEdges[outputEdge].size() << std::endl;

			for (auto& readId : outputEdges[outputEdge])
			{
				fout << _readSeqs.seqName(readId) << std::endl;
			}
			fout << std::endl;
		}
	}
}

std::vector<Contig> GraphProcessor::edgesPaths() const
{
	std::vector<Contig> paths;
	for (auto& edge : _graph.iterEdges())
	{
		GraphPath path = {edge};
		paths.emplace_back(path, edge->edgeId, false,
						   edge->length(), edge->meanCoverage);
		paths.back().repetitive = edge->repetitive;
	}
	this->generateContigSequences(paths);
	return paths;
}

void GraphProcessor::outputDot(bool contigs, const std::string& filename)
{
	if (contigs)
	{
		this->outputEdgesDot(_contigs, filename);
	}
	else
	{
		this->outputEdgesDot(this->edgesPaths(), filename);
	}
}

void GraphProcessor::outputGfa(bool contigs, const std::string& filename)
{
	if (contigs)
	{
		this->outputEdgesGfa(_contigs, filename);
	}
	else
	{
		this->outputEdgesGfa(this->edgesPaths(), filename);
	}
}

void GraphProcessor::outputFasta(bool contigs, const std::string& filename)
{
	if (contigs)
	{
		this->outputEdgesFasta(_contigs, filename);
	}
	else
	{
		this->outputEdgesFasta(this->edgesPaths(), filename);
	}
}



void GraphProcessor::outputEdgesFasta(const std::vector<Contig>& paths,
									  const std::string& filename)
{
	static const size_t FASTA_SLICE = 80;

	std::ofstream fout(filename);
	if (!fout.is_open()) throw std::runtime_error("Can't open " + filename);
	
	for (auto& contig : paths)
	{
		if (!contig.id.strand()) continue;

		std::string nameTag = contig.circular ? "circular" : "linear";
		fout << ">" << nameTag << "_" << contig.id.signedId() << std::endl;
		for (size_t c = 0; c < contig.sequence.length(); c += FASTA_SLICE)
		{
			fout << contig.sequence.substr(c, FASTA_SLICE) << std::endl;
		}
	}
}

void GraphProcessor::outputEdgesGfa(const std::vector<Contig>& paths,
							    	const std::string& filename)
{
	std::ofstream fout(filename);
	if (!fout.is_open()) throw std::runtime_error("Can't open " + filename);

	fout << "H\tVN:Z:1.0\n";
	for (auto& contig : paths)
	{
		if (!contig.id.strand()) continue;

		size_t kmerCount = contig.sequence.size() * contig.meanCoverage;
		fout << "S\t" << contig.name() << "\t"<< contig.sequence << "\tKC:i:" <<
			kmerCount << std::endl;
	}

	for (auto& contigLeft : paths)
	{
		for (auto& contigRight : paths)
		{
			if (contigLeft.path.back()->nodeRight != 
				contigRight.path.front()->nodeLeft) continue;

			std::string leftSign = contigLeft.id.strand() ? "+" :"-";
			std::string leftName = contigLeft.nameUnsigned();

			std::string rightSign = contigRight.id.strand() ? "+" :"-";
			std::string rightName = contigRight.nameUnsigned();

			fout << "L\t" << leftName << "\t" << leftSign << "\t" <<
				rightName << "\t" << rightSign << "\t0M\n";
		}
	}
}

void GraphProcessor::outputEdgesDot(const std::vector<Contig>& paths,
									const std::string& filename)
{
	std::ofstream fout(filename);
	if (!fout.is_open()) throw std::runtime_error("Can't open " + filename);

	fout << "digraph {\n";
	fout << "node [shape = circle, label = \"\"]\n";
	
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

	for (auto& contig : paths)
	{
		std::stringstream lengthStr;
		if (contig.length < 5000)
		{
			lengthStr << std::fixed << std::setprecision(1) 
				<< (float)contig.length / 1000 << "k";
		}
		else
		{
			lengthStr << contig.length / 1000 << "k";
		}
		lengthStr << " " << contig.meanCoverage << "x";

		if (contig.repetitive)
		{
			std::string color = idToColor(contig.id);

			fout << "\"" << nodeToId(contig.path.front()->nodeLeft) 
				 << "\" -> \"" << nodeToId(contig.path.back()->nodeRight)
				 << "\" [label = \"id " << contig.id.signedId() << 
				 "\\l" << lengthStr.str() << "\", color = \"" 
				 << color << "\" " << " penwidth = 3] ;\n";
		}
		else
		{
			fout << "\"" << nodeToId(contig.path.front()->nodeLeft)
				 << "\" -> \"" << nodeToId(contig.path.back()->nodeRight)
				 << "\" [label = \"id " << contig.id.signedId()
				 << "\\l" << lengthStr.str() << "\", color = \"black\"] ;\n";
		}
	}

	fout << "}\n";

}
