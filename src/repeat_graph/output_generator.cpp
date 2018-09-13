//(c) 2016-2017 by Authors
//This file is a part of ABruijn program.
//Released under the BSD license (see LICENSE file)

#include "output_generator.h"
#include "../sequence/consensus_generator.h"
#include <iomanip>


//Generates FASTA from the given graph paths
std::vector<FastaRecord> OutputGenerator::
	generatePathSequences(const std::vector<UnbranchingPath>& paths) const
{
	std::vector<ContigPath> contigParts;

	for (auto& contig : paths)
	{
		//As each edge might correspond to multiple sequences,
		//we need to select them so as to minimize the
		//number of original contigs (that were used to build the graph)
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
		contigPath.name = contig.name();
		int32_t prevFlank = 0;
		int32_t prevSubLength = 0;

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
			}
			if (bestSegment->length() == 0) continue;

			DnaSequence sequence;
			switch (bestSegment->segType)
			{
				case SequenceSegment::Asm:
					sequence = _asmSeqs.getSeq(bestSegment->seqId);
					break;
				case SequenceSegment::Read:
					sequence = _readSeqs.getSeq(bestSegment->seqId);
					break;
			}

			//make the consecutive sequences overlapping if possible,
			//so the consensus module can correct possibly imprecise
			//ends of the edges sequences
			int32_t leftFlank = std::min((int32_t)Config::get("max_separation"),
										 bestSegment->start);
			if (i == 0) 
			{
				leftFlank = 0;
			}
			int32_t rightFlank = std::min((int32_t)Config::get("max_separation"),
										  (int32_t)sequence.length() - 
										  		bestSegment->end);
			if (i == contig.path.size() - 1) 
			{
				rightFlank = 0;
			}

			int32_t curSubLength = bestSegment->length() + leftFlank + rightFlank;
			contigPath.sequences
				.push_back(sequence.substr(bestSegment->start - leftFlank,
										   curSubLength));

			assert(leftFlank >= 0);
			assert(rightFlank >= 0);

			if (i != 0)
			{
				int32_t adjustedNextFlank = 
					std::min(leftFlank, prevSubLength - prevFlank);
				int32_t adjustedPrevFlank = 
					std::min(prevFlank, curSubLength - leftFlank);
				int32_t overlapLen = adjustedPrevFlank + adjustedNextFlank;
				OverlapRange ovlp;
				ovlp.curBegin = prevSubLength - overlapLen;
				ovlp.curEnd = prevSubLength;
				ovlp.curLen = prevSubLength;
				ovlp.extBegin = 0;
				ovlp.extEnd = overlapLen;
				ovlp.extLen = curSubLength;
				contigPath.overlaps.push_back(ovlp);

				assert(ovlp.curBegin >= 0);
			}
			prevFlank = rightFlank;
			prevSubLength = curSubLength;
		}
		contigParts.push_back(contigPath);
	}

	//finally, generate a consensus
	ConsensusGenerator gen;
	return gen.generateConsensuses(contigParts, /*verbose*/ false);
}

//dumps repeat information for the consecutive analysis
void OutputGenerator::dumpRepeats(const std::vector<UnbranchingPath>& paths,
								  const std::string& outFile)
{
	std::ofstream fout(outFile);
	if (!fout.is_open()) throw std::runtime_error("Can't open " + outFile);

	for (auto& contig : paths)
	{
		if (!contig.path.front()->isRepetitive() ||
			contig.path.front()->selfComplement) continue;

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
		if (!isSimple || inputs.size() != outputs.size() ||
			inputs.size() < 2) continue;

		std::unordered_set<GraphEdge*> innerEdges(contig.path.begin(), 
												  contig.path.end());

		std::unordered_set<FastaRecord::Id> allReads;
		std::unordered_map<GraphEdge*, 
						   std::unordered_set<FastaRecord::Id>> inputEdges;
		std::unordered_map<GraphEdge*, 
						   std::unordered_set<FastaRecord::Id>> outputEdges;


		fout << "#Repeat " << contig.id.signedId() << "\t"
			<< inputs.size() << "\n";

		//classifying reads into inner, input, output
		for (auto& readAln : _aligner.getAlignments())
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
		fout << "\n#All reads\t" << allReads.size() << "\n";
		for (auto& readId : allReads)
		{
			fout << _readSeqs.seqName(readId) << "\n";
		}
		fout << "\n";

		for (auto& inputEdge : inputs)
		{
			//TODO: more accurate version!
			int ctgId = 0;
			for (auto& ctg : paths)
			{
				for (auto& edge : ctg.path) 
				{
					if (edge == inputEdge) ctgId = ctg.id.signedId();
				}
			}

			fout << "#Input " << ctgId << "\t" 
				<< inputEdges[inputEdge].size() << "\n";

			for (auto& readId : inputEdges[inputEdge])
			{
				fout << _readSeqs.seqName(readId) << "\n";
			}
			fout << "\n";
		}

		for (auto& outputEdge : outputs)
		{
			int ctgId = 0;
			for (auto& ctg : paths)
			{
				for (auto& edge : ctg.path) 
				{
					if (edge == outputEdge) ctgId = ctg.id.signedId();
				}
			}

			fout << "#Output " << ctgId << "\t" 
				<< outputEdges[outputEdge].size() << "\n";

			for (auto& readId : outputEdges[outputEdge])
			{
				fout << _readSeqs.seqName(readId) << "\n";
			}
			fout << "\n";
		}
	}
}

void OutputGenerator::outputFasta(const std::vector<UnbranchingPath>& paths,
								  const std::string& filename)
{
	std::vector<UnbranchingPath> posStrandPaths;
	for (auto& path : paths)
	{
		if (path.id.strand()) posStrandPaths.push_back(path);
	}
	SequenceContainer::writeFasta(this->generatePathSequences(posStrandPaths), 
								  filename);
}

void OutputGenerator::outputGfa(const std::vector<UnbranchingPath>& paths,
							    const std::string& filename)
{
	auto sequences = this->generatePathSequences(paths);

	Logger::get().debug() << "Writing Gfa";
	FILE* fout = fopen(filename.c_str(), "w");
	if (!fout) throw std::runtime_error("Can't open " + filename);

	fprintf(fout, "H\tVN:Z:1.0\n");
	for (size_t i = 0; i < paths.size(); ++i)
	{
		if (!paths[i].id.strand()) continue;

		//size_t kmerCount = sequences[i].sequence.length() * paths[i].meanCoverage;
		fprintf(fout, "S\t%s\t%s\tDP:i:%d\n", paths[i].name().c_str(), 
				sequences[i].sequence.str().c_str(), (int)paths[i].meanCoverage);
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

			fprintf(fout, "L\t%s\t%s\t%s\t%s\t0M\n", leftName.c_str(), 
					leftSign.c_str(), rightName.c_str(), rightSign.c_str());
		}
	}
}

void OutputGenerator::outputDot(const std::vector<UnbranchingPath>& paths,
								const std::string& filename)
{
	Logger::get().debug() << "Writing Dot";

	std::ofstream fout(filename);
	if (!fout.is_open()) throw std::runtime_error("Can't open " + filename);

	fout << "digraph {\n";
	fout << "nodesep = 0.5;\n";
	fout << "node [shape = circle, label = \"\", height = 0.3];\n";
	
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

	for (auto& node : _graph.iterNodes())
	{
		if (node->isTelomere())
		{
			fout << "\"" << nodeToId(node) 
				<< "\" [style = \"filled\", fillcolor = \"grey\"];\n";
		}
	}

	//coloring repeat clusters
	const std::string COLORS[] = {"red", "darkgreen", "blue", "goldenrod", 
								  "cadetblue1", "darkorchid", "aquamarine1", 
								  "darkgoldenrod1", "deepskyblue1", 
								  "darkolivegreen3"};
	std::vector<GraphEdge*> dfsStack;
	std::unordered_set<GraphEdge*> visited;
	std::unordered_map<GraphEdge*, std::string> edgeColors;
	size_t colorId = 0;
	for (auto& edge: _graph.iterEdges())
	{
		if (!edge->isRepetitive() || visited.count(edge)) continue;
		dfsStack.push_back(edge);
		while(!dfsStack.empty())
		{
			auto curEdge = dfsStack.back(); 
			dfsStack.pop_back();
			if (visited.count(curEdge)) continue;
			edgeColors[curEdge] = COLORS[colorId];
			edgeColors[_graph.complementEdge(curEdge)] = COLORS[colorId];
			visited.insert(curEdge);
			visited.insert(_graph.complementEdge(curEdge));
			for (auto adjEdge: curEdge->adjacentEdges())
			{
				if (adjEdge->isRepetitive() && !visited.count(adjEdge))
				{
					dfsStack.push_back(adjEdge);
				}
			}
		}
		colorId = (colorId + 1) % (sizeof(COLORS) / sizeof(COLORS[0]));
	}

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
			std::string color = edgeColors[contig.path.front()];
			std::string direction = contig.path.front()->selfComplement ?
									", dir = both" : "";
			fout << "\"" << nodeToId(contig.path.front()->nodeLeft) 
				 << "\" -> \"" << nodeToId(contig.path.back()->nodeRight)
				 << "\" [label = \"id " << contig.id.signedId() << 
				 "\\l" << lengthStr.str() << "\", color = \"" 
				 << color << "\" " << ", penwidth = 3" << direction << "] ;\n";
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
