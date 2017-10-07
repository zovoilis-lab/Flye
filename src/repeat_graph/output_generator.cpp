//(c) 2016-2017 by Authors
//This file is a part of ABruijn program.
//Released under the BSD license (see LICENSE file)

#include "output_generator.h"
#include "../sequence/contig_generator.h"
#include <iomanip>

void OutputGenerator::generateContigs()
{
	GraphProcessor proc(_graph, _asmSeqs, _readSeqs);
	_contigs = proc.getUnbranchingPaths();
	this->generateContigSequences(_contigs);

	for (auto& ctg : _contigs)
	{
		if (ctg.id.strand())
		{
			Logger::get().debug() << "Contig " << ctg.id.signedId() 
							<< ": " << ctg.edgesStr();
		}
	}

	Logger::get().info() << "Generated " << _contigs.size() / 2 << " contigs";
}

void OutputGenerator::
	generateContigSequences(std::vector<UnbranchingPath>& contigs) const
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
			if (bestSegment->length() == 0) continue;

			//auto name = (!bestSegment->readSequence) ? 
			//				_asmSeqs.seqName(bestSegment->seqId) :
			//				_readSeqs.seqName(bestSegment->seqId);
			//Logger::get().debug() << "\tChosen: " << name << "\t" 
			//	<< contigStart << "\t" << bestSegment->start << "\t" << bestSegment->end;
			//contigStart += bestSegment->end - bestSegment->start;

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

			int32_t substrLen = bestSegment->end - bestSegment->start
										   	+ leftFlank + rightFlank;
			contigPath.sequences
				.push_back(sequence.substr(bestSegment->start - leftFlank,
										   substrLen));

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
		contig.sequence = fastaRec.sequence.str();
	}
}

void OutputGenerator::dumpRepeats(const std::vector<GraphAlignment>& readAlignments,
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

std::vector<UnbranchingPath> OutputGenerator::edgesPaths() const
{
	std::vector<UnbranchingPath> paths;
	for (auto& edge : _graph.iterEdges())
	{
		GraphPath path = {edge};
		paths.emplace_back(path, edge->edgeId, false,
						   edge->length(), edge->meanCoverage);
		paths.back().repetitive = edge->repetitive;
	}
	//this->generateContigSequences(paths);
	return paths;
}

void OutputGenerator::outputDot(bool contigs, const std::string& filename)
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

void OutputGenerator::outputGfa(bool contigs, const std::string& filename)
{
	if (contigs)
	{
		this->outputEdgesGfa(_contigs, filename);
	}
	else
	{
		auto paths = this->edgesPaths();
		this->generateContigSequences(paths);
		this->outputEdgesGfa(paths, filename);
	}
}

void OutputGenerator::outputFasta(bool contigs, const std::string& filename)
{
	if (contigs)
	{
		this->outputEdgesFasta(_contigs, filename);
	}
	else
	{
		auto paths = this->edgesPaths();
		this->generateContigSequences(paths);
		this->outputEdgesFasta(paths, filename);
	}
}

void OutputGenerator::outputEdgesFasta(const std::vector<UnbranchingPath>& paths,
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

void OutputGenerator::outputEdgesGfa(const std::vector<UnbranchingPath>& paths,
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

void OutputGenerator::outputEdgesDot(const std::vector<UnbranchingPath>& paths,
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
