//(c) 2016-2017 by Authors
//This file is a part of ABruijn program.
//Released under the BSD license (see LICENSE file)

#include "output_generator.h"
#include "../sequence/contig_generator.h"
#include <iomanip>

#undef NDEBUG
#include <cassert>

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

	this->generateContigSequences(_contigs);
	Logger::get().debug() << "Final graph contiain " << _contigs.size() / 2 
		<< " egdes";
}

void OutputGenerator::extendContigs(const std::vector<GraphAlignment>& readAln,
									const std::string& outFile)
{
	std::unordered_set<GraphEdge*> coveredRepeats;
	std::unordered_map<GraphEdge*, bool> repeatDirections;
	std::unordered_set<GraphEdge*> tandems;

	auto extendPath = [&readAln, this, &coveredRepeats, 
					   &repeatDirections, &tandems]
		(UnbranchingPath& upath, bool direction)
	{
		bool extendFwd = !upath.path.back()->nodeRight->outEdges.empty();
		if (!extendFwd) return std::string();

		//first, choose the longest aligned read from this edge
		int32_t maxExtension = 0;
		GraphAlignment bestAlignment;
		for (auto& path : readAln)
		{
			for (size_t i = 0; i < path.size(); ++i)
			{
				if (path[i].edge == upath.path.back() &&
					i < path.size() - 1)
				{
					size_t j = i + 1;
					while (j < path.size() && 
						   path[j].edge->repetitive) ++j;
					if (i == j) continue;

					int32_t alnLen = path[j - 1].overlap.curEnd - 
									 path[i + 1].overlap.curBegin;
					if (alnLen > maxExtension)
					{
						maxExtension = alnLen;
						bestAlignment.clear();
						std::copy(path.begin() + i + 1, path.begin() + j, 
								  std::back_inserter(bestAlignment));
					}
					break;
				}
			}
		}
		if (maxExtension == 0) return std::string();

		//check if some repeats are traversed multiple by a single read
		std::unordered_map<GraphEdge*, int> multiplicity;
		for (auto& aln : bestAlignment)
		{
			++multiplicity[aln.edge];
		}
		for (auto& edgeMult : multiplicity)
		{
			if (edgeMult.second > 1) 
			{
				tandems.insert(edgeMult.first);
				tandems.insert(_graph.complementEdge(edgeMult.first));
			}
		}

		//Logger::get().debug() << "Ctg " << upath.name() 
		//	<< " extension " << maxExtension;

		//check if we are not traverse non-tandem repeats in prohibited 
		//directions
		bool canExtend = true;
		for (auto& aln : bestAlignment)
		{
			++multiplicity[aln.edge];
			if (repeatDirections.count(aln.edge))
			{
				if (!tandems.count(aln.edge) &&
					repeatDirections.at(aln.edge) != direction)
				{
					canExtend = false;
					break;
				}
			}
		}
		if (!canExtend) return std::string();

		//if we traverse some repeats for the first time, record the direction
		for (auto& aln : bestAlignment)
		{
			//Logger::get().debug() << "\t" << aln.edge->edgeId.signedId();
			if (aln.overlap.extLen - aln.overlap.extEnd < 
					Constants::maxSeparation)
			{
				repeatDirections[aln.edge] = direction;
				repeatDirections[_graph.complementEdge(aln.edge)] = !direction;

				coveredRepeats.insert(aln.edge);
				coveredRepeats.insert(_graph.complementEdge(aln.edge));
			}
		}

		FastaRecord::Id seqId = bestAlignment.front().overlap.curId;
		int32_t start = bestAlignment.front().overlap.curBegin;
		int32_t end = bestAlignment.back().overlap.curEnd;

		auto extSeq = _readSeqs.getSeq(seqId).substr(start, end - start);
		return direction ? extSeq.str() : extSeq.complement().str();
	};

	std::vector<UnbranchingPath> extendedContigs;
	std::unordered_map<FastaRecord::Id, UnbranchingPath*> idToPath;
	for (auto& ctg : _contigs)
	{
		idToPath[ctg.id] = &ctg;
	}
	for (auto& ctg : _contigs)
	{
		if (ctg.repetitive || !ctg.id.strand()) continue;
		if (!idToPath.count(ctg.id.rc())) continue;	//self-complement

		std::string rightExt = extendPath(ctg, true);
		std::string leftExt = extendPath(*idToPath[ctg.id.rc()], false);
		extendedContigs.push_back(ctg);
		extendedContigs.back().sequence = leftExt + 
							extendedContigs.back().sequence + rightExt;
	}

	int numCovered = 0;
	for (auto& ctg : _contigs)
	{
		if (!ctg.repetitive || !ctg.id.strand()) continue;

		bool covered = true;
		for (auto& edge : ctg.path)
		{
			if (!coveredRepeats.count(edge)) covered = false;
		}
		if (!covered)
		{
			extendedContigs.push_back(ctg);
		}
		else
		{
			++numCovered;
			Logger::get().debug() << "Covered: " << ctg.name();
		}
	}
	Logger::get().debug() << "Covered " << numCovered << " repetitive contigs";

	Logger::get().info() << "Generated " << extendedContigs.size() 
		<< " extended contigs";
	this->outputEdgesFasta(extendedContigs, outFile);
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

			auto& sequence = (!bestSegment->readSequence) ? 
							  _asmSeqs.getSeq(bestSegment->seqId) :
							  _readSeqs.getSeq(bestSegment->seqId);

			int32_t leftFlank = std::min(5000,
										 bestSegment->start);
			if (i == 0) 
			{
				leftFlank = 0;
			}
			int32_t rightFlank = std::min(5000,
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
		if (!isSimple || inputs.size() != outputs.size() ||
			inputs.empty()) continue;

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
	std::vector<FastaRecord> records;
	for (auto& contig : paths)
	{
		if (!contig.id.strand()) continue;

		std::string nameTag = contig.circular ? "circular" : "linear";
		nameTag += "_" + std::to_string(contig.id.signedId());
		records.emplace_back(DnaSequence(contig.sequence), nameTag, 
							 FastaRecord::ID_NONE);
	}
	SequenceContainer::writeFasta(records, filename);
}

void OutputGenerator::outputEdgesGfa(const std::vector<UnbranchingPath>& paths,
							    	 const std::string& filename)
{
	Logger::get().debug() << "Writing Gfa";
	FILE* fout = fopen(filename.c_str(), "w");
	if (!fout) throw std::runtime_error("Can't open " + filename);

	fprintf(fout, "H\tVN:Z:1.0\n");
	for (auto& contig : paths)
	{
		if (!contig.id.strand()) continue;

		size_t kmerCount = contig.sequence.size() * contig.meanCoverage;
		fprintf(fout, "S\t%s\t%s\tKC:i:%d\n", contig.name().c_str(), 
				contig.sequence.c_str(), (int)kmerCount);
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

void OutputGenerator::outputEdgesDot(const std::vector<UnbranchingPath>& paths,
									 const std::string& filename)
{
	Logger::get().debug() << "Writing Dot";

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
