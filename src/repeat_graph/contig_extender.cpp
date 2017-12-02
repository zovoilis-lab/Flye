//(c) 2016-2017 by Authors
//This file is a part of ABruijn program.
//Released under the BSD license (see LICENSE file)

#include "contig_extender.h"
#include "output_generator.h"

void ContigExtender::generateUnbranchingPaths()
{
	GraphProcessor proc(_graph, _asmSeqs, _readSeqs);
	_unbranchingPaths = proc.getUnbranchingPaths();

	_edgeToPath.clear();
	for (auto& path : _unbranchingPaths)
	{
		if (path.id.strand())
		{
			Logger::get().debug() << "UPath " << path.id.signedId() 
							<< ": " << path.edgesStr();
		}

		for (auto& edge : path.path)
		{
			_edgeToPath[edge] = &path;
		}
	}

	Logger::get().debug() << "Final graph contiain " 
		<< _unbranchingPaths.size() / 2 << " egdes";
}


void ContigExtender::generateContigs()
{
	OutputGenerator outGen(_graph, _aligner, _asmSeqs, _readSeqs);
	auto coreSeqs = outGen.generatePathSequences(_unbranchingPaths);
	std::unordered_map<UnbranchingPath*, FastaRecord*> upathsSeqs;
	for (size_t i = 0; i < _unbranchingPaths.size(); ++i)
	{
		upathsSeqs[&_unbranchingPaths[i]] = &coreSeqs[i];
	}

	std::vector<const GraphAlignment*> interestingAlignments;
	for (auto& aln : _aligner.getAlignments())
	{
		if (aln.size() > 1) interestingAlignments.push_back(&aln);
	}

	std::unordered_set<GraphEdge*> coveredRepeats;
	std::unordered_map<const GraphEdge*, bool> repeatDirections;
	auto canTraverse = [&repeatDirections] (const GraphEdge* edge)
	{
		return edge->selfComplement || 
			   !repeatDirections.count(edge) || 
			   repeatDirections.at(edge);
	};

	typedef std::pair<GraphPath, std::string> PathAndSeq;
	auto extendPathRight =
		[this, &coveredRepeats, &repeatDirections, &upathsSeqs, 
		 &canTraverse, &interestingAlignments] (UnbranchingPath& upath)
	{

		bool extendFwd = !upath.path.back()->nodeRight->outEdges.empty();
		if (!extendFwd) return PathAndSeq();

		//first, choose the longest aligned read from this edge
		int32_t maxExtension = 0;
		GraphAlignment bestAlignment;
		for (auto pathPtr : interestingAlignments)
		{
			const GraphAlignment& path = *pathPtr;
			for (size_t i = 0; i < path.size(); ++i)
			{
				if (path[i].edge == upath.path.back() &&
					i < path.size() - 1)
				{
					size_t j = i + 1;
					while (j < path.size() && path[j].edge->repetitive &&
						   canTraverse(path[j].edge)) ++j;
					if (j == i + 1) break;

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
		if (maxExtension == 0) return PathAndSeq();

		//record repeat directions
		for (auto& aln : bestAlignment)
		{
			repeatDirections[aln.edge] = true;
			repeatDirections[_graph.complementEdge(aln.edge)] = false;
			coveredRepeats.insert(aln.edge);
			coveredRepeats.insert(_graph.complementEdge(aln.edge));
		}

		GraphPath extendedPath;
		for (auto& aln : bestAlignment) extendedPath.push_back(aln.edge);

		//generate extension sequence
		std::string extendedSeq;
		if (!bestAlignment.empty())
		{
			auto lastUpath = this->asUPaths({bestAlignment.back().edge}).front();
			int32_t overhang = bestAlignment.back().overlap.extLen - 
							   bestAlignment.back().overlap.extEnd;

			Logger::get().debug() << "Ctg " << upath.id.signedId() <<
				" overhang " << overhang << " upath " << lastUpath->id.signedId();
			if (overhang > Constants::maxSeparation && !lastUpath->isLoop())
			{
				bestAlignment.pop_back();
			}
			if (!bestAlignment.empty())
			{
				FastaRecord::Id readId = bestAlignment.front().overlap.curId;
				int32_t readStart = bestAlignment.front().overlap.curBegin;
				int32_t readEnd = bestAlignment.back().overlap.curEnd;
				extendedSeq = _readSeqs.getSeq(readId)
					.substr(readStart, readEnd - readStart).str();
			}
			if (overhang > Constants::maxSeparation && !lastUpath->isLoop())
			{
				extendedSeq += upathsSeqs[lastUpath]->sequence.str();
			}
		}
		
		return PathAndSeq(extendedPath, extendedSeq);
	};

	std::unordered_map<FastaRecord::Id, UnbranchingPath*> idToPath;
	for (auto& ctg : _unbranchingPaths)
	{
		idToPath[ctg.id] = &ctg;
	}

	for (auto& upath : _unbranchingPaths)
	{
		if (upath.repetitive || !upath.id.strand()) continue;
		if (!idToPath.count(upath.id.rc())) continue;	//self-complement

		auto rightExt = extendPathRight(upath);
		auto leftExt = extendPathRight(*idToPath[upath.id.rc()]);
		leftExt.first = _graph.complementPath(leftExt.first);
		leftExt.second = DnaSequence(leftExt.second).complement().str();

		Contig contig(upath);
		auto leftPaths = this->asUPaths(leftExt.first);
		auto rightPaths = this->asUPaths(rightExt.first);
		for (auto& path : leftPaths)
		{
			for (auto& edge : path->path)
			{
				contig.graphEdges.path.insert(contig.graphEdges.path.end() - 1, 
											  edge);
			}
			contig.graphPaths.insert(contig.graphPaths.end() - 1, path);
		}
		for (auto& path : rightPaths)
		{
			for (auto& edge : path->path)
			{
				contig.graphEdges.path.push_back(edge);
			}
			contig.graphPaths.push_back(path);
		}

		auto coreSeq = upathsSeqs[&upath]->sequence.str();
		contig.sequence = DnaSequence(leftExt.second + coreSeq + rightExt.second);

		_contigs.push_back(std::move(contig));
	}

	//add repetitive contigs that were not covered by the extended paths
	int numCovered = 0;
	for (auto& upath : _unbranchingPaths)
	{
		if (!upath.repetitive || !upath.id.strand()) continue;

		bool covered = false;
		for (auto& edge : upath.path)
		{
			if (coveredRepeats.count(edge)) covered = true;
		}
		if (!covered)
		{
			_contigs.emplace_back(upath);
			_contigs.back().sequence = upathsSeqs[&upath]->sequence;
		}
		else
		{
			++numCovered;
			Logger::get().debug() << "Covered: " << upath.id.signedId();
		}
	}

	Logger::get().debug() << "Covered " << numCovered << " repetitive contigs";
	Logger::get().info() << "Generated " << _contigs.size() << " contigs";
}

std::vector<UnbranchingPath*> ContigExtender::asUPaths(const GraphPath& path)
{
	std::vector<UnbranchingPath*> upathRepr;
	for (size_t i = 0; i < path.size(); ++i)
	{
		UnbranchingPath* upath = _edgeToPath.at(path[i]);
		if (upathRepr.empty() || upathRepr.back() != upath ||
			path[i - 1] == path[i])
		{
			upathRepr.push_back(upath);
		}
	}

	return upathRepr;
}

void ContigExtender::outputStatsTable(const std::string& filename)
{
	std::ofstream fout(filename);
	if (!fout) throw std::runtime_error("Can't write " + filename);

	fout << "contig_id\tlength\tcoverage\tcircular\trepeat"
		<< "\tmult\ttelomere\tgraph_path\n";

	char YES_NO[] = {'N', 'Y'};

	for (auto& ctg : _contigs)
	{
		std::string pathStr;
		for (auto& upath : ctg.graphPaths)
		{
			pathStr += std::to_string(upath->id.signedId()) + ",";
		}
		pathStr.pop_back();

		int minMult = std::numeric_limits<int>::max();
		for (auto& edge : ctg.graphEdges.path) 
		{
			if (edge->multiplicity > 0) 
			{
				minMult = std::min(minMult, edge->multiplicity);
			}
		}

		fout << ctg.graphEdges.name() << "\t" << ctg.sequence.length() << "\t" 
			<< ctg.graphEdges.meanCoverage << "\t"
			<< YES_NO[ctg.graphEdges.circular]
			<< "\t" << YES_NO[ctg.graphEdges.repetitive] << "\t"
			<< minMult << "\t" << "L:" 
			<< YES_NO[ctg.graphEdges.path.front()->nodeLeft->isTelomere()]
			<< ",R:" 
			<< YES_NO[ctg.graphEdges.path.back()->nodeRight->isTelomere()]
			<< "\t" << pathStr << "\n";

		Logger::get().debug() << "Contig: " << ctg.graphEdges.id.signedId()
			<< ": " << pathStr;
	}
}

void ContigExtender::outputContigs(const std::string& filename)
{
	std::vector<FastaRecord> contigsFasta;
	for (auto& ctg : _contigs)
	{
		contigsFasta.emplace_back(ctg.sequence, ctg.graphEdges.name(), 
					   			  FastaRecord::ID_NONE);
	}
	SequenceContainer::writeFasta(contigsFasta, filename);
}
