//(c) 2016-2017 by Authors
//This file is a part of ABruijn program.
//Released under the BSD license (see LICENSE file)

#include "contig_extender.h"

void ContigExtender::generateUnbranchingPaths()
{
	GraphProcessor proc(_graph, _asmSeqs, _readSeqs);
	_unbranchingPaths = proc.getUnbranchingPaths();
	//this->generateContigSequences(_contigs);

	for (auto& ctg : _unbranchingPaths)
	{
		if (ctg.id.strand())
		{
			Logger::get().debug() << "UPath " << ctg.id.signedId() 
							<< ": " << ctg.edgesStr();
		}
	}

	//this->generateContigSequences(_contigs);
	Logger::get().debug() << "Final graph contiain " 
		<< _unbranchingPaths.size() / 2 << " egdes";
}

/*
void ContigExtender::extendContigs(std::vector<UnbranchingPath>& contigs)
{
	std::unordered_set<GraphEdge*> coveredRepeats;
	std::unordered_map<GraphEdge*, bool> repeatDirections;
	std::unordered_set<GraphEdge*> tandems;
	auto readAln = _aligner.getAlignments();

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
}*/

