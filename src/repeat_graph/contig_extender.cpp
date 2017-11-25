//(c) 2016-2017 by Authors
//This file is a part of ABruijn program.
//Released under the BSD license (see LICENSE file)

#include "contig_extender.h"

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
	std::unordered_set<GraphEdge*> coveredRepeats;
	std::unordered_map<GraphEdge*, bool> repeatDirections;
	//std::unordered_set<GraphEdge*> tandems;

	auto extendPathRight =
		[this, &coveredRepeats, &repeatDirections] (UnbranchingPath& upath)
	{
		auto& readAln = _aligner.getAlignments();

		bool extendFwd = !upath.path.back()->nodeRight->outEdges.empty();
		if (!extendFwd) return GraphPath();

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
		if (maxExtension == 0) return GraphPath();

		//check if some repeats are traversed multiple times by a single read
		/*std::unordered_map<GraphEdge*, int> multiplicity;
		for (auto& aln : bestAlignment) ++multiplicity[aln.edge];
		for (auto& edgeMult : multiplicity)
		{
			if (edgeMult.second > 1) 
			{
				tandems.insert(edgeMult.first);
				tandems.insert(_graph.complementEdge(edgeMult.first));
			}
		}*/

		//check if we are not traverse non-tandem repeats in prohibited directions
		bool canExtend = true;
		for (auto& aln : bestAlignment)
		{
			if (repeatDirections.count(aln.edge))
			{
				//if (!tandems.count(aln.edge) &&
				if (!aln.edge->selfComplement && !repeatDirections.at(aln.edge))
				{
					canExtend = false;
					break;
				}
			}
		}
		if (!canExtend) return GraphPath();

		//if we traverse some repeats for the first time, record the direction
		for (auto& aln : bestAlignment)
		{
			repeatDirections[aln.edge] = true;
			repeatDirections[_graph.complementEdge(aln.edge)] = false;

			coveredRepeats.insert(aln.edge);
			coveredRepeats.insert(_graph.complementEdge(aln.edge));
		}

		GraphPath extendedPath;
		for (auto& aln : bestAlignment) extendedPath.push_back(aln.edge);
		
		return extendedPath;
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
		leftExt = _graph.complementPath(leftExt);

		Contig contig(upath);
		auto leftPaths = this->asUPaths(leftExt);
		auto rightPaths = this->asUPaths(rightExt);
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
		_contigs.push_back(std::move(contig));
	}

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
		}
		else
		{
			++numCovered;
			Logger::get().debug() << "Covered: " << upath.id.signedId();
		}
	}

	Logger::get().debug() << "Covered " << numCovered << " repetitive contigs";
	Logger::get().info() << "Generated " << _contigs.size() << " contigs";

	for (auto& ctg : _contigs)
	{
		std::string pathStr;
		for (auto& upath : ctg.graphPaths)
		{
			pathStr += std::to_string(upath->id.signedId()) + " -> ";
		}
		pathStr.erase(pathStr.size() - 4);

		Logger::get().debug() << "Contig: " << ctg.graphEdges.id.signedId()
			<< ": " << pathStr;
	}
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


std::vector<UnbranchingPath> ContigExtender::getContigPaths()
{
	std::vector<UnbranchingPath> upaths;
	for (auto& ctg : _contigs) upaths.push_back(ctg.graphEdges);
	return upaths;
}
