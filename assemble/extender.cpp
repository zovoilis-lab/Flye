//(c) 2016 by Authors
//This file is a part of ABruijn program.
//Released under the BSD license (see LICENSE file)

#include <limits>
#include <algorithm>
#include <iomanip>
#include <stack>

#include "logger.h"
#include "extender.h"

namespace
{
	template <class T>
	T median(const std::vector<T>& vals)
	{
		if (vals.empty()) return T();
		std::vector<T> tmp = vals;
		std::sort(tmp.begin(), tmp.end());
		//NOTE: there's a bug in libstdc++ nth_element, 
		//that sometimes leads to a segfault
		//std::nth_element(tmp.begin(), tmp.begin() + tmp.size() / 2, 
		//				 tmp.end());
		return tmp[tmp.size() / 2];
	}

	template <class T>
	T mean(const std::vector<T>& vals)
	{
		if (vals.empty()) return T();
		T sum = 0;
		for (const T& val : vals) sum += val;
		return sum / vals.size();
	}
	
	template <class T>
	T vecMax(const std::vector<T>& vals)
	{
		if (vals.empty()) return T();
		T max = 0;
		for (const T& val : vals) max = std::max(max, val);
		return max;
	}
}

ContigPath Extender::extendContig(FastaRecord::Id startRead)
{
	ContigPath contigPath;
	FastaRecord::Id curRead = startRead;
	_chromosomeStart.clear();
	_overlapsStart = true;
	bool rightExtension = true;

	std::unordered_set<FastaRecord::Id> curPathVisited;	//for debugging only

	Logger::get().debug() << "Start Read: " << 
				_seqContainer.seqName(startRead);

	while(true)
	{
		FastaRecord::Id extRead = this->stepRight(curRead);

		if (extRead != FastaRecord::ID_NONE && _chromosomeStart.empty())
		{
			for (auto& ovlp : _ovlpDetector.getOverlapIndex().at(extRead))
			{
				if (this->extendsRight(ovlp.reverse()) && 
					!this->isBranching(ovlp.extId) &&
					!this->isBranching(ovlp.extId.rc())) 
				{
					_chromosomeStart.insert(ovlp.extId);
				}
			}
		}

		if (extRead != FastaRecord::ID_NONE) 
		{
			Logger::get().debug() << "Extension: " << 
				    	_seqContainer.seqName(extRead);
			if (curPathVisited.count(extRead)) 
				Logger::get().debug() << "Visited in this path";
			if (_visitedReads.count(extRead)) 
				Logger::get().debug() << "Visited globally";
		}
		else
		{
			Logger::get().debug() << "No extension found"; 
		}

		if (!_overlapsStart && _chromosomeStart.count(extRead) && 
			rightExtension)
		{
			Logger::get().debug() << "Circular contig";
			contigPath.circular = true;
			contigPath.reads.push_back(extRead);
			break;
		}

		if (_visitedReads.count(extRead) || extRead == FastaRecord::ID_NONE)
		{
			if (rightExtension && !contigPath.reads.empty())
			{
				Logger::get().debug() << "Changing direction";
				rightExtension = false;
				extRead = contigPath.reads.front().rc();
				std::reverse(contigPath.reads.begin(), contigPath.reads.end());
				for (size_t i = 0; i < contigPath.reads.size(); ++i) 
				{
					contigPath.reads[i] = contigPath.reads[i].rc();
				}
				contigPath.reads.pop_back();
			}
			else
			{
				Logger::get().debug() << "Linear contig";
				break;
			}
		}

		_visitedReads.insert(extRead);
		_visitedReads.insert(extRead.rc());
		curPathVisited.insert(extRead);
		curPathVisited.insert(extRead.rc());

		contigPath.reads.push_back(extRead);
		curRead = extRead;
	}

	Logger::get().debug() << "Assembled contig with " << contigPath.reads.size()
						  << " reads";
	return contigPath;
}


void Extender::assembleContigs()
{
	Logger::get().info() << "Extending reads";
	const int MIN_CONTIG = 15;
	const int MIN_EXTENSIONS = std::max(_coverage / 10, 1);

	_visitedReads.clear();
	for (auto& indexPair : _seqContainer.getIndex())
	{	
		if (_visitedReads.count(indexPair.first) ||
			_chimDetector.isChimeric(indexPair.first) ||
			this->countRightExtensions(indexPair.first) < MIN_EXTENSIONS ||
			this->countRightExtensions(indexPair.first.rc()) < MIN_EXTENSIONS ||
			this->isBranching(indexPair.first)) 
		{
			continue;
		}

		ContigPath path = this->extendContig(indexPair.first);

		if (path.reads.size() >= MIN_CONTIG || _contigPaths.empty())
		{
			//marking visited reads
			for (auto& readId : path.reads)
			{
				if (!this->isBranching(readId) &&
					!this->isBranching(readId.rc()))
				{
					for (auto& ovlp : _ovlpDetector.getOverlapIndex().at(readId))
					{
						_visitedReads.insert(ovlp.extId);
						_visitedReads.insert(ovlp.extId.rc());
					}
				}
			}
			
			_contigPaths.push_back(std::move(path));
		}
	}

	Logger::get().info() << "Assembled " << _contigPaths.size() << " contigs";
}


void Extender::coveredReads(const std::unordered_set<FastaRecord::Id>& allReads,
					   		FastaRecord::Id startRead, 
						    std::unordered_set<FastaRecord::Id>& result)
{
	std::unordered_set<FastaRecord::Id> visited;
	std::stack<FastaRecord::Id> dfsStack;
	dfsStack.push(startRead);

	while (!dfsStack.empty())
	{
		FastaRecord::Id curRead = dfsStack.top();
		dfsStack.pop();
		
		if (visited.count(curRead)) continue;

		visited.insert(curRead);
		for (auto& ovlp : _ovlpDetector.getOverlapIndex().at(curRead))
		{
			if (allReads.count(ovlp.extId) && this->coversRight(ovlp.reverse()))
			{
				dfsStack.push(ovlp.extId);
			}
		}
	}

	visited.erase(startRead);
	result = std::move(visited);
}


int Extender::rightMultiplicity(FastaRecord::Id readId)
{
	std::unordered_set<FastaRecord::Id> extensions;
	auto& overlaps = _ovlpDetector.getOverlapIndex().at(readId);
	for (auto& ovlp : overlaps)
	{
		if (this->extendsRight(ovlp)) extensions.insert(ovlp.extId);
	}
	if (extensions.size() < 2) return 1;

	//Logger::get().debug() << "Multiplicity of " << _seqContainer.seqName(readId);
	std::unordered_map<FastaRecord::Id, 
					   std::unordered_set<FastaRecord::Id>> coveredByNode;
	std::unordered_map<FastaRecord::Id, int> globalCoverage;
	std::unordered_set<FastaRecord::Id> clusterIds;

	for (auto& extCandidate : extensions)
	{
		this->coveredReads(extensions, extCandidate, 
						   coveredByNode[extCandidate]);
	}

	while(true)
	{
		FastaRecord::Id maxUniqueCoveredId = FastaRecord::ID_NONE;
		int maxUniqueCovered = 0;
		for (auto& extCandidate : extensions)
		{
			//ignore reads that already belong to clusters
			if (globalCoverage.count(extCandidate)) continue;

			//count number of reads the candidate extends
			int extUniqueCovered = 0;
			for (auto read : coveredByNode[extCandidate])
			{
				if (!globalCoverage.count(read)) ++extUniqueCovered;
			}

			//keeping maximum
			if (extUniqueCovered > maxUniqueCovered)
			{
				maxUniqueCovered = extUniqueCovered;
				maxUniqueCoveredId = extCandidate;
			}
		}

		if (maxUniqueCoveredId == FastaRecord::ID_NONE) break;

		//Logger::get().debug() << "\tCl: " << _seqContainer.seqName(maxUniqueCoveredId);
		clusterIds.insert(maxUniqueCoveredId);
		for (auto readId : coveredByNode[maxUniqueCoveredId]) 
		{
			//Logger::get().debug() << "\t\t " << _seqContainer.seqName(readId);
			++globalCoverage[readId];
		}
	}

	//counting cluster sizes
	int numClusters = 0;
	std::string strClusters;
	for (auto& clustId : clusterIds)
	{
		int clustSize = 0;
		for (auto readId : coveredByNode[clustId]) 
		{
			if (globalCoverage[readId] < 2)
			{
				++clustSize;
			}
		}
		++numClusters;
		strClusters += std::to_string(clustSize) + " ";
	}
	//Logger::get().debug() << "\tClusters: " << strClusters;

	return std::max(1, numClusters);
}


//makes one extension to the right
FastaRecord::Id Extender::stepRight(FastaRecord::Id readId)
{
	auto& overlaps = _ovlpDetector.getOverlapIndex().at(readId);
	std::unordered_set<FastaRecord::Id> extensions;

	bool locOverlapsStart = false;
	for (auto& ovlp : overlaps)
	{
		if (this->extendsRight(ovlp)) 
		{
			if (_chromosomeStart.count(ovlp.extId))
			{
				locOverlapsStart = true;
				if (!_overlapsStart) return ovlp.extId;	//circular chromosome
			}

			if (this->countRightExtensions(ovlp.extId) > 0)
			{
				extensions.insert(ovlp.extId);
			}
		}
	}
	if (!_chromosomeStart.empty() && _overlapsStart && !locOverlapsStart) 
	{
		_overlapsStart = false;
	}
	
	//rank extension candidates
	std::unordered_map<FastaRecord::Id, 
					   std::tuple<int, int, int, int, int>> supportIndex;

	for (auto& extCandidate : extensions)
	{
		int leftSupport = 0;
		int rightSupport = 0;
		int ovlpSize = 0;
		for (auto& ovlp : _ovlpDetector.getOverlapIndex().at(extCandidate))
		{
			if (ovlp.extId == readId) ovlpSize = ovlp.curRange();
			if (!extensions.count(ovlp.extId)) continue;

			if (this->extendsRight(ovlp)) ++rightSupport;
			if (this->extendsRight(ovlp.reverse())) ++leftSupport;
		}

		int minSupport = std::min(leftSupport, rightSupport);

		int resolvesRepeat = (!this->isBranching(readId) && 
							  this->majorClusterAgreement(extCandidate.rc(), 
														  readId.rc())) ||
							 (!this->isBranching(extCandidate.rc()) &&
							   this->majorClusterAgreement(readId, 
							   							   extCandidate));
		int stepAhead = this->resolvableRepeat(extCandidate);

		supportIndex[extCandidate] = std::make_tuple(resolvesRepeat, stepAhead, minSupport,
													 rightSupport, ovlpSize);

		Logger::get().debug() << "\t" 
				    << _seqContainer.seqName(extCandidate) << "\trr:" << resolvesRepeat
					<< "\tstep:" << stepAhead << "\tcons:(" << leftSupport << "," 
					<< rightSupport << ")\tovlp:" << ovlpSize;
	}
	

	auto bestSupport = std::make_tuple(0, 0, 0, 0, 0);
	auto bestExtension = FastaRecord::ID_NONE;
	for (auto& extCandidate : extensions)
	{
		if (supportIndex[extCandidate] > bestSupport)
		{
			bestSupport = supportIndex[extCandidate];
			bestExtension = extCandidate;
		}
	}

	if (bestExtension != FastaRecord::ID_NONE && std::get<0>(bestSupport) == 0)
	{
		Logger::get().debug() << "Can't resolve repeat";
		return FastaRecord::ID_NONE;
	}

	return bestExtension;
}

bool Extender::resolvableRepeat(FastaRecord::Id readId)
{
	if (this->isBranching(readId) && this->isBranching(readId.rc())) 
	{
		return false;
	}

	for (auto& ovlp : _ovlpDetector.getOverlapIndex().at(readId))
	{
		if (this->extendsRight(ovlp) &&
			this->countRightExtensions(ovlp.extId) > 0)
		{
			if ((!this->isBranching(readId) && 
				 this->majorClusterAgreement(ovlp.extId.rc(), readId.rc())) ||
				(!this->isBranching(ovlp.extId.rc()) &&
				 this->majorClusterAgreement(readId, ovlp.extId)))
			{
				return true;
			}
		}
	}
	return false;
}

bool Extender::majorClusterAgreement(FastaRecord::Id leftRead,
									 FastaRecord::Id rightRead)
{
	auto& overlaps = _ovlpDetector.getOverlapIndex().at(rightRead.rc());
	std::unordered_set<FastaRecord::Id> extensions;
	for (auto& ovlp : overlaps)
	{
		if (this->extendsRight(ovlp)) extensions.insert(ovlp.extId);
	}

	size_t maxCovered = 0;
	std::unordered_set<FastaRecord::Id> coveredSet;
	for (auto& extCandidate : extensions)
	{
		std::unordered_set<FastaRecord::Id> resCovered;
		this->coveredReads(extensions, extCandidate, resCovered);
		if (resCovered.size() > maxCovered)
		{
			maxCovered = resCovered.size();
			coveredSet = resCovered;
		}
	}

	return coveredSet.empty() || coveredSet.count(leftRead.rc());
}

int Extender::countRightExtensions(FastaRecord::Id readId)
{
	int count = 0;
	for (auto& ovlp : _ovlpDetector.getOverlapIndex().at(readId))
	{
		if (this->extendsRight(ovlp)) ++count;
	}
	return count;
}

bool Extender::isBranching(FastaRecord::Id readId)
{
	if (!_readsMultiplicity.count(readId))
	{
		_readsMultiplicity[readId] = this->rightMultiplicity(readId);
	}
	return _readsMultiplicity[readId] > 1;
}

bool Extender::extendsRight(const OverlapRange& ovlp)
{
	return !_chimDetector.isChimeric(ovlp.extId) && 
		   ovlp.rightShift > _maximumJump;
}

bool Extender::coversRight(const OverlapRange& ovlp)
{
	return ovlp.rightShift > 0;
}
