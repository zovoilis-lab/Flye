//(c) 2016 by Authors
//This file is a part of ABruijn program.
//Released under the BSD license (see LICENSE file)

#include <limits>
#include <algorithm>
#include <iomanip>
#include <stack>

#include "logger.h"
#include "extender.h"
#include "config.h"


namespace
{
	int robustStd(const std::vector<int>& shifts)
	{
		if (shifts.size() < 2) return 0;

		auto _mean = [](const std::vector<int>& vals)
		{
			int sum = 0;
			for(int v : vals) sum += v;
			return (float)sum / vals.size();
		};

		auto _std = [_mean](const std::vector<int>& vals)
		{
			float mean = _mean(vals);
			float sumStd = 0.0f;
			for (float v : vals) sumStd += (v - mean) * (v - mean);
			return sqrt(sumStd / vals.size());  
		};

		float mean = _mean(shifts);
		float std = _std(shifts);

		//filter outliers
		std::vector<int> inliers;
		for (float s : shifts)
		{
			if (fabsf(s - mean) < 2.0f * std) inliers.push_back(s);
		}

		return _std(inliers);
	}
}

ContigPath Extender::extendContig(FastaRecord::Id startRead)
{
	ContigPath contigPath;
	FastaRecord::Id curRead = startRead;
	_chromosomeStart.clear();
	_overlapsStart = true;
	_rightExtension = true;

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
			if (_visitedReads.count(extRead)) 
				Logger::get().debug() << "Visited globally";
		}
		else
		{
			Logger::get().debug() << "No extension found"; 
		}

		if (!_overlapsStart && _chromosomeStart.count(extRead) && 
			_rightExtension)
		{
			Logger::get().debug() << "Circular contig";
			contigPath.circular = true;
			contigPath.reads.push_back(extRead);
			break;
		}

		if (_visitedReads.count(extRead) || extRead == FastaRecord::ID_NONE)
		{
			if (_rightExtension && !contigPath.reads.empty())
			{
				Logger::get().debug() << "Changing direction";
				_rightExtension = false;
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

	uint64_t lenSum = 0;
	for (auto indexPair : _seqContainer.getIndex()) 
	{
		if (this->countRightExtensions(indexPair.first) > 0)
		{
			lenSum += _seqContainer.seqLen(indexPair.first);
		}
	}
	_minimumShift = lenSum / _seqContainer.getIndex().size() 
						/ Constants::shiftToReadLen;

	//left for debugging
	//for (auto& indexPair : _seqContainer.getIndex())
	//{
	//	_readsMultiplicity[indexPair.first] = 
	//				this->rightMultiplicity(indexPair.first);
	//}
	//

	const int MIN_EXTENSIONS = std::max(_coverage / 
										Constants::minExtensionsRate, 1);
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

		if (path.reads.size() >= Constants::minReadsInContig)
		{
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
	if (extensions.size() > 2UL * _coverage) return 2;

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
		if (clusterIds.empty())
		{
			_maxClusters[readId] = coveredByNode[maxUniqueCoveredId];
			_maxClusters[readId].insert(readId);
		}
		clusterIds.insert(maxUniqueCoveredId);
		if (clusterIds.size() > 1) return 2;

		for (auto readId : coveredByNode[maxUniqueCoveredId]) 
		{
			//Logger::get().debug() << "\t\t " << _seqContainer.seqName(readId);
			++globalCoverage[readId];
		}
	}
	//Logger::get().debug() << "\tClusters: " << strClusters;

	return std::max(1UL, clusterIds.size());
}


//makes one extension to the right
FastaRecord::Id Extender::stepRight(FastaRecord::Id readId)
{
	auto& overlaps = _ovlpDetector.getOverlapIndex().at(readId);
	std::unordered_set<FastaRecord::Id> extensions;

	bool locOverlapsStart = false;
	std::vector<int> extensionShifts;
	for (auto& ovlp : overlaps)
	{
		if (this->extendsRight(ovlp)) 
		{
			if (_chromosomeStart.count(ovlp.extId))
			{
				locOverlapsStart = true;
				Logger::get().debug() << "Bumped into start";
				//circular chromosome
				if (!_overlapsStart && _rightExtension) return ovlp.extId;	
			}

			if (this->countRightExtensions(ovlp.extId) > 0)
			{
				extensions.insert(ovlp.extId);
				extensionShifts.push_back(ovlp.rightShift);
			}
		}
	}
	if (!_chromosomeStart.empty() && _overlapsStart && !locOverlapsStart) 
	{
		_overlapsStart = false;
	}
	
	Logger::get().debug() << "Shift: " << robustStd(extensionShifts);
	//rank extension candidates
	std::unordered_map<FastaRecord::Id, 
					   std::tuple<int, int, int, int, int>> supportIndex;

	int totalSupport = 0;
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

		totalSupport += leftSupport + rightSupport;
		int minSupport = std::min(leftSupport, rightSupport);
		int resolvesRepeat = this->resolvesRepeat(readId, extCandidate);
		int stepAhead = this->stepAhead(extCandidate);

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

	if (bestExtension != FastaRecord::ID_NONE && 
		robustStd(extensionShifts) < _minimumShift && 
		extensions.size() > (size_t)_coverage / Constants::minExtensionsRate)
	{
		Logger::get().debug() << "End of linear chromosome";
		return FastaRecord::ID_NONE;
	}

	return bestExtension;
}

bool Extender::stepAhead(FastaRecord::Id readId)
{
	if (this->isBranching(readId) && this->isBranching(readId.rc())) 
	{
		return false;
	}

	std::unordered_set<FastaRecord::Id> candidates;
	for (auto& ovlp : _ovlpDetector.getOverlapIndex().at(readId))
	{
		if (this->extendsRight(ovlp) &&
			this->countRightExtensions(ovlp.extId) > 0)
		{
			candidates.insert(ovlp.extId);
		}
	}

	int goodReads = 0;
	for (auto candidate : candidates)
	{
		if (this->resolvesRepeat(readId, candidate))
		{
			++goodReads;
		}
		if ((float)goodReads / candidates.size() > Constants::minGoodReads)
		{
			return true;
		}
	}

	return false;
}


bool Extender::resolvesRepeat(FastaRecord::Id leftRead, 
							  FastaRecord::Id rightRead)
{
	return (!this->isBranching(leftRead) && 
			this->majorClusterAgreement(leftRead, rightRead)) ||
			(!this->isBranching(rightRead.rc()) &&
			this->majorClusterAgreement(rightRead.rc(), leftRead.rc()));
}

bool Extender::majorClusterAgreement(FastaRecord::Id leftRead,
									 FastaRecord::Id rightRead)
{
	if (!_readsMultiplicity.count(leftRead))
	{
		_readsMultiplicity[leftRead] = 
				this->rightMultiplicity(leftRead);
	}

	return _maxClusters[leftRead].empty() || 
		   _maxClusters[leftRead].count(rightRead);
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
		   ovlp.rightShift > Constants::maxumumJump;
}

bool Extender::coversRight(const OverlapRange& ovlp)
{
	return ovlp.rightShift > 0;
}
