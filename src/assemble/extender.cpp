//(c) 2016 by Authors
//This file is a part of ABruijn program.
//Released under the BSD license (see LICENSE file)

#include <limits>
#include <algorithm>
#include <iomanip>
#include <stack>

#include "config.h"
#include "logger.h"
#include "extender.h"


ContigPath Extender::extendContig(FastaRecord::Id startRead)
{
	ContigPath contigPath;
	contigPath.reads.push_back(startRead);
	_rightExtension = true;

	Logger::get().debug() << "Start Read: " << 
				_readsContainer.seqName(startRead);

	FastaRecord::Id currentRead = startRead;
	std::vector<int> numOverlaps;
	while(true)
	{
		auto overlaps = _ovlpContainer.lazySeqOverlaps(currentRead);
		std::vector<OverlapRange> extensions;

		//std::vector<int> extensionShifts;
		for (auto& ovlp : overlaps)
		{
			if (this->extendsRight(ovlp)) 
			{
				extensions.push_back(ovlp);
				//extensionShifts.push_back(ovlp.rightShift);
			}
		}
		numOverlaps.push_back(overlaps.size());

		int64_t sum = 0;
		for (auto& ovlp : extensions) sum += ovlp.rightShift;
		int32_t meanShift = !extensions.empty() ? sum / extensions.size() : 0;

		std::sort(extensions.begin(), extensions.end(), 
				  [meanShift](const OverlapRange& a, const OverlapRange& b)
					 {return abs(a.rightShift - meanShift) < 
							 abs(b.rightShift - meanShift);});

		bool foundExtension = false;
		int numVisited = 0;
		for (auto& ovlp : extensions)
		{
			if (_visitedReads.count(ovlp.extId)) ++numVisited;
		}
		//if (numVisited <= (int)extensions.size() / 2)
		//{
			for (auto& ovlp : extensions)
			{
				if (!_chimDetector.isChimeric(ovlp.extId) &&
					this->countRightExtensions(ovlp.extId) > 0)
				{
					foundExtension = true;
					currentRead = ovlp.extId;
					_assembledSequence += ovlp.rightShift;
					break;
				}
			}
		//}
		Logger::get().debug() << extensions.size() << " " << numVisited;

		if (foundExtension) 
		{
			Logger::get().debug() << "Extension: " << 
				    	_readsContainer.seqName(currentRead);

			contigPath.reads.push_back(currentRead);
			_progress.setValue(_assembledSequence);

			if (_visitedReads.count(currentRead))
			{
				Logger::get().debug() << "Already visited"; 
			}
		}
		else
		{
			Logger::get().debug() << "No extension found"; 
		}

		if (!foundExtension || _visitedReads.count(currentRead))
		{
			if (_rightExtension && !contigPath.reads.empty())
			{
				Logger::get().debug() << "Changing direction";
				_rightExtension = false;
				currentRead = contigPath.reads.front().rc();
				std::reverse(contigPath.reads.begin(), contigPath.reads.end());
				for (size_t i = 0; i < contigPath.reads.size(); ++i) 
				{
					contigPath.reads[i] = contigPath.reads[i].rc();
				}
				//contigPath.reads.pop_back();
			}
			else
			{
				Logger::get().debug() << "Linear contig";
				break;
			}
		}

		_visitedReads.insert(currentRead);
		_visitedReads.insert(currentRead.rc());
	}

	int64_t meanOvlps = 0;
	for (int num : numOverlaps) meanOvlps += num;
	Logger::get().debug() << "Mean overlaps: " << meanOvlps / numOverlaps.size();
	Logger::get().debug() << "Assembled contig with " << contigPath.reads.size()
						  << " reads";
	return contigPath;
}


void Extender::assembleContigs()
{
	Logger::get().info() << "Extending reads";

	/*
	uint64_t lenSum = 0;
	for (auto indexPair : _readsContainer.getIndex()) 
	{
		lenSum += _readsContainer.seqLen(indexPair.first);
	}
	_minimumShift = lenSum / _readsContainer.getIndex().size() 
						/ Constants::shiftToReadLen;*/

	//const int MIN_EXTENSIONS = std::max(_coverage / 
	//									Constants::minExtensionsRate, 1);

	_visitedReads.clear();
	int extended = 0;
	for (auto& indexPair : _readsContainer.getIndex())
	{
		if (extended++ > Constants::extensionTries) break;
		if (_visitedReads.count(indexPair.first)) continue;

		ContigPath path;

		if (_chimDetector.isChimeric(indexPair.first))
		{
			if (!_visitedReads.count(indexPair.first))
			{
				Logger::get().debug() << _chimDetector.isChimeric(indexPair.first) << " "
					<< this->countRightExtensions(indexPair.first) << " "
					<< this->countRightExtensions(indexPair.first.rc());
			}

			path.reads.push_back(indexPair.first);
		}
		else
		{
			path = this->extendContig(indexPair.first);
		}

		for (auto& readId : path.reads)
		{
			for (auto& ovlp : _ovlpContainer.lazySeqOverlaps(readId))
			{
				_visitedReads.insert(ovlp.extId);
				_visitedReads.insert(ovlp.extId.rc());
			}
		}
		
		if (path.reads.size() >= Constants::minReadsInContig)
		{
			_contigPaths.push_back(std::move(path));
		}
	}

	_progress.setDone();
	Logger::get().info() << "Assembled " << _contigPaths.size() << " contigs";
}


//makes one extension to the right
FastaRecord::Id Extender::stepRight(FastaRecord::Id readId)
{
	auto overlaps = _ovlpContainer.lazySeqOverlaps(readId);
	std::vector<OverlapRange> extensions;
	//Logger::get().debug() << "Ovlps: " << overlaps.size();
	//Logger::get().debug() << "Index: " << _ovlpContainer.getOverlapIndex().size();

	//bool locOverlapsStart = false;
	//std::vector<int> extensionShifts;
	for (auto& ovlp : overlaps)
	{
		if (this->extendsRight(ovlp)) 
		{
			/*
			if (_chromosomeStart.count(ovlp.extId))
			{
				locOverlapsStart = true;
				Logger::get().debug() << "Bumped into start";
				//circular chromosome
				if (!_overlapsStart && _rightExtension) return ovlp.extId;	
			}*/
			extensions.push_back(ovlp);

			/*if (this->countRightExtensions(ovlp.extId) > 0)
			{
				extensions.insert(ovlp.extId);
				extensionShifts.push_back(ovlp.rightShift);
			}*/
		}
	}

	int64_t sum = 0;
	for (auto& ovlp : extensions) sum += ovlp.rightShift;
	int32_t meanShift = !extensions.empty() ? sum / extensions.size() : 0;

	std::sort(extensions.begin(), extensions.end(), 
			  [meanShift](const OverlapRange& a, const OverlapRange& b)
			  	 {return abs(a.rightShift - meanShift) < 
				 		 abs(b.rightShift - meanShift);});

	for (auto& ovlp : extensions)
	{
		if (!_chimDetector.isChimeric(ovlp.extId) &&
			this->countRightExtensions(ovlp.extId) > 0) return ovlp.extId;
	}

	return FastaRecord::ID_NONE;
	/*
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
	*/
}

int Extender::countRightExtensions(FastaRecord::Id readId)
{
	int count = 0;
	for (auto& ovlp : _ovlpContainer.lazySeqOverlaps(readId))
	{
		if (this->extendsRight(ovlp)) ++count;
	}
	return count;
}

bool Extender::extendsRight(const OverlapRange& ovlp)
{
	return ovlp.rightShift > Constants::maximumJump;
}
