//(c) 2016 by Authors
//This file is a part of ABruijn program.
//Released under the BSD license (see LICENSE file)

#include <limits>
#include <algorithm>
#include <iomanip>
#include <stack>

#include "../common/config.h"
#include "../common/logger.h"
#include "extender.h"


ContigPath Extender::extendContig(FastaRecord::Id startRead)
{
	Logger::get().debug() << "Start Read: " << 
				_readsContainer.seqName(startRead);

	_usedReads.insert(startRead);
	_usedReads.insert(startRead.rc());
	_rightExtension = true;
	FastaRecord::Id currentRead = startRead;
	std::vector<int> numOverlaps;
	ContigPath contigPath;
	contigPath.reads.push_back(startRead);

	auto leftExtendsStart = [startRead, this](const FastaRecord::Id readId)
	{
		for (auto& ovlp : _ovlpContainer.lazySeqOverlaps(readId))
		{
			if (ovlp.extId == startRead &&
				ovlp.leftShift > Constants::maximumJump) 
			{
				return true;
			}
		}
		return false;
	};

	while(true)
	{
		auto& overlaps = _ovlpContainer.lazySeqOverlaps(currentRead);
		std::vector<OverlapRange> extensions;
		for (auto& ovlp : overlaps)
		{
			if (this->extendsRight(ovlp)) 
			{
				extensions.push_back(ovlp);
			}
		}
		numOverlaps.push_back(overlaps.size());

		//getting mean shift, sorting extensions according to it
		int64_t sum = 0;
		for (auto& ovlp : extensions) sum += ovlp.rightShift;
		int32_t meanShift = !extensions.empty() ? sum / extensions.size() : 0;
		std::sort(extensions.begin(), extensions.end(), 
				  [meanShift](const OverlapRange& a, const OverlapRange& b)
					 {return abs(a.rightShift - meanShift) < 
							 abs(b.rightShift - meanShift);});

		//checking if read overlaps with one of the used reads
		bool foundExtension = false;
		bool overlapsVisited = false;
		for (auto& ovlp : extensions)
		{
			if (_usedReads.count(ovlp.extId)) 
			{
				overlapsVisited = true;
				foundExtension = true;
				currentRead = ovlp.extId;
				break;
			}
		}

		//getting extension
		int minExtensions = (int)extensions.size() / 
							Constants::maxCoverageDropRate;
		//Logger::get().debug() << extensions.size();
		if (!overlapsVisited)
		{
			for (auto& ovlp : extensions)
			{
				//Logger::get().debug() << "\t" << this->countRightExtensions(ovlp.extId);
				if (!_chimDetector.isChimeric(ovlp.extId) &&
					this->countRightExtensions(ovlp.extId) > minExtensions &&
					!leftExtendsStart(ovlp.extId))
				{
					foundExtension = true;
					currentRead = ovlp.extId;
					_assembledSequence += ovlp.rightShift;
					break;
				}
			}
		}

		overlapsVisited |= _coveredReads.count(currentRead);
		if (foundExtension) 
		{
			Logger::get().debug() << "Extension: " << 
				    	_readsContainer.seqName(currentRead);

			contigPath.reads.push_back(currentRead);
			_progress.setValue(_assembledSequence);

			if (overlapsVisited)
			{
				Logger::get().debug() << "Already visited"; 
			}
		}
		else
		{
			Logger::get().debug() << "No extension found"; 
		}

		if (!foundExtension || overlapsVisited)
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
			}
			else
			{
				break;
			}
		}

		_usedReads.insert(currentRead);
		_usedReads.insert(currentRead.rc());
	}

	int64_t meanOvlps = 0;
	for (int num : numOverlaps) meanOvlps += num;
	Logger::get().debug() << "Mean overlaps: " << meanOvlps / numOverlaps.size();
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

	_coveredReads.clear();
	_usedReads.clear();

	int numFails = 0;
	for (auto& indexPair : _readsContainer.getIndex())
	{
		if (numFails++ > Constants::extensionTries) break;
		if (_coveredReads.count(indexPair.first) ||
			_chimDetector.isChimeric(indexPair.first) ||
			this->countRightExtensions(indexPair.first) < 
				Constants::minExtensions ||
			this->countRightExtensions(indexPair.first.rc()) < 
				Constants::minExtensions)
		{
			continue;
		}

		ContigPath path = this->extendContig(indexPair.first);
		
		std::unordered_set<FastaRecord::Id> rightExtended;
		std::unordered_set<FastaRecord::Id> leftExtended;
		for (auto& readId : path.reads)
		{
			//so each read is covered by at least two others, from left and right
			for (auto& ovlp : _ovlpContainer.lazySeqOverlaps(readId))
			{
				if (ovlp.leftShift > Constants::maximumJump)
				{
					leftExtended.insert(ovlp.extId);
					rightExtended.insert(ovlp.extId.rc());
				}
				else if (ovlp.rightShift < -Constants::maximumJump)	
				{
					rightExtended.insert(ovlp.extId);
					leftExtended.insert(ovlp.extId.rc());
				}
			}
		}
		for (auto& read : rightExtended)
		{
			if (leftExtended.count(read)) _coveredReads.insert(read);
		}
		///
		
		if (path.reads.size() >= Constants::minReadsInContig)
		{
			numFails = 0;
			Logger::get().debug() << "Assembled contig with " 
				<< path.reads.size() << " reads";
			_contigPaths.push_back(std::move(path));
		}
	}

	_progress.setDone();
	Logger::get().info() << "Assembled " << _contigPaths.size() << " contigs";
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
