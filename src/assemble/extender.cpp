//(c) 2016 by Authors
//This file is a part of ABruijn program.
//Released under the BSD license (see LICENSE file)

#include <limits>
#include <algorithm>
#include <iomanip>
#include <stack>

#include "../common/config.h"
#include "../common/logger.h"
#include "../common/parallel.h"
#include "extender.h"


Extender::ExtensionInfo Extender::extendContig(FastaRecord::Id startRead)
{

	std::unordered_set<FastaRecord::Id> currentReads;
	currentReads.insert(startRead);
	currentReads.insert(startRead.rc());

	bool rightExtension = true;
	FastaRecord::Id currentRead = startRead;
	std::vector<int> numOverlaps;
	ExtensionInfo exInfo;
	exInfo.reads.push_back(startRead);

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
		auto overlaps = _ovlpContainer.lazySeqOverlaps(currentRead);
		std::vector<OverlapRange> extensions;
		int innerOverlaps = 0;
		for (auto& ovlp : overlaps)
		{
			if (_innerReads.contains(ovlp.extId)) ++innerOverlaps;
			if (this->extendsRight(ovlp)) extensions.push_back(ovlp);
		}
		numOverlaps.push_back(extensions.size());

		//sort from longes to shortest overlap
		std::sort(extensions.begin(), extensions.end(), 
				  [](const OverlapRange& a, const OverlapRange& b)
					 {return a.curRange() > b.curRange();});

		//checking if read overlaps with one of the used reads
		bool foundExtension = false;
		bool overlapsVisited = false;

		//getting extension
		int minExtensions = (int)extensions.size() / 
							Constants::maxCoverageDropRate;
		FastaRecord::Id bestSuspicious = FastaRecord::ID_NONE;
		if (!overlapsVisited)
		{
			for (auto& ovlp : extensions)
			{
				if(leftExtendsStart(ovlp.extId)) continue;

				//try to find a good one
				if (!_chimDetector.isChimeric(ovlp.extId) &&
					!this->isRightRepeat(ovlp.extId) &&
					this->countRightExtensions(ovlp.extId) > minExtensions)
				{
					foundExtension = true;
					currentRead = ovlp.extId;
					break;
				}
				//or not so good
				else if(bestSuspicious == FastaRecord::ID_NONE)
				{
					bestSuspicious = ovlp.extId;
				}
			}

			if (!foundExtension && bestSuspicious != FastaRecord::ID_NONE)
			{
				++exInfo.numSuspicious;
				foundExtension = true;
				currentRead = bestSuspicious;
			}
		}

		overlapsVisited |= _innerReads.contains(currentRead);
		overlapsVisited |= currentReads.count(currentRead);
		if (foundExtension) 
		{
			//Logger::get().debug() << "Extension: " << 
			//	    	_readsContainer.seqName(currentRead) << " " << innerOverlaps
			//			<< " " << extensions.size();
			exInfo.reads.push_back(currentRead);
		}
		else
		{
			rightExtension ? exInfo.leftTip = true : exInfo.rightTip = true;
		}

		if (!foundExtension || overlapsVisited)
		{
			//right extension done, try to extend left from start read
			if (rightExtension && !exInfo.reads.empty())
			{
				exInfo.stepsToTurn = exInfo.reads.size();
				rightExtension = false;
				currentRead = exInfo.reads.front().rc();
				std::reverse(exInfo.reads.begin(), exInfo.reads.end());
				for (size_t i = 0; i < exInfo.reads.size(); ++i) 
				{
					exInfo.reads[i] = exInfo.reads[i].rc();
				}
			}
			//done with extension
			else
			{
				break;
			}
		}

		currentReads.insert(currentRead);
		currentReads.insert(currentRead.rc());
	}

	int64_t meanOvlps = 0;
	for (int num : numOverlaps) meanOvlps += num;
	exInfo.meanOverlaps = meanOvlps / numOverlaps.size();

	return exInfo;
}


void Extender::assembleContigs()
{
	Logger::get().info() << "Extending reads";
	_chimDetector.estimateGlobalCoverage();
	_innerReads.clear();
	cuckoohash_map<FastaRecord::Id, size_t> coveredReads;
	
	std::mutex indexMutex;
	auto processRead = [this, &indexMutex, &coveredReads] 
		(FastaRecord::Id startRead)
	{
		if (coveredReads.contains(startRead)) return true;

		int numInnerOvlp = 0;
		for (auto& ovlp : _ovlpContainer.lazySeqOverlaps(startRead))
		{
			if (_innerReads.contains(ovlp.extId)) ++numInnerOvlp;
		}
		if (numInnerOvlp > 0) return true;

		if (_chimDetector.isChimeric(startRead) ||
			this->isRightRepeat(startRead) ||
			this->isRightRepeat(startRead.rc())) return true;
		
		//Good to go!
		ExtensionInfo exInfo = this->extendContig(startRead);
		if (exInfo.reads.size() < Constants::minReadsInContig) return true;

		//Exclusive part - updating the overall assembly
		std::lock_guard<std::mutex> guard(indexMutex);
		
		int innerCount = 0;
		for (auto& readId : exInfo.reads)
		{
			if (_innerReads.contains(readId)) ++innerCount;
		}
		if (innerCount > Constants::maxInnerFraction)
		{
			Logger::get().debug() << "Discarded contig with "
				<< exInfo.reads.size() << " reads and "
				<< innerCount << " inner overlaps";
			return false;
		}

		Logger::get().debug() << "Assembled contig" 
			<< "\n\tWith " << exInfo.reads.size() << " reads"
			<< "\n\tStart read: " << _readsContainer.seqName(startRead)
			<< "\n\tAt position: " << exInfo.stepsToTurn
			<< "\n\tleftTip: " << exInfo.leftTip 
			<< " rightTip: " << exInfo.rightTip
			<< "\n\tSuspicios: " << exInfo.numSuspicious
			<< "\n\tMean extensions: " << exInfo.meanOverlaps
			<< "\n\tInner reads: " << innerCount;
		
		//update inner read index
		std::unordered_set<FastaRecord::Id> rightExtended;
		std::unordered_set<FastaRecord::Id> leftExtended;
		for (auto& readId : exInfo.reads)
		{
			//repetitive read - will bring to much "off-target" reads
			if (this->isRightRepeat(readId) ||
				this->isRightRepeat(readId.rc())) continue;

			//so each read is covered from the left and right
			for (auto& ovlp : _ovlpContainer.lazySeqOverlaps(readId))
			{
				coveredReads.insert(ovlp.extId, true);

				//contained inside the other read - probably repetitive
				if (ovlp.leftShift > Constants::maximumJump &&
					ovlp.rightShift < -Constants::maximumJump) continue;

				if (ovlp.leftShift > Constants::maximumJump)
				{
					leftExtended.insert(ovlp.extId);
					rightExtended.insert(ovlp.extId.rc());
				}
				if (ovlp.rightShift < -Constants::maximumJump)	
				{
					rightExtended.insert(ovlp.extId);
					leftExtended.insert(ovlp.extId.rc());
				}
			}
		}
		for (auto& read : rightExtended)
		{
			if (leftExtended.count(read)) _innerReads.insert(read, true);
		}

		Logger::get().debug() << "Inner: " << 
			_innerReads.size() << " covered: " << coveredReads.size()
			<< " total: "<< _readsContainer.getIndex().size();
		
		_readLists.push_back(std::move(exInfo));
		return true;
	};

	std::function<void(const FastaRecord::Id&)> threadWorker = 
		[processRead] (const FastaRecord::Id& readId)
	{
		processRead(readId);
	};
	std::vector<FastaRecord::Id> allReads;
	for (auto& hashPair : _readsContainer.getIndex())
	{
		allReads.push_back(hashPair.first);
	}
	processInParallel(allReads, threadWorker,
					  Parameters::get().numThreads, true);

	this->convertToContigs();
	Logger::get().info() << "Assembled " << _contigPaths.size()
		<< " draft contigs";
}

void Extender::convertToContigs()
{
	for (auto& exInfo : _readLists)
	{
		ContigPath path;
		path.name = "contig_" + std::to_string(_contigPaths.size() + 1);

		for (size_t i = 0; i < exInfo.reads.size() - 1; ++i)
		{
			bool found = false;
			OverlapRange readsOvlp;

			for (auto& ovlp : _ovlpContainer.lazySeqOverlaps(exInfo.reads[i]))
			{
				if (ovlp.extId == exInfo.reads[i + 1]) 
				{
					readsOvlp = ovlp;
					found = true;
					break;
				}
			}
			if (!found) throw std::runtime_error("Ovlp not found!");

			path.sequences.push_back(_readsContainer.getSeq(exInfo.reads[i]));
			path.overlaps.push_back(readsOvlp);
		}
		path.sequences.push_back(_readsContainer.getSeq(exInfo.reads.back()));
		_contigPaths.push_back(std::move(path));
	}
}

int Extender::countRightExtensions(FastaRecord::Id readId) const
{
	int count = 0;
	for (auto& ovlp : _ovlpContainer.lazySeqOverlaps(readId))
	{
		if (this->extendsRight(ovlp)) ++count;
	}
	return count;
}

bool Extender::isRightRepeat(FastaRecord::Id readId) const
{
	int maxExtensions = _chimDetector.getOverlapCoverage() * 10;
	return this->countRightExtensions(readId) >= maxExtensions;
}

bool Extender::extendsRight(const OverlapRange& ovlp) const
{
	return ovlp.rightShift > Constants::maximumJump;
}
