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


Extender::ReadsList Extender::extendContig(FastaRecord::Id startRead)
{
	Logger::get().debug() << "Start Read: " << 
				_readsContainer.seqName(startRead);

	std::unordered_set<FastaRecord::Id> currentReads;
	currentReads.insert(startRead);
	currentReads.insert(startRead.rc());

	_rightExtension = true;
	FastaRecord::Id currentRead = startRead;
	std::vector<int> numOverlaps;
	ReadsList contigPath;
	contigPath.push_back(startRead);

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
		int innerOverlaps = 0;
		for (auto& ovlp : overlaps)
		{
			if (_innerReads.count(ovlp.extId)) ++innerOverlaps;
			if (this->extendsRight(ovlp)) extensions.push_back(ovlp);
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
		bool mayStop = _innerReads.empty() || 
			innerOverlaps > (int)extensions.size() / Constants::maxCoverageDropRate;
		for (auto& ovlp : extensions)
		{
			if (mayStop && (currentReads.count(ovlp.extId) ||
				_innerReads.count(ovlp.extId))) 
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
		
		FastaRecord::Id bestSuspicious = FastaRecord::ID_NONE;

		if (!overlapsVisited)
		{
			for (auto& ovlp : extensions)
			{
				//Logger::get().debug() << "\t" << this->countRightExtensions(ovlp.extId);
				if(leftExtendsStart(ovlp.extId)) continue;

				if (!_chimDetector.isChimeric(ovlp.extId) &&
					this->countRightExtensions(ovlp.extId) > minExtensions)
				{
					foundExtension = true;
					currentRead = ovlp.extId;
					_assembledSequence += ovlp.rightShift;
					break;
				}
				else if(bestSuspicious == FastaRecord::ID_NONE)
				{
					bestSuspicious = ovlp.extId;
				}
			}

			if (!foundExtension && bestSuspicious != FastaRecord::ID_NONE)
			{
				Logger::get().debug() << "Suspicious!";
				foundExtension = true;
				currentRead = bestSuspicious;
			}
		}

		overlapsVisited |= _innerReads.count(currentRead);
		overlapsVisited |= currentReads.count(currentRead);
		if (foundExtension) 
		{
			Logger::get().debug() << "Extension: " << 
				    	_readsContainer.seqName(currentRead) << " " << innerOverlaps
						<< " " << extensions.size();

			contigPath.push_back(currentRead);
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
			if (_rightExtension && !contigPath.empty())
			{
				Logger::get().debug() << "Changing direction";
				_rightExtension = false;
				currentRead = contigPath.front().rc();
				std::reverse(contigPath.begin(), contigPath.end());
				for (size_t i = 0; i < contigPath.size(); ++i) 
				{
					contigPath[i] = contigPath[i].rc();
				}
			}
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
	Logger::get().debug() << "Mean overlaps: " << meanOvlps / numOverlaps.size();
	return contigPath;
}


void Extender::assembleContigs()
{
	Logger::get().info() << "Extending reads";
	_chimDetector.estimateGlobalCoverage();
	_innerReads.clear();
	std::unordered_set<FastaRecord::Id> coveredReads;
	//_coveredReads.clear();

	int numChecked = 0;
	int readsToCheck = Constants::startReadsPercent * 
							_readsContainer.getIndex().size();
	for (auto& indexPair : _readsContainer.getIndex())
	{
		if (numChecked++ > readsToCheck) break;
		/*Logger::get().debug() << _readsContainer.seqName(indexPair.first) << "\t"
			<< _coveredReads.count(indexPair.first) << "\t"
			<< _chimDetector.isChimeric(indexPair.first) << "\t"
			<< this->countRightExtensions(indexPair.first) << "\t"
			<< this->countRightExtensions(indexPair.first.rc());*/

		if (coveredReads.count(indexPair.first) ||
			_chimDetector.isChimeric(indexPair.first) ||
			this->countRightExtensions(indexPair.first) < 
				Constants::minExtensions ||
			this->countRightExtensions(indexPair.first.rc()) < 
				Constants::minExtensions)
		{
			continue;
		}
		numChecked = 0;

		auto& overlaps = _ovlpContainer.lazySeqOverlaps(indexPair.first);
		int numOvlp = 0;
		for (auto& ovlp : overlaps)
		{
			if (_innerReads.count(ovlp.extId)) ++numOvlp;
		}
		//float usedIndex = (float)numOvlp / overlaps.size();
		//Logger::get().debug() << "Used: " << usedIndex;
		if (numOvlp > 0) continue;

		ReadsList path = this->extendContig(indexPair.first);
		
		std::unordered_set<FastaRecord::Id> rightExtended;
		std::unordered_set<FastaRecord::Id> leftExtended;
		for (auto& readId : path)
		{
			//so each read is covered from the left and right
			for (auto& ovlp : _ovlpContainer.lazySeqOverlaps(readId))
			{
				coveredReads.insert(ovlp.extId);
				coveredReads.insert(ovlp.extId.rc());
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
			if (leftExtended.count(read)) _innerReads.insert(read);
		}
		///
		
		if (path.size() >= Constants::minReadsInContig)
		{
			Logger::get().debug() << "Assembled contig with " 
				<< path.size() << " reads";
			_readLists.push_back(std::move(path));
		}
	}

	this->convertToContigs();
	_progress.setDone();
	Logger::get().info() << "Assembled " << _contigPaths.size() << " draft contigs";
}

void Extender::convertToContigs()
{
	for (auto& readsList : _readLists)
	{
		ContigPath path;
		path.name = "contig_" + std::to_string(_contigPaths.size() + 1);

		for (size_t i = 0; i < readsList.size() - 1; ++i)
		{
			bool found = false;
			OverlapRange readsOvlp;

			for (auto& ovlp : _ovlpContainer.lazySeqOverlaps(readsList[i]))
			{
				if (ovlp.extId == readsList[i + 1]) 
				{
					readsOvlp = ovlp;
					found = true;
					break;
				}
			}
			if (!found) throw std::runtime_error("Ovlp not found!");

			path.sequences.push_back(_readsContainer.getSeq(readsList[i]));
			path.overlaps.push_back(readsOvlp);
		}
		path.sequences.push_back(_readsContainer.getSeq(readsList.back()));
		_contigPaths.push_back(std::move(path));
	}
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
