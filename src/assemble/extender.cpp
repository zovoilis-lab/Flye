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

/*
namespace
{
	template<typename Out>
	void split(const std::string &s, char delim, Out result) {
		    std::stringstream ss;
			    ss.str(s);
				    std::string item;
					    while (std::getline(ss, item, delim)) {
							        *(result++) = item;
									    }
	}


	std::vector<std::string> split(const std::string &s, char delim) {
		    std::vector<std::string> elems;
			    split(s, delim, std::back_inserter(elems));
				    return elems;
	}
}*/

ContigPath Extender::extendContig(FastaRecord::Id startRead)
{
	Logger::get().debug() << "Start Read: " << 
				_readsContainer.seqName(startRead);

	std::unordered_set<FastaRecord::Id> currentReads;
	currentReads.insert(startRead);
	currentReads.insert(startRead.rc());

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
		for (auto& ovlp : extensions)
		{
			if (currentReads.count(ovlp.extId)) 
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
				else if(bestSuspicious == FastaRecord::ID_NONE &&
						this->countRightExtensions(ovlp.extId) > 0)
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
		if (foundExtension) 
		{
			Logger::get().debug() << "Extension: " << 
				    	_readsContainer.seqName(currentRead) << " " << innerOverlaps
						<< " " << extensions.size();

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

	/*
	int falsePositive = 0;
	int falseNegative = 0;
	int allReads = 0;
	for (auto& indexPair : _readsContainer.getIndex())
	{
		bool chimeric = _chimDetector.isChimeric(indexPair.first);
		auto tokens = split(indexPair.second.description, '_');
		if (tokens.size() < 3) continue;

		std::string chimStr = tokens[2];
		bool truChim = chimStr.back() == '1';

		++allReads;
		if (truChim && !chimeric) ++falseNegative;
		if (!truChim && chimeric) 
		{
			Logger::get().info() << indexPair.second.description;
			++falsePositive;
		}
	}

	Logger::get().info() << "FN: " << falseNegative << " " 
		<< (float)falseNegative / allReads;
	Logger::get().info() << "FP: " << falsePositive << " " 
		<< (float)falsePositive / allReads;

	exit(0);*/

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

		auto& overlaps = _ovlpContainer.lazySeqOverlaps(indexPair.first);
		int numOvlp = 0;
		for (auto& ovlp : overlaps)
		{
			if (_innerReads.count(ovlp.extId)) ++numOvlp;
		}
		//float usedIndex = (float)numOvlp / overlaps.size();
		//Logger::get().debug() << "Used: " << usedIndex;
		if (numOvlp > 0) continue;

		ContigPath path = this->extendContig(indexPair.first);
		
		std::unordered_set<FastaRecord::Id> rightExtended;
		std::unordered_set<FastaRecord::Id> leftExtended;
		for (auto& readId : path.reads)
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
		
		if (path.reads.size() >= Constants::minReadsInContig)
		{
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
