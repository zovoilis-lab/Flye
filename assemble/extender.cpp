//(c) 2016 by Authors
//This file is a part of ABruijn program.
//Released under the BSD license (see LICENSE file)

#include <limits>
#include <cassert>
#include <algorithm>
#include <queue>
#include <iomanip>

#include "logger.h"
#include "extender.h"

namespace
{
	template <class T>
	T median(const std::deque<T>& vals)
	{
		if (vals.empty()) return T();
		std::deque<T> tmp = vals;
		std::nth_element(tmp.begin(), tmp.begin() + tmp.size() / 2, 
						 tmp.end());
		return tmp[tmp.size() / 2];
	}
}

ContigPath Extender::extendRead(FastaRecord::Id startRead)
{
	const size_t MAX_HISTORY = 3;

	ContigPath contigPath;
	FastaRecord::Id curRead = startRead;
	contigPath.reads.push_back(curRead);
	_visitedReads.insert(curRead);
	_visitedReads.insert(curRead.rc());
	bool rightExtension = true;

	Logger::get().debug() << "Start Read: " << 
				_seqContainer.getIndex().at(startRead).description;

	std::unordered_set<FastaRecord::Id> curPathVisited;
	_coverageHistory.clear();
	_coverageHistory.push_back(this->countRightExtensions(startRead));

	while(true)
	{
		FastaRecord::Id extRead = this->stepRight(curRead, startRead);

		if (extRead == startRead)	//circular
		{
			Logger::get().debug() << "Circular contig";
			contigPath.circular = true;
			break;
		}

		if (extRead != FastaRecord::ID_NONE) 
		{
			Logger::get().debug() << "Extension: " << 
				    	_seqContainer.getIndex().at(extRead).description;
			if (curPathVisited.count(extRead)) 
				Logger::get().debug() << "Visited in this path";
			if (_visitedReads.count(extRead)) 
				Logger::get().debug() << "Visited globally";
		}
		else
		{
			Logger::get().debug() << "No extension found"; 
		}

		if (_visitedReads.count(extRead) || extRead == FastaRecord::ID_NONE)
		{
			if (rightExtension)
			{
				Logger::get().debug() << "Changing direction";
				rightExtension = false;
				extRead = startRead.rc();
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

		_coverageHistory.push_back(this->countRightExtensions(extRead));
		if (_coverageHistory.size() > MAX_HISTORY) _coverageHistory.pop_front();

		_visitedReads.insert(extRead);
		_visitedReads.insert(extRead.rc());
		curPathVisited.insert(extRead);
		curPathVisited.insert(extRead.rc());

		curRead = extRead;
		contigPath.reads.push_back(curRead);
	}



	Logger::get().debug() << "Made " << contigPath.reads.size() - 1 
						  << " extensions";
	return contigPath;
}

void Extender::assembleContigs()
{
	const int MIN_CONTIG_SIZE = 20;	//TODO: a smarter filter
	Logger::get().info() << "Extending reads";
	_visitedReads.clear();

	for (auto& indexPair : _seqContainer.getIndex())
	{	
		if (_visitedReads.count(indexPair.first) ||
			_chimDetector.isChimeric(indexPair.first) ||
			this->countRightExtensions(indexPair.first) == 0) continue;
			//this->isBranching(indexPair.first)) continue;

		//additionaly, should not overlap with any of visited reads
		bool ovlpsVisited = false;
		for (auto& ovlp : _ovlpDetector.getOverlapIndex().at(indexPair.first))
		{
			if (_visitedReads.count(ovlp.extId)) ovlpsVisited = true;
		}
		if (ovlpsVisited) continue;

		ContigPath path = this->extendRead(indexPair.first);
		if (path.reads.size() > MIN_CONTIG_SIZE)
		{
		
			//marking visited reads
			std::unordered_set<FastaRecord::Id> leftSupported;
			std::unordered_set<FastaRecord::Id> rightSupported;
			for (auto& readId : path.reads)
			{
				for (auto& ovlp : _ovlpDetector.getOverlapIndex().at(readId))
				{
					if (this->isProperRightExtension(ovlp))
					{
						leftSupported.insert(ovlp.extId);
						rightSupported.insert(ovlp.extId.rc());
					}
					if (this->isProperLeftExtension(ovlp))
					{
						rightSupported.insert(ovlp.extId);
						leftSupported.insert(ovlp.extId.rc());
					}
				}
			}
			for (auto readId : leftSupported)
			{
				if (rightSupported.count(readId))
				{
					_visitedReads.insert(readId);
					_visitedReads.insert(readId.rc());
				}
			}
			
			_contigPaths.push_back(std::move(path));
		}
	}

	Logger::get().info() << "Assembled " << _contigPaths.size() << " contigs";
}

float Extender::branchIndex(FastaRecord::Id readId)
{
	int curEstimate = median(_coverageHistory);
	return (double)this->countRightExtensions(readId) / curEstimate;
	/*
	auto& overlaps = _ovlpDetector.getOverlapIndex().at(readId);
	int totalNewExt = 0;;
	int curExtensions = 0;
	for (auto& ovlp : overlaps)
	{
		if (this->isProperRightExtension(ovlp))
		{
			++curExtensions;
			totalNewExt += this->countRightExtensions(ovlp.extId);
		}
	}
	if (curExtensions == 0 || totalNewExt == 0) return 0.0f;

	float avgNew = (float)totalNewExt / curExtensions;
	float ratio = curExtensions / avgNew;
	return ratio;*/
}

//makes one extension to the right
FastaRecord::Id Extender::stepRight(FastaRecord::Id readId, 
									FastaRecord::Id startReadId)
{
	auto& overlaps = _ovlpDetector.getOverlapIndex().at(readId);
	std::unordered_set<FastaRecord::Id> extensions;

	for (auto& ovlp : overlaps)
	{
		assert(ovlp.curId != ovlp.extId);
		if (this->isProperRightExtension(ovlp)) 
		{
			if (ovlp.extId == startReadId) return startReadId;
			if (this->countRightExtensions(ovlp.extId) > 0)
			{
				extensions.insert(ovlp.extId);
			}
		}
	}
	
	//rank extension candidates
	std::unordered_map<FastaRecord::Id, 
					   std::tuple<int, int, int, int>> supportIndex;

	Logger::get().debug() << "Branch index " << this->branchIndex(readId);
	for (auto& extCandidate : extensions)
	{
		int leftSupport = 0;
		int rightSupport = 0;
		int ovlpSize = 0;
		int ovlpShift = 0;
		for (auto& ovlp : _ovlpDetector.getOverlapIndex().at(extCandidate))
		{
			if (ovlp.extId == readId) ovlpSize = ovlp.curRange();
			if (ovlp.extId == readId) ovlpShift = -ovlp.rightShift;
			if (!extensions.count(ovlp.extId)) continue;

			if (this->isProperRightExtension(ovlp)) ++rightSupport;
			if (this->isProperLeftExtension(ovlp)) ++leftSupport;
		}
		int minSupport = std::min(leftSupport, rightSupport);
		int endsRepeat = 1 - this->isBranching(extCandidate);
		if (!this->isBranching(readId))
		{
			supportIndex[extCandidate] = std::make_tuple(endsRepeat, minSupport, 
													 	 rightSupport, ovlpSize);
		}
		else
		{
			int startsRepeat = 1 - this->isBranching(extCandidate.rc());
			supportIndex[extCandidate] = std::make_tuple(startsRepeat, minSupport, 
														 rightSupport, ovlpSize);
		}
		Logger::get().debug() << "\t" 
				    << _seqContainer.getIndex().at(extCandidate).description
					<< "\t" << leftSupport << "\t" << rightSupport << "\t"
					<< std::fixed << std::setprecision(2) 
					<< this->branchIndex(extCandidate) << "\t" << ovlpSize
					<< "\t" << ovlpShift;
	}

	auto bestSupport = std::make_tuple(0, 0, 0, 0);
	auto bestExtension = FastaRecord::ID_NONE;
	for (auto& extCandidate : extensions)
	{
		if (extCandidate == startReadId) return startReadId;

		if (supportIndex[extCandidate] > bestSupport)
		{
			bestSupport = supportIndex[extCandidate];
			bestExtension = extCandidate;
		}
	}


	if (bestExtension != FastaRecord::ID_NONE)
	{
		if (_chimDetector.isChimeric(bestExtension))
			Logger::get().debug() << "Chimeric extension! ";
		if (this->isBranching(bestExtension))
			Logger::get().debug() << "Branching extension! ";
	}
	return bestExtension;
}

int Extender::countRightExtensions(FastaRecord::Id readId)
{
	int count = 0;
	for (auto& ovlp : _ovlpDetector.getOverlapIndex().at(readId))
	{
		if (this->isProperRightExtension(ovlp)) ++count;
	}
	return count;
}

bool Extender::isBranching(FastaRecord::Id readId)
{
	return this->branchIndex(readId) > 2.0f;
}

//Checks if read is extended to the right
bool Extender::isProperRightExtension(const OverlapRange& ovlp)
{
	return !_chimDetector.isChimeric(ovlp.extId) && 
		   ovlp.rightShift > _maximumJump;
}

//Checks if read is extended to the left
bool Extender::isProperLeftExtension(const OverlapRange& ovlp)
{
	return !_chimDetector.isChimeric(ovlp.extId) &&
		   ovlp.leftShift < -_maximumJump;
}
