//(c) 2016 by Authors
//This file is a part of ABruijn program.
//Released under the BSD license (see LICENSE file)

#include <set>
#include <iostream>
#include <cassert>
#include <algorithm>

#include "utility.h"
#include "overlap.h"

void OverlapDetector::findAllOverlaps()
{
	LOG_PRINT("Finding overlaps");
	_overlapMatrix.clear();

	//const size_t MATRIX_DIM = _seqContainer.getIndex().size();
	//_overlapMatrix.resize(MATRIX_DIM);
	//for (size_t i = 0; i < MATRIX_DIM; ++i) _overlapMatrix[i].assign(MATRIX_DIM, 0);

	for (auto& seqHash : _seqContainer.getIndex())
	{
		auto detectedOverlaps = this->getReadOverlaps(seqHash.first);
		_overlapIndex[seqHash.first];	//empty vector by default
		for (auto ovlp : detectedOverlaps)
		{
			//detected overlap
			_overlapMatrix.insert(std::make_tuple(ovlp.curId, ovlp.extId));
			_overlapIndex[ovlp.curId].push_back(ovlp);

			//in opposite direction
			ovlp.reverse();
			_overlapMatrix.insert(std::make_tuple(ovlp.curId, ovlp.extId));
			_overlapIndex[ovlp.curId].push_back(ovlp);

			//on reverse strands
			auto curLen = _seqContainer.seqLen(ovlp.curId);
			auto extLen = _seqContainer.seqLen(ovlp.extId);
			ovlp.complement(curLen, extLen);
			_overlapMatrix.insert(std::make_tuple(ovlp.curId, ovlp.extId));
			_overlapIndex[ovlp.curId].push_back(ovlp);

			//opposite again
			ovlp.reverse();
			_overlapMatrix.insert(std::make_tuple(ovlp.curId, ovlp.extId));
			_overlapIndex[ovlp.curId].push_back(ovlp);
		}
	}
	_overlapMatrix.clear();
}

//pre-filtering (maybe it's not needed)
bool OverlapDetector::goodStart(int32_t curPos, int32_t extPos, 
								int32_t curLen, int32_t extLen)
{	
	return 	(extPos < _maximumOverhang && curPos < curLen - _minimumOverlap) ||
		   	(curPos < _maximumOverhang && extPos < extLen - _minimumOverlap);
}

//TODO: check if we need different stop results
OverlapDetector::JumpRes 
OverlapDetector::jumpTest(int32_t curPrev, int32_t curNext,
						  int32_t extPrev, int32_t extNext)
{
	static const int CLOSE_FRAC = 8;
	static const int FAR_FRAC = 2;
	if (curNext - curPrev < _maximumJump && extNext - extPrev < _maximumJump)
	{
		if (abs((curNext - curPrev) - (extNext - extPrev)) 
			< _maximumJump / CLOSE_FRAC)
			return J_CLOSE;
		if (abs((curNext - curPrev) - (extNext - extPrev)) 
			< _maximumJump / FAR_FRAC)
			return J_FAR;
	}
	return J_END;
}


//Check if it is a proper overlap
bool OverlapDetector::overlapTest(const OverlapRange& ovlp, int32_t curLen, 
								  int32_t extLen)
{
	if (ovlp.curRange() < _minimumOverlap || ovlp.extRange() < _minimumOverlap)
		return false;
	if (abs(ovlp.curRange() - ovlp.extRange()) > _maximumJump)
		return false;
	if (std::min(ovlp.curBegin, ovlp.extBegin) > _maximumOverhang)
		return false;
	if (std::min(curLen - ovlp.curEnd, extLen - ovlp.extEnd) > _maximumOverhang)
		return false;
	return true;
}

//Getting all possible overlaps
//based on the shared kmers (common jump-paths)
std::vector<OverlapRange> 
OverlapDetector::getReadOverlaps(FastaRecord::Id currentReadId)
{
	auto& readIndex = _vertexIndex.getIndexByRead();
	auto& kmerIndex = _vertexIndex.getIndexByKmer();
	if (!readIndex.count(currentReadId)) return std::vector<OverlapRange>();
	
	std::unordered_map<FastaRecord::Id, 
					   std::vector<OverlapRange>> activePaths;
		
	auto curLen = _seqContainer.seqLen(currentReadId);
	//for all kmers in this read
	for (auto& curKmerPos : readIndex.at(currentReadId))
	{
		int32_t curPos = curKmerPos.position;
		//for all other occurences of this kmer (extension candidates)
		for (auto& extReadPos : kmerIndex.at(curKmerPos.kmer))
		{
			//don't want self-overlaps
			if (extReadPos.readId == currentReadId ||
				extReadPos.readId == currentReadId.rc()) continue;
			//maybe it's already processed
			if (_overlapMatrix.count(std::make_tuple(extReadPos.readId,
							   						 currentReadId))) continue;

			int32_t extPos = extReadPos.position;
			auto& extPaths = activePaths[extReadPos.readId];
			auto extLen = _seqContainer.getIndex().at(extReadPos.readId)
														.sequence.length();

			//searching for longest possible extension
			size_t maxCloseId = 0;
			size_t maxFarId = 0;
			bool extendsClose = false;
			bool extendsFar = false;
			std::set<size_t> eraseMarks;
			for (size_t pathId = 0; pathId < extPaths.size(); ++pathId)
			{
				JumpRes jumpResult = this->jumpTest(extPaths[pathId].curEnd, curPos,
													extPaths[pathId].extEnd, extPos);
				int32_t jumpLength = curPos - extPaths[pathId].curEnd;

				switch (jumpResult)
				{
					//TODO: check, if this affects the performance
					//maybe we really need to erase "dead ends"
					case J_END:
					case J_INCONS:
						break;
					case J_CLOSE:
						eraseMarks.insert(pathId);
						extendsClose = true;
						if (jumpLength > curPos - extPaths[maxCloseId].curEnd)
						{
							maxCloseId = pathId;	
						}
						break;
					case J_FAR:
						extendsFar = true;
						if (jumpLength > curPos - extPaths[maxFarId].curEnd)
						{
							maxFarId = pathId;	
						}
						break;
				}
			}
			//update the best close extension
			if (extendsClose)
			{
				eraseMarks.erase(maxCloseId);
				extPaths[maxCloseId].curEnd = curPos;
				extPaths[maxCloseId].extEnd = extPos;
			}
			//update the best far extension, keep the old path as a copy
			if (extendsFar)
			{
				extPaths.push_back(extPaths[maxFarId]);
				extPaths.back().curEnd = curPos;
				extPaths.back().extEnd = extPos;
			}
			//if no extensions possible (or there are no active paths), start a new path
			if (!extendsClose && !extendsFar && 
				this->goodStart(curPos, extPos, curLen, extLen))
			{
				extPaths.push_back(OverlapRange(currentReadId, extReadPos.readId,
												curPos, extPos));
			}
			//cleaning up
			for (auto itEraseId = eraseMarks.rbegin(); 
				 itEraseId != eraseMarks.rend(); ++itEraseId)
			{
				extPaths[*itEraseId] = extPaths.back();
				extPaths.pop_back();
			}
		} //end loop over kmer occurences in other reads
	} //end loop over kmers in the current read
	
	std::vector<OverlapRange> detectedOverlaps;
	for (auto& ap : activePaths)
	{
		auto extLen = _seqContainer.seqLen(ap.first);
		OverlapRange maxOverlap(0, 0, 0, 0);
		bool passedTest = false;
		for (auto& ovlp : ap.second)
		{
			if (this->overlapTest(ovlp, curLen, extLen))
			{
				passedTest = true;
				if (maxOverlap.curRange() < ovlp.curRange()) maxOverlap = ovlp;
			}
		}
		if (passedTest)
		{
			this->addOverlapShifts(maxOverlap);
			detectedOverlaps.push_back(std::move(maxOverlap));
		}
	}

	return detectedOverlaps;
}

namespace
{
	template<typename T>
	T median(std::vector<T>& vec)
	{
		std::nth_element(vec.begin(), vec.begin() + vec.size() / 2, 
						 vec.end());
		return vec[vec.size() / 2];
	}
}

void OverlapDetector::addOverlapShifts(OverlapRange& ovlp)
{
	//get shared kmers inside the overlap
	std::vector<int32_t> ovlpShifts;
	for (auto& curKmer : _vertexIndex.getIndexByRead().at(ovlp.curId))
	{
		if (ovlp.curBegin <= curKmer.position &&
			curKmer.position <= ovlp.curEnd)
		{
			for (auto& extKmer : _vertexIndex.getIndexByKmer().at(curKmer.kmer))
			{
				if (extKmer.readId == ovlp.extId &&
				    ovlp.extBegin <= extKmer.position &&
					extKmer.position <= ovlp.extEnd)
				{
					ovlpShifts.push_back(curKmer.position - extKmer.position);
				}
			}
		}
	}

	ovlp.leftShift = median(ovlpShifts);
	ovlp.rightShift = _seqContainer.seqLen(ovlp.extId) - 
					  _seqContainer.seqLen(ovlp.curId) + ovlp.leftShift;
}
