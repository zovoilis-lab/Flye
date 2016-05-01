//(c) 2016 by Authors
//This file is a part of ABruijn program.
//Released under the BSD license (see LICENSE file)

#include <set>
#include <iostream>
#include <cassert>

#include "utility.h"
#include "overlap.h"

void OverlapDetector::findAllOverlaps(const VertexIndex& vertexIndex,
									  const SequenceContainer& seqContainer)
{
	LOG_PRINT("Finding overlaps");
	_overlapMatrix.clear();

	const size_t MATRIX_DIM = seqContainer.getIndex().size();
	_overlapMatrix.resize(MATRIX_DIM);
	for (size_t i = 0; i < MATRIX_DIM; ++i) _overlapMatrix[i].assign(MATRIX_DIM, 0);

	for (auto& seqHash : seqContainer.getIndex())
	{
		auto detectedOverlaps = this->getReadOverlaps(seqHash.first, 
													  vertexIndex, seqContainer);
		_overlapIndex[seqHash.first];	//empty vector by default
		for (auto ovlp : detectedOverlaps)
		{
			//detected overlap
			_overlapMatrix[ovlp.curId][ovlp.extId] = true;
			_overlapIndex[ovlp.curId].push_back(ovlp);

			//in opposite direction
			ovlp.reverse();
			_overlapMatrix[ovlp.curId][ovlp.extId] = true;
			_overlapIndex[ovlp.curId].push_back(ovlp);

			//on reverse strands
			auto curLen = seqContainer.seqLen(ovlp.curId);
			auto extLen = seqContainer.seqLen(ovlp.extId);
			ovlp.complement(curLen, extLen);
			_overlapMatrix[ovlp.curId][ovlp.extId] = true;
			_overlapIndex[ovlp.curId].push_back(ovlp);

			//opposite again
			ovlp.reverse();
			_overlapMatrix[ovlp.curId][ovlp.extId] = true;
			_overlapIndex[ovlp.curId].push_back(ovlp);
		}
	}
	_overlapMatrix.clear();

	for (auto& seqPair : _overlapIndex)
	{
		for (auto& ovlp : seqPair.second)
		{
			bool found = false;
			for (auto& extOvlp : _overlapIndex[ovlp.extId])
			{
				if (extOvlp.extId == seqPair.first) found = true;
			}
			if (!found) DEBUG_PRINT("not found!");
		}
	}
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
OverlapDetector::getReadOverlaps(FastaRecord::Id currentReadId, 
						 		 const VertexIndex& vertexIndex,
								 const SequenceContainer& seqContainer)
{
	auto& readIndex = vertexIndex.getIndexByRead();
	auto& kmerIndex = vertexIndex.getIndexByKmer();
	if (!readIndex.count(currentReadId)) return std::vector<OverlapRange>();
	
	std::unordered_map<FastaRecord::Id, 
					   std::vector<OverlapRange>> activePaths;
		
	auto curLen = seqContainer.seqLen(currentReadId);
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
			if (_overlapMatrix[extReadPos.readId]
							  [currentReadId]) continue;

			int32_t extPos = extReadPos.position;
			auto& extPaths = activePaths[extReadPos.readId];
			auto extLen = seqContainer.getIndex().at(extReadPos.readId)
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
		auto extLen = seqContainer.seqLen(ap.first);
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
			detectedOverlaps.push_back(maxOverlap);
		}
	}

	return detectedOverlaps;
}
