#include "overlap.h"
#include <set>
#include <iostream>
#include "utility.h"

void OverlapDetector::findAllOverlaps(const VertexIndex& vertexIndex,
									  const SequenceContainer& seqContainer)
{
	DEBUG_PRINT("Finding overlaps");
	for (auto& seqHash : seqContainer.getIndex())
	{
		this->getReadOverlaps(seqHash.first, vertexIndex, seqContainer);
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
	//if (curPrev == curNext)	//??
	//	return J_INCONS;
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
	//return J_INCONS;
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

//Getting all possible overlaps to the right direction
//based on the shared kmers (common jump-paths)
void OverlapDetector::getReadOverlaps(FastaRecord::ReadIdType currentReadId, 
						 			  const VertexIndex& vertexIndex,
									  const SequenceContainer& seqContainer)
{
	auto& readIndex = vertexIndex.getIndexByRead();
	auto& kmerIndex = vertexIndex.getIndexByKmer();
	if (!readIndex.count(currentReadId)) return;
	
	std::unordered_map<FastaRecord::ReadIdType, 
					   std::vector<OverlapRange>> activePaths;
		
	auto curLen = seqContainer.getIndex().at(currentReadId)
											.sequence.length();
	//for all kmers in this read
	for (auto& curKmerPos : readIndex.at(currentReadId))
	{
		auto curPos = curKmerPos.position;
		//for all other occurences of this kmer (extension candidates)
		for (auto& extReadPos : kmerIndex.at(curKmerPos.kmer))
		{
			if (extReadPos.readId == currentReadId) continue;

			auto extPos = extReadPos.position;
			auto& extPaths = activePaths[extReadPos.readId];
			auto extLen = seqContainer.getIndex().at(extReadPos.readId)
														.sequence.length();

			if (extPaths.empty() && this->goodStart(curPos, extPos, curLen, extLen))
			{
				extPaths.push_back(OverlapRange(curPos, extPos));
				continue;
			}

			//searching for longest possible extension
			size_t maxCloseId = 0;
			size_t maxFarId = 0;
			bool extendsClose = false;
			bool extendsFar = false;
			std::set<size_t> eraseMarks;
			for (size_t ovlpId = 0; ovlpId < extPaths.size(); ++ovlpId)
			{
				int jumpResult = this->jumpTest(extPaths[ovlpId].curEnd, curPos,
												extPaths[ovlpId].extEnd, extPos);
				int32_t jumpLength = curPos - extPaths[ovlpId].curEnd;

				switch (jumpResult)
				{
					//TODO: check, if this affects the preformance
					//maybe we really need to erase "dead ends"
					case J_END:
					case J_INCONS:
						//eraseMarks.insert(ovlpId);
						break;
					case J_CLOSE:
						eraseMarks.insert(ovlpId);
						extendsClose = true;
						if (jumpLength > curPos - extPaths[maxCloseId].curEnd)
							maxCloseId = ovlpId;	
						break;
					case J_FAR:
						extendsFar = true;
						if (jumpLength > curPos - extPaths[maxFarId].curEnd)
							maxFarId = ovlpId;	
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
			//if no extensions possible, start a new path
			if (!extendsClose && !extendsFar && 
				this->goodStart(curPos, extPos, curLen, extLen))
				extPaths.push_back(OverlapRange(curPos, extPos));
			//cleaning up
			for (auto itEraseId = eraseMarks.rbegin(); 
				 itEraseId != eraseMarks.rend(); ++itEraseId)
			{
				extPaths[*itEraseId] = extPaths.back();
				extPaths.pop_back();
			}
		}
	} //end loop over kmers in current read
	
	std::unordered_map<FastaRecord::ReadIdType, OverlapRange> detectedOverlaps;
	for (auto& ap : activePaths)
	{
		auto extLen = seqContainer.getIndex().at(ap.first)
											.sequence.length();
		for (auto& ovlp : ap.second)
			if (this->overlapTest(ovlp, curLen, extLen) &&
				detectedOverlaps[ap.first].curRange() < ovlp.curRange())
				detectedOverlaps[ap.first] = ovlp;
	}

	std::cout << "Overlaps for " 
			  << seqContainer.getIndex().at(currentReadId).description 
			  << std::endl;
	for (auto& ovlp : detectedOverlaps)
	{
		std::cout << "\twtih " 
				  << seqContainer.getIndex().at(ovlp.first).description 
				  << " : " << ovlp.second.curBegin << "," << ovlp.second.curEnd 
				  << " : " << ovlp.second.extBegin << "," << ovlp.second.extEnd 
				  << std::endl;
	}
	//return detectedOverlaps;
}
