#include "overlap.h"
#include <set>

void OverlapDetector::findAllOverlaps(const VertexIndex& vertexIndex)
{
}

namespace
{
	struct OverlapRange
	{
		OverlapRange(uint32_t curInit = 0, uint32_t extInit = 0): 
			curBegin(curInit), curEnd(curInit), 
			extBegin(extInit), extEnd(extInit){}

		//current read
		uint32_t curBegin;
		uint32_t curEnd;
		//extension read
		uint32_t extBegin;
		uint32_t extEnd;
	};
}

//Getting all possible overlaps to the right direction
//based on the shared kmers (common jump-paths)
void OverlapDetector::getReadOverlaps(FastaRecord::ReadIdType currentReadId, 
						 			  const VertexIndex& vertexIndex)
{
	auto& readIndex = vertexIndex.getIndexByRead();
	auto& kmerIndex = vertexIndex.getIndexByKmer();
	if (!readIndex.count(currentReadId)) return;
	
	std::unordered_map<FastaRecord::ReadIdType, 
					   std::vector<OverlapRange>> activePaths;
		
	//for all kmers in this read
	for (auto& curKmerPos : readIndex.at(currentReadId))
	{
		//for all other occurences of this kmer (extension candidates)
		for (auto& extReadPos : kmerIndex.at(curKmerPos.kmer))
		{
			auto curPos = curKmerPos.position;
			auto extPos = extReadPos.position;
			auto& extPaths = activePaths[extReadPos.readId];

			if (extPaths.empty() && this->goodStart(curPos, extPos))
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
				uint32_t extLength = curPos - extPaths[ovlpId].extEnd;

				switch (jumpResult)
				{
					//TODO: check, if this affects the preformance
					//maybe we really need to erase "dead ends"
					//case 0:
					//	eraseMarks.insert(ovlpId);
					//	break;
					case 2:
						eraseMarks.insert(ovlpId);
						extendsClose = true;
						if (extLength > curPos - extPaths[maxCloseId].curEnd)
							maxCloseId = ovlpId;	
						break;
					case 3:
						extendsFar = true;
						if (extLength > curPos - extPaths[maxFarId].curEnd)
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
				extPaths[-1].curEnd = curPos;
				extPaths[-1].extEnd = extPos;
			}
			//if no extensions possible, start a new path
			if (!extendsClose && !extendsFar && this->goodStart(curPos, extPos))
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
	
	for (auto ap : activePaths)
	{
		//for each overlap
		//check if it's good
		//check if it's maximal
		//record it
	}
}
