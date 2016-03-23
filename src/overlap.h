#pragma once

#include "kmer_index.h"
#include "fasta.h"

class OverlapDetector
{
public:
	OverlapDetector(int maximumJump, int minimumOverlap,
					int maximumOverhang):
		_maximumJump(maximumJump), _minimumOverlap(minimumOverlap),
		_maximumOverhang(maximumOverhang) 
	{}

	void findAllOverlaps(const VertexIndex& vertexIndex, 
						 const SequenceContainer& seqContainer);

	struct OverlapRange
	{
		OverlapRange(int32_t curInit = 0, int32_t extInit = 0): 
			curBegin(curInit), curEnd(curInit), 
			extBegin(extInit), extEnd(extInit){}
		int32_t curRange() const {return curEnd - curBegin;}
		int32_t extRange() const {return extEnd - extBegin;}

		//current read
		int32_t curBegin;
		int32_t curEnd;
		//extension read
		int32_t extBegin;
		int32_t extEnd;
	};
private:

	void getReadOverlaps(FastaRecord::ReadIdType readId, 
						 const VertexIndex& vertexIndex,
						 const SequenceContainer& seqContainer);

	bool  goodStart(int32_t currentPos, int32_t extensionPos, 
				   int32_t currentLength, int32_t extensionLength);
	bool  overlapTest(const OverlapRange& ovlp, int32_t curLen, int32_t extLen);

	enum JumpRes {J_END, J_INCONS, J_CLOSE, J_FAR};
	JumpRes jumpTest(int32_t currentPrev, int32_t currentNext,
				     int32_t extensionPrev, int32_t extensionNext);

	int _maximumJump;
	int _minimumOverlap;
	int _maximumOverhang;
};
