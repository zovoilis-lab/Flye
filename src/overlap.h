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
	
private:
	void getReadOverlaps(FastaRecord::ReadIdType readId, 
						 const VertexIndex& vertexIndex,
						 const SequenceContainer& seqContainer);

	bool  goodStart(int32_t currentPos, int32_t extensionPos, 
				   int32_t currentLength);
	bool  overlapTest(int32_t curStart, int32_t curEnd, int32_t curLen,
					  int32_t extStart, int32_t extEnd);

	enum JumpRes {J_END, J_INCONS, J_CLOSE, J_FAR};
	JumpRes jumpTest(int32_t currentPrev, int32_t currentNext,
				     int32_t extensionPrev, int32_t extensionNext);

	int _maximumJump;
	int _minimumOverlap;
	int _maximumOverhang;
};
