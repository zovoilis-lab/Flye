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
		OverlapRange(FastaRecord::ReadIdType curId, FastaRecord::ReadIdType extId, 
					 int32_t curInit, int32_t extInit): 
			curId(curId), curBegin(curInit), curEnd(curInit), 
			extId(extId), extBegin(extInit), extEnd(extInit){}
		int32_t curRange() const {return curEnd - curBegin;}
		int32_t extRange() const {return extEnd - extBegin;}

		//current read
		FastaRecord::ReadIdType curId;
		int32_t curBegin;
		int32_t curEnd;
		//extension read
		FastaRecord::ReadIdType extId;
		int32_t extBegin;
		int32_t extEnd;
	};

	typedef std::unordered_map<FastaRecord::ReadIdType, 
					   std::vector<OverlapRange>> OverlapIndex;
	const OverlapIndex& getOverlapIndex() const {return _overlapIndex;}
private:
	std::vector<OverlapRange> getReadOverlaps(FastaRecord::ReadIdType readId, 
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

	OverlapIndex _overlapIndex;
};
