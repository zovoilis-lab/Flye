#pragma once

#include "kmer_index.h"
#include "fasta.h"

class OverlapDetector
{
public:
	OverlapDetector(unsigned int maximumJump, unsigned int minimumOverlap,
					unsigned int maximumOverhang):
		_maximumJump(maximumJump), _minimumOverlap(minimumOverlap),
		_maximumOverhang(maximumOverhang) 
	{}

	void findAllOverlaps(const VertexIndex& vertexIndex);
	
private:
	void getReadOverlaps(FastaRecord::ReadIdType readId, 
						 const VertexIndex& vertexIndex);
	bool goodStart(uint32_t currentPos, uint32_t extensionPos);
	int  jumpTest(uint32_t currentPrev, uint32_t currentNext,
				  uint32_t extensionPrev, uint32_t extensionNext);

	unsigned int _maximumJump;
	unsigned int _minimumOverlap;
	unsigned int _maximumOverhang;
};
