#pragma once

#include "kmer_index.h"

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
	unsigned int _maximumJump;
	unsigned int _minimumOverlap;
	unsigned int _maximumOverhang;
};
