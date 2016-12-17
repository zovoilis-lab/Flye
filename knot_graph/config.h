//(c) 2016 by Authors
//This file is a part of ABruijn program.
//Released under the BSD license (see LICENSE file)

#pragma once

#include <unistd.h>

namespace Constants
{
	//global constants
	const int maximumJump = 1500;
	const int maximumOverhang = 500;

	//kmer enumeration
	const int hardMinCoverageRate = 10;
	const int repeatCoverageRate = 10;

	//overlap constants
	const int closeJumpRate = 100;
	const int farJumpRate = 2;
	const int overlapDivergenceRate = 5;

	//chimera detection
	const int maxCoverageDropRate = 5;
	const int chimeraWindow = 500;

	//extension
	const int minReadsInContig = 10;
	const int minExtensionsRate = 10;
	const int shiftToReadLen = 40;
	const float minGoodReads = 0.05f;
}

struct Parameters
{
	static int minimumOverlap;
	static size_t kmerSize;
	static size_t numThreads;
};
