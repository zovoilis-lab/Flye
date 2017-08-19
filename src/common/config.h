//(c) 2016 by Authors
//This file is a part of ABruijn program.
//Released under the BSD license (see LICENSE file)

#pragma once

#include <unistd.h>

namespace Constants
{
	const int maximumJump = 100;
	const int maximumOverhang = 100;

	//kmer enumeration
	const int hardMinCoverageRate = 10;
	const int repeatCoverageRate = 10;

	//overlap constants
	const int closeJumpRate = 10;
	const int farJumpRate = 2;
	const int overlapDivergenceRate = 5;
	const int gapJump = 100;
	const int penaltyWindow = 10;

	//chimera detection
	const int maxCoverageDropRate = 5;
	const int chimeraWindow = 100;

	//extension
	const int minReadsInContig = 2;
	const int minExtensions = 2;
	const float startReadsPercent = 0.5;

	//repeat graph
	const int maxSeparation = 500;
	const int trustedEdgeLength = 10000;
	const float minRepeatResSupport = 0.25f;

	//index construction parameters
	const int assembleKmerSample = 1;
	const int repeatGraphKmerSample = 5;
	const int repeatGraphMaxKmer = 1000;
	const int readAlignKmerSample = 1;
	const int readAlignMaxKmer = 1000;
}

struct Parameters
{
	static Parameters& get()
	{
		static Parameters param;
		return param;
	}

	int minimumOverlap;
	size_t kmerSize;
	size_t numThreads;
};
