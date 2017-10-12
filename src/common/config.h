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
	const int repeatCoverageRate = 500;

	//overlap constants
	const int closeJumpRate = 10;
	const int farJumpRate = 2;
	const int overlapDivergenceRate = 5;
	const int penaltyWindow = 10;

	//chimera detection
	const int maxCoverageDropRate = 10;
	const int chimeraWindow = 100;

	//extension
	const int minReadsInContig = 2;
	const int maxInnerFraction = 10;

	//repeat graph
	const int maxSeparation = 500;
	const int trustedEdgeLength = 10000;
	const float minRepeatResSupport = 0.5f;
	const int outPathsRatio = 5;
	const int readCovRate = 100;

	//index construction parameters
	const int assembleKmerSample = 5;
	const int assembleGap = 50;
	const int repeatGraphKmerSample = 5;
	const int repeatGraphGap = 50;
	const int repeatGraphMaxKmer = 500;
	const int readAlignKmerSample = 5;
	const int readAlignGap = 50;
	const int readAlignMaxKmer = 500;
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
