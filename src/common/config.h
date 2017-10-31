//(c) 2016 by Authors
//This file is a part of ABruijn program.
//Released under the BSD license (see LICENSE file)

#pragma once

#include <unistd.h>

namespace Constants
{
	//global constants
	const int maximumJump = 1500;
	const int maximumOverhang = 1500;

	//kmer enumeration
	const int hardMinCoverageRate = 10;
	const int repeatCoverageRate = 10;

	//overlap constants
	const int closeJumpRate = 100;
	const int farJumpRate = 2;
	const int overlapDivergenceRate = 5;
	const int penaltyWindow = 100;

	//chimera detection
	const int maxCoverageDropRate = 5;
	const int chimeraWindow = 100;

	//extension
	const int minReadsInContig = 4;
	const int maxInnerReads = 10;
	const float maxInnerFraction = 0.25;

	//repeat graph
	const int maxSeparation = 500;
	const int trustedEdgeLength = 10000;
	const float minRepeatResSupport = 0.5f;
	const int outPathsRatio = 5;
	const int readCovRate = 10;

	//index construction parameters
	const int assembleKmerSample = 1;
	const int assembleGap = 500;
	const int repeatGraphKmerSample = 5;
	const int repeatGraphGap = 100;
	const int repeatGraphMaxKmer = 500;
	const int readAlignKmerSample = 1;
	const int readAlignGap = 500;
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
