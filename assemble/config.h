//(c) 2016 by Authors
//This file is a part of ABruijn program.
//Released under the BSD license (see LICENSE file)

#pragma once

namespace Constants
{
	//global constants
	const int maxumumJump = 1500;
	const int maxumumOverhang = 500;
	const int minimumOverlap = 5000;

	//kmer enumeration
	const int hardMinCoverageRate = 10;
	const int repeatCoverageRate = 10;

	//overlap constants
	const int closeJumpRate = 100;
	const int farJumpRate = 2;

	//chimer detection
	const int maxCoverageDropRate = 5;

	//extension
	const int minReadsInContig = 10;
	const int minExtensionsRate = 10;
	const float minGoodReads = 0.05f;
}
