//(c) 2016 by Authors
//This file is a part of ABruijn program.
//Released under the BSD license (see LICENSE file)

#pragma once

#include <vector>

#include "../sequence/overlap.h"

#include "matrix.h"

struct ContigPath
{
	ContigPath() {}

	std::string name;
	std::vector<DnaSequence> sequences;
	std::vector<OverlapRange> overlaps;
};

class ContigGenerator
{
public:
	void generateContigs(const std::vector<ContigPath>& contigs);
	void outputContigs(const std::string& fileName);
	FastaRecord generateLinear(const ContigPath& path);
	
private:
	struct AlignmentInfo
	{
		std::string alnOne;
		std::string alnTwo;

		int32_t startOne;
		int32_t startTwo;
	};
	
	std::vector<AlignmentInfo> generateAlignments(const ContigPath& path);
	std::pair<int32_t, int32_t> getSwitchPositions(const AlignmentInfo& aln,
												   int32_t prevSwitch);

	std::vector<FastaRecord> _contigs;
};
