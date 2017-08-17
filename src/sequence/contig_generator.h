//(c) 2016 by Authors
//This file is a part of ABruijn program.
//Released under the BSD license (see LICENSE file)

#pragma once

#include <vector>

#include "../sequence/sequence_container.h"
#include "../sequence/overlap.h"

#include "matrix.h"


class ContigGenerator
{
public:
	typedef std::vector<FastaRecord::Id> SeqVec;

	ContigGenerator(const SequenceContainer& seqContainer):
		_seqContainer(seqContainer) {}
		//_overlapContainer(overlapContainer) {}
	
	void generateContigs(const std::vector<ContigPath>& contigs);
	void outputContigs(const std::string& fileName);
	
private:
	const SequenceContainer& _seqContainer;
	//OverlapContainer& _overlapContainer;

	struct AlignmentInfo
	{
		std::string alnOne;
		std::string alnTwo;
		int32_t startOne;
		int32_t startTwo;
	};
	
	std::vector<AlignmentInfo> generateAlignments(const ContigPath& path);
	std::pair<int32_t, int32_t> getSwitchPositions(AlignmentInfo aln,
												   int32_t prevSwitch);
	FastaRecord generateLinear(const ContigPath& path);

	std::vector<FastaRecord> _contigs;
};
