//(c) 2016 by Authors
//This file is a part of ABruijn program.
//Released under the BSD license (see LICENSE file)

#pragma once

#include <vector>

#include "../sequence/sequence_container.h"
#include "../sequence/overlap.h"
#include "extender.h"

#include "matrix.h"

class ContigGenerator
{
public:
	ContigGenerator(const Extender& extender, 
					const SequenceContainer& seqContainer,
	 				OverlapContainer& overlapContainer):
		_extender(extender), 
		_seqContainer(seqContainer),
		_overlapContainer(overlapContainer) {}
	
	void generateContigs();
	void outputContigs(const std::string& fileName);
	
private:
	const Extender& _extender;
	const SequenceContainer& _seqContainer;
	OverlapContainer& _overlapContainer;

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
	std::vector<FastaRecord> generateCircular(const ContigPath& path);
	std::vector<FastaRecord> generateLinear(const ContigPath& path);

	std::vector<std::vector<FastaRecord>> _contigs;
};
