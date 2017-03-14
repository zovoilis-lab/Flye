//(c) 2016 by Authors
//This file is a part of ABruijn program.
//Released under the BSD license (see LICENSE file)

#pragma once

#include <vector>

#include "sequence_container.h"
#include "extender.h"
#include "overlap.h"

#include "matrix.h"

class ContigGenerator
{
public:
	ContigGenerator(const Extender& extender, 
					const OverlapDetector& overlapDetector):
		_extender(extender), 
		_overlapDetector(overlapDetector),
		_vertexIndex(VertexIndex::get()), 
		_seqContainer(SequenceContainer::get()) {}
	
	void generateContigs();
	void outputContigs(const std::string& fileName);
	
private:
	const Extender& _extender;
	const OverlapDetector& _overlapDetector;
	const VertexIndex& _vertexIndex;
	const SequenceContainer& _seqContainer;

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
