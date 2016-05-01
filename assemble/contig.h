//(c) 2016 by Authors
//This file is a part of ABruijn program.
//Released under the BSD license (see LICENSE file)

#pragma once

#include <vector>

#include "fasta.h"
#include "extender.h"
#include "overlap.h"

class ContigGenerator
{
public:
	ContigGenerator(int maxJump, const Extender& extender, 
					const OverlapDetector& overlapDetector,
					const VertexIndex& vertexIndex,
					const SequenceContainer& seqContainer):
		_maximumJump(maxJump),
		_extender(extender), _overlapDetector(overlapDetector),
		_vertexIndex(vertexIndex), _seqContainer(seqContainer) {}
	
	void generateContigs();
	void outputContigs(const std::string& fileName);
	
private:
	int _maximumJump;

	const Extender& _extender;
	const OverlapDetector& _overlapDetector;
	const VertexIndex& _vertexIndex;
	const SequenceContainer& _seqContainer;

	std::pair<int32_t, int32_t> getSwitchPositions(FastaRecord::Id leftRead, 
												   FastaRecord::Id rightRead,
												   int32_t prevSwitch);

	std::vector<std::vector<FastaRecord>> _contigs;
};
