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
	ContigGenerator(int maxJump, const Extender& extender, 
					const OverlapDetector& overlapDetector):
		_maximumJump(maxJump),
		_extender(extender), 
		_overlapDetector(overlapDetector),
		_vertexIndex(VertexIndex::get()), 
		_seqContainer(SequenceContainer::get()) {}
	
	void generateContigs();
	void outputContigs(const std::string& fileName);
	
private:
	int _maximumJump;

	const Extender& _extender;
	const OverlapDetector& _overlapDetector;
	const VertexIndex& _vertexIndex;
	const SequenceContainer& _seqContainer;

	void initAlignmentMatrices();
	void pairwiseAlignment(const std::string& seqOne, const std::string& seqTwo,
						   std::string& outOne, std::string& outTwo);
	std::pair<int32_t, int32_t> getSwitchPositions(FastaRecord::Id leftRead, 
												   FastaRecord::Id rightRead,
												   int32_t prevSwitch);
	std::vector<FastaRecord> generateCircular(const ContigPath& path);
	std::vector<FastaRecord> generateLinear(const ContigPath& path);

	std::vector<std::vector<FastaRecord>> _contigs;
	Matrix<int32_t>	_scoreMatrix;
	Matrix<char>	_backtrackMatrix;	
};
