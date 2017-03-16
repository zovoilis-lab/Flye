//(c) 2016 by Authors
//This file is a part of ABruijn program.
//Released under the BSD license (see LICENSE file)

#pragma once

#include "../sequence/overlap.h"
#include "../sequence/sequence_container.h"
#include <unordered_map>

class ChimeraDetector
{
public:
	ChimeraDetector(int coverage,
					const SequenceContainer& readContainer,
					const OverlapContainer& ovlpContainer):
		_ovlpContainer(ovlpContainer), 
		_seqContainer(readContainer),
		_coverage(coverage)
	{}

	bool isChimeric(FastaRecord::Id readId);
	int getCoverage() const {return _coverage;}

private:
	//int  estimateOverlapCoverage();
	bool testReadByCoverage(FastaRecord::Id readId);

	const OverlapContainer& _ovlpContainer;
	const SequenceContainer& _seqContainer;

	std::unordered_map<FastaRecord::Id, bool> _chimeras;
	float _coverage;
};
