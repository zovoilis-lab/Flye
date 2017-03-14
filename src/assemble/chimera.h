//(c) 2016 by Authors
//This file is a part of ABruijn program.
//Released under the BSD license (see LICENSE file)

#pragma once

#include "overlap.h"
#include "sequence_container.h"
#include <unordered_set>

class ChimeraDetector
{
public:
	ChimeraDetector(int coverage,
					const OverlapDetector& ovlpDetector):
		_ovlpDetector(ovlpDetector), 
		_seqContainer(SequenceContainer::get()),
		_coverage(coverage)
	{}

	void detectChimeras();
	bool isChimeric(FastaRecord::Id readId) const
		{return _chimeras.count(readId) != 0;}
	int getCoverage() const {return _coverage;}

private:
	int  estimateOverlapCoverage();
	bool testReadByCoverage(FastaRecord::Id readId);
	bool testReadByClusters(FastaRecord::Id readId);
	bool testSelfOverlap(FastaRecord::Id readId);

	const OverlapDetector& _ovlpDetector;
	const SequenceContainer& _seqContainer;

	std::unordered_set<FastaRecord::Id> _chimeras;
	float _coverage;
};
