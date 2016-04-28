//(c) 2016 by Authors
//This file is a part of ABruijn program.
//Released under the BSD license (see LICENSE file)

#pragma once

#include "overlap.h"
#include "fasta.h"
#include <unordered_set>

class ChimeraDetector
{
public:
	ChimeraDetector(int maxOverhang, int maxJump, 
					const OverlapDetector& ovlpDetector,
					const SequenceContainer& seqContainer):
		_ovlpDetector(ovlpDetector), _seqContainer(seqContainer),
		_maximumOverhang(maxOverhang), _maximumJump(maxJump),
		_coverage(0)
	{}

	void detectChimeras();
	bool isChimeric(FastaRecord::ReadIdType readId) const
		{return _chimeras.count(readId) != 0;}
	int getCoverage() const {return _coverage;}

private:
	bool testRead(FastaRecord::ReadIdType readId);

	const OverlapDetector& _ovlpDetector;
	const SequenceContainer& _seqContainer;

	std::unordered_set<FastaRecord::ReadIdType> _chimeras;
	int _maximumOverhang;
	int _maximumJump;
	float _coverage;
};
