#pragma once

#include "overlap.h"
#include "fasta.h"
#include <unordered_set>

class ChimeraDetector
{
public:
	ChimeraDetector(int maxOverhang, int maxJump, int coverage):
		_maximumOverhang(maxOverhang), _maximumJump(maxJump),
		_coverage(coverage) {}

	void detectChimeras(const OverlapDetector& ovlpDetector,
						const SequenceContainer& seqContainer);
	bool isChimeric(FastaRecord::ReadIdType readId) const
		{return _chimeras.count(readId) != 0;}

private:
	std::unordered_set<FastaRecord::ReadIdType> _chimeras;
	int _maximumOverhang;
	int _maximumJump;
	int _coverage;
};
