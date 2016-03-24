#pragma once

#include "overlap.h"
#include "fasta.h"

class ChimeraDetector
{
public:
	ChimeraDetector(int maxOverhang, int maxJump):
		_maximumOverhang(maxOverhang), _maximumJump(maxJump) {}

	void detectChimeras(const OverlapDetector& ovlpDetector,
						const SequenceContainer& seqContainer);

private:
	int _maximumOverhang;
	int _maximumJump;
};
