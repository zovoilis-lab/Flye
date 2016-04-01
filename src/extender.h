#pragma once

#include "fasta.h"
#include "overlap.h"
#include "chimera.h"

class Extender
{
public:
	Extender(const OverlapDetector& ovlpDetector, 
			 const ChimeraDetector& chimDetector,
			 const SequenceContainer& seqContainer):
		_ovlpDetector(ovlpDetector), _chimDetector(chimDetector),
		_seqContainer(seqContainer)
	{}

	typedef std::vector<FastaRecord::ReadIdType> ReadPath;

	std::vector<ReadPath> extendPaths();

private:
	const OverlapDetector& _ovlpDetector;
	const ChimeraDetector& _chimDetector;
	const SequenceContainer& _seqContainer;

	void extendRead(FastaRecord::ReadIdType readId);
	void stepRight(FastaRecord::ReadIdType readId);
};
