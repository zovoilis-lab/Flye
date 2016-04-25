//(c) 2016 by Authors
//This file is a part of ABruijn program.
//Released under the BSD license (see LICENSE file)

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

	void extendReads();
	const std::vector<ReadPath>& getContigPaths() const
		{return _contigPaths;}

private:
	const OverlapDetector& _ovlpDetector;
	const ChimeraDetector& _chimDetector;
	const SequenceContainer& _seqContainer;

	FastaRecord::ReadIdType stepRight(FastaRecord::ReadIdType readId, 
									  FastaRecord::ReadIdType startReadId);
	bool isProperRightExtension(const OverlapDetector::OverlapRange& ovlp);
	int  countRightExtensions(FastaRecord::ReadIdType readId);

	std::vector<ReadPath> _contigPaths;
};
