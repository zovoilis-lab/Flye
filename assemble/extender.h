//(c) 2016 by Authors
//This file is a part of ABruijn program.
//Released under the BSD license (see LICENSE file)

#pragma once

#include "fasta.h"
#include "overlap.h"
#include "chimera.h"

struct ContigPath
{
	ContigPath(): 
		circular(false) {}
	ContigPath(ContigPath && other): 
		circular(other.circular) 
	{reads.swap(other.reads);}

	std::vector<FastaRecord::ReadIdType> reads;
	bool circular;
};

class Extender
{
public:
	Extender(const OverlapDetector& ovlpDetector, 
			 const ChimeraDetector& chimDetector,
			 const SequenceContainer& seqContainer):
		_ovlpDetector(ovlpDetector), _chimDetector(chimDetector),
		_seqContainer(seqContainer)
	{}

	//typedef std::vector<FastaRecord::ReadIdType> ReadPath;

	void assembleContigs();
	const std::vector<ContigPath>& getContigPaths() const
		{return _contigPaths;}

private:
	const OverlapDetector& _ovlpDetector;
	const ChimeraDetector& _chimDetector;
	const SequenceContainer& _seqContainer;

	FastaRecord::ReadIdType stepRight(FastaRecord::ReadIdType readId, 
									  FastaRecord::ReadIdType startReadId);
	ContigPath extendRead(FastaRecord::ReadIdType readId);

	bool  isProperRightExtension(const OverlapRange& ovlp);
	bool  isProperLeftExtension(const OverlapRange& ovlp);
	float branchIndex(FastaRecord::ReadIdType readId);
	int   countRightExtensions(FastaRecord::ReadIdType readId);

	std::vector<ContigPath> _contigPaths;
	std::unordered_set<FastaRecord::ReadIdType>	_visitedReads;
};
