//(c) 2016 by Authors
//This file is a part of ABruijn program.
//Released under the BSD license (see LICENSE file)

#pragma once

#include <deque>

#include "sequence_container.h"
#include "overlap.h"
#include "chimera.h"

struct ContigPath
{
	ContigPath(): 
		circular(false) {}
	ContigPath(ContigPath && other): 
		circular(other.circular) 
	{reads.swap(other.reads);}

	std::vector<FastaRecord::Id> reads;
	bool circular;
};

class Extender
{
public:
	Extender(const OverlapDetector& ovlpDetector, 
			 const ChimeraDetector& chimDetector,
			 int maxJump, int coverage):
		_ovlpDetector(ovlpDetector), 
		_chimDetector(chimDetector),
		_seqContainer(SequenceContainer::get()), 
		_maximumJump(maxJump),
		_coverage(coverage)
	{}

	void assembleContigs();
	const std::vector<ContigPath>& getContigPaths() const
		{return _contigPaths;}

private:
	const OverlapDetector&   _ovlpDetector;
	const ChimeraDetector&   _chimDetector;
	const SequenceContainer& _seqContainer;
	const int _maximumJump;
	const int _coverage;

	FastaRecord::Id stepRight(FastaRecord::Id readId);
	ContigPath extendContig(FastaRecord::Id startingRead);

	int   rightMultiplicity(FastaRecord::Id readId);
	bool  isBranching(FastaRecord::Id readId);
	int   countRightExtensions(FastaRecord::Id readId);
	bool  extendsRight(const OverlapRange& ovlp);
	bool  coversRight(const OverlapRange& ovlp);
	bool  stepAhead(FastaRecord::Id);
	bool  resolvesRepeat(FastaRecord::Id leftRead, FastaRecord::Id rightRead);

	void  coveredReads(const std::unordered_set<FastaRecord::Id>& allReads,
					   FastaRecord::Id startRead, 
					   std::unordered_set<FastaRecord::Id>& result);
	bool majorClusterAgreement(FastaRecord::Id leftRead,
							   FastaRecord::Id rightRead);

	std::vector<ContigPath> 				 _contigPaths;
	std::unordered_set<FastaRecord::Id>		 _visitedReads;
	std::unordered_set<FastaRecord::Id>		 _chromosomeStart;
	bool 									 _overlapsStart;
	std::unordered_map<FastaRecord::Id, int> _readsMultiplicity;
	std::unordered_map<FastaRecord::Id,
					   std::unordered_set<FastaRecord::Id>> _maxClusters;
};
