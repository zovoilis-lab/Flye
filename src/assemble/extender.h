//(c) 2016 by Authors
//This file is a part of ABruijn program.
//Released under the BSD license (see LICENSE file)

#pragma once

#include <deque>

#include "../sequence/sequence_container.h"
#include "../sequence/overlap.h"
#include "chimera.h"

struct ContigPath
{
	ContigPath(): 
		circular(false) {}
	ContigPath(const ContigPath& other):
		reads(other.reads),
		circular(other.circular)
	{}
	ContigPath(ContigPath && other): 
		circular(other.circular) 
	{reads.swap(other.reads);}

	std::vector<FastaRecord::Id> reads;
	bool circular;
};

class Extender
{
public:
	Extender(const SequenceContainer& readsContainer, 
			 OverlapContainer& ovlpContainer,
			 int coverage):
		_readsContainer(readsContainer), 
		_ovlpContainer(ovlpContainer),
		_coverage(coverage),
		_chimDetector(coverage, readsContainer, ovlpContainer)
	{}

	void assembleContigs();
	const std::vector<ContigPath>& getContigPaths() const
		{return _contigPaths;}

private:
	const SequenceContainer& _readsContainer;
	OverlapContainer& _ovlpContainer;
	const int 		  _coverage;
	ChimeraDetector   _chimDetector;

	FastaRecord::Id stepRight(FastaRecord::Id readId);
	ContigPath extendContig(FastaRecord::Id startingRead);

	//int   rightMultiplicity(FastaRecord::Id readId);
	//bool  isBranching(FastaRecord::Id readId);
	int   countRightExtensions(FastaRecord::Id readId);
	bool  extendsRight(const OverlapRange& ovlp);
	//bool  coversRight(const OverlapRange& ovlp);
	//bool  stepAhead(FastaRecord::Id);
	//bool  resolvesRepeat(FastaRecord::Id leftRead, FastaRecord::Id rightRead);

	//void  coveredReads(const std::unordered_set<FastaRecord::Id>& allReads,
	//				   FastaRecord::Id startRead, 
	//				   std::unordered_set<FastaRecord::Id>& result);
	//bool majorClusterAgreement(FastaRecord::Id leftRead,
	//						   FastaRecord::Id rightRead);

	std::vector<ContigPath> 				 _contigPaths;
	std::unordered_set<FastaRecord::Id>		 _visitedReads;
	//std::unordered_set<FastaRecord::Id>		 _chromosomeStart;
	//bool 									 _overlapsStart;
	bool 									 _rightExtension;
	int 									 _minimumShift;
	std::unordered_map<FastaRecord::Id, int> _readsMultiplicity;
	std::unordered_map<FastaRecord::Id,
					   std::unordered_set<FastaRecord::Id>> _maxClusters;
};
