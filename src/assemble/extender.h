//(c) 2016 by Authors
//This file is a part of ABruijn program.
//Released under the BSD license (see LICENSE file)

#pragma once

#include <deque>

#include "../sequence/sequence_container.h"
#include "../sequence/overlap.h"
#include "../sequence/contig_generator.h"
#include "chimera.h"

class Extender
{
public:
	Extender(const SequenceContainer& readsContainer, 
			 OverlapContainer& ovlpContainer,
			 int coverage, int genomeSize):
		_readsContainer(readsContainer), 
		_ovlpContainer(ovlpContainer),
		_chimDetector(readsContainer, ovlpContainer),
		_coverage(coverage), _genomeSize(genomeSize),
		_progress(genomeSize)
	{}

	void assembleContigs();
	const std::vector<ContigPath>& getContigPaths() const
		{return _contigPaths;}

private:
	typedef std::vector<FastaRecord::Id> ReadsList;

	const SequenceContainer& _readsContainer;
	OverlapContainer& _ovlpContainer;
	ChimeraDetector   _chimDetector;
	const int 		  _coverage;
	const int 		  _genomeSize;
	ProgressPercent   _progress;

	FastaRecord::Id stepRight(FastaRecord::Id readId);
	ReadsList extendContig(FastaRecord::Id startingRead);

	int   countRightExtensions(FastaRecord::Id readId);
	bool  extendsRight(const OverlapRange& ovlp);
	void  convertToContigs();

	std::vector<ReadsList> 		_readLists;
	std::vector<ContigPath> 	_contigPaths;
	//std::unordered_set<FastaRecord::Id>		 _coveredReads;
	std::unordered_set<FastaRecord::Id>		 _innerReads;
	//std::unordered_set<FastaRecord::Id>		 _usedReads;
	bool _rightExtension;
	int  _assembledSequence;
};
