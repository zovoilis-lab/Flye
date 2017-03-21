//(c) 2016 by Authors
//This file is a part of ABruijn program.
//Released under the BSD license (see LICENSE file)

#pragma once

#include <string>
#include <unordered_map>
#include <unordered_set>
#include <map>
#include <vector>
#include <iostream>

#include <cuckoohash_map.hh>

#include "kmer.h"
#include "sequence_container.h"

class VertexIndex
{
public:
	~VertexIndex()
	{
		this->clear();
	}
	VertexIndex(const SequenceContainer& seqContainer):
		_seqContainer(seqContainer) {}

	VertexIndex(const VertexIndex&) = delete;
	void operator=(const VertexIndex&) = delete;

	struct ReadPosition
	{
		ReadPosition(FastaRecord::Id readId, int32_t position):
			readId(readId), position(position) {}
		FastaRecord::Id readId;
		int32_t position;
	};

	typedef std::vector<ReadPosition> ReadVector;
	typedef std::map<int, int> KmerDistribution;

	void countKmers(size_t hardThreshold);
	void buildIndex(int minCoverage, int maxCoverage, int filterRatio);
	void clear();

	const ReadVector& byKmer(Kmer kmer) const
		{return *_kmerIndex[kmer];}

	bool isSolid(Kmer kmer) const
		{return _kmerIndex.contains(kmer);}

	int numSolid() const 
		{return _kmerIndex.size();}

	bool isRepetitive(Kmer kmer) const
		{return _repetitiveKmers.count(kmer);}

	const KmerDistribution& getKmerHist() const
		{return _kmerDistribution;}


private:
	const SequenceContainer& _seqContainer;

	void addFastaSequence(const FastaRecord& fastaRecord);

	cuckoohash_map<Kmer, ReadVector*> _kmerIndex;
	std::unordered_set<Kmer>		_repetitiveKmers;
	KmerDistribution 			 	_kmerDistribution;
	cuckoohash_map<Kmer, size_t>  	_kmerCounts;
};
