//(c) 2016 by Authors
//This file is a part of ABruijn program.
//Released under the BSD license (see LICENSE file)

#pragma once

#include <string>
#include <unordered_map>
#include <map>
#include <vector>
#include <iostream>

#include <cuckoohash_map.hh>

#include "kmer.h"
#include "sequence_container.h"

class VertexIndex
{
public:
	static VertexIndex& get()
	{
		static VertexIndex instance;
		return instance;
	}
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

	void countKmers(const SequenceContainer& seqContainer,
					size_t hardThreshold);
	void buildIndex(int minCoverage, int maxCoverage);

	void setKmerSize(unsigned int size);
	unsigned int getKmerSize() const {return _kmerSize;}

	const ReadVector& byKmer(Kmer kmer) const
		{return *_kmerIndex[kmer];}
	bool isSolid(Kmer kmer) const
		{return _kmerIndex.contains(kmer);}
	const KmerDistribution& getKmerHist() const
		{return _kmerDistribution;}

private:
	VertexIndex();
	void addFastaSequence(const FastaRecord& fastaRecord);

	cuckoohash_map<Kmer, ReadVector*> _kmerIndex;
	unsigned int 					_kmerSize;
	KmerDistribution 			 	_kmerDistribution;
	cuckoohash_map<Kmer, size_t>  	_kmerCounts;
};
