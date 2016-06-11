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

#include "fasta.h"

static_assert(sizeof(size_t) == 8, "32-bit architectures are not supported");

class Kmer
{
public:
	Kmer(const std::string& dnaString);

	void reverseComplement();
	void appendRight(char dnaSymbol);
	void appendLeft(char dnaSymbol);
	std::string dnaRepresentation() const;
	typedef size_t KmerRepr;

	bool operator==(const Kmer& other) const
		{return this->_representation == other._representation;}
	size_t hash() const
		{return _representation;}

private:
	KmerRepr _representation;
};


namespace std
{
	template <>
	struct hash<Kmer>
	{
		std::size_t operator()(const Kmer& kmer) const
		{
			return std::hash<size_t>()(kmer.hash());
		}
	};
}


class VertexIndex
{
public:
	static VertexIndex& getInstance()
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
	struct KmerPosition
	{
		KmerPosition(Kmer kmer, int32_t position):
			kmer(kmer), position(position) {}
		Kmer kmer;
		int32_t position;
	};

	typedef std::unordered_map<Kmer, std::vector<ReadPosition>> KmerIndex;
	typedef std::unordered_map<FastaRecord::Id, 
					   		   std::vector<KmerPosition>> ReadIndex;
	typedef std::map<int, int> KmerDistribution;

	void		 countKmers(const SequenceContainer& seqContainer,
							size_t hardThreshold);
	void 		 buildIndex(const SequenceContainer& seqContainer,
							int minCoverage, int maxCoverage);
	void 		 setKmerSize(unsigned int size);

	unsigned int getKmerSize() const 
		{return _kmerSize;}

	const KmerIndex&  getIndexByKmer() const
		{return _kmerIndex;}
	const ReadIndex&  getIndexByRead() const
		{return _readIndex;}
	const KmerDistribution& getKmerHist() const
		{return _kmerDistribution;}

private:
	VertexIndex();
	void addFastaSequence(const FastaRecord& fastaRecord);

	unsigned int _kmerSize;
	KmerIndex 	 _kmerIndex;
	ReadIndex 	 _readIndex;

	KmerDistribution 			 _kmerDistribution;
	cuckoohash_map<size_t, size_t> _kmerCounts;
};
