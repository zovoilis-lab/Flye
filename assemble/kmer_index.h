//(c) 2016 by Authors
//This file is a part of ABruijn program.
//Released under the BSD license (see LICENSE file)

#pragma once

#include <string>
#include <unordered_map>
#include <map>
#include <vector>
#include <iostream>

#include <bf.h>

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
	typedef uint64_t KmerRepr;

	bool operator==(const Kmer& other) const
		{return this->_representation == other._representation;}
	size_t hash() const
		{return _representation;}

	struct KmerHash
	{
		std::size_t operator()(const Kmer& kmer) const
		{
			return std::hash<KmerRepr>()(kmer.hash());
		}
	};

private:
	KmerRepr _representation;
};


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

	typedef std::unordered_map<Kmer, std::vector<ReadPosition>, 
					   		   Kmer::KmerHash> KmerIndex;
	typedef std::unordered_map<FastaRecord::Id, 
					   		   std::vector<KmerPosition>> ReadIndex;
	typedef std::map<int, int> KmerDistribution;

	void		 buildKmerIndex(const SequenceContainer& seqContainer,
								size_t hardThreshold);
	void 		 setKmerSize(unsigned int size);
	void 		 applyKmerThresholds(unsigned int minCoverage, 
							 		 unsigned int maxCoverage);
	unsigned int getKmerSize() const 
		{return _kmerSize;}
	void 		 buildReadIndex();

	const KmerIndex&  getIndexByKmer() const
		{return _kmerIndex;}
	const ReadIndex&  getIndexByRead() const
		{return _readIndex;}
	const KmerDistribution& getKmerHist() const
		{return _kmerDistribution;}

	void outputCounts() const;
private:
	VertexIndex();
	void addFastaSequence(const FastaRecord& fastaRecord);

	unsigned int _kmerSize;
	KmerIndex 	 _kmerIndex;
	ReadIndex 	 _readIndex;
	KmerDistribution _kmerDistribution;
};

class CountingBloom
{
public:
	CountingBloom(size_t maxCount, size_t numHash, size_t width):
		_maxCount(maxCount)
	{
		for (size_t i = 0; i < maxCount; ++i)
		{
			_filters.emplace_back(bf::make_hasher(numHash), width);
		}
	}
	template <typename T>
	void add(const T& obj)
	{
		for (size_t i = 0; i < _maxCount; ++i)
		{
			if (!_filters[i].lookup(obj))
			{
				_filters[i].add(obj);
				break;
			}
		}
	}
	template <typename T>
	size_t count(const T& obj, size_t minCount = 0)
	{
		for (size_t i = minCount; i < _maxCount; ++i)
		{
			if (!_filters[i].lookup(obj)) return i;
		}
		return _maxCount;
	}

private:
	std::vector<bf::basic_bloom_filter> _filters;
	size_t _maxCount;
};
