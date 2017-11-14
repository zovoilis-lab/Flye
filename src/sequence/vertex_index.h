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
#include "../common/config.h"

class VertexIndex
{
public:
	~VertexIndex()
	{
		this->clear();
	}
	VertexIndex(const SequenceContainer& seqContainer):
		_seqContainer(seqContainer), _outputProgress(false) 
	{}

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
	typedef std::map<size_t, size_t> KmerDistribution;
	union CountOrReadVector
	{
		size_t count;
		ReadVector* rv;
	};

	class KmerPosIterator
	{
	public:
		KmerPosIterator(const ReadVector* rv, size_t index, bool revComp, 
						const SequenceContainer& seqContainer):
			rv(rv), index(index), revComp(revComp), 
			seqContainer(seqContainer) 
		{}

		bool operator==(const KmerPosIterator& other) const
		{
			return rv == other.rv && index == other.index;
		}
		bool operator!=(const KmerPosIterator& other) const
		{
			return !(*this == other);
		}

		ReadPosition operator*() const
		{
			ReadPosition pos = (*rv)[index];
			if (revComp)
			{
				pos.readId = pos.readId.rc();
				int32_t seqLen = seqContainer.seqLen(pos.readId);
				pos.position = seqLen - pos.position - 
							   Parameters::get().kmerSize;
			}
			return pos;
		}

		KmerPosIterator& operator++()
		{
			++index;
			return *this;
		}
	
	private:
		const ReadVector* rv;
		size_t index;
		bool revComp;
		const SequenceContainer& seqContainer;
	};

	class IterHelper
	{
	public:
		IterHelper(const ReadVector& rv, bool revComp, 
				   const SequenceContainer& seqContainer): 
			rv(rv), revComp(revComp), seqContainer(seqContainer) {}

		KmerPosIterator begin()
		{
			return KmerPosIterator(&rv, 0, revComp, seqContainer);
		}

		KmerPosIterator end()
		{
			return KmerPosIterator(&rv, rv.size(), revComp, seqContainer);
		}

	private:
		const ReadVector& rv;
		bool revComp;
		const SequenceContainer& seqContainer;
	};


	void countKmers(size_t hardThreshold);
	void buildIndex(int minCoverage, int maxCoverage, int filterRatio);
	void clear();

	IterHelper iterKmerPos(Kmer kmer) const
	{
		bool revComp = kmer.standardForm();
		return IterHelper(*_kmerIndex[kmer].rv, revComp,
						  _seqContainer);
	}

	bool isSolid(Kmer kmer) const
	{
		kmer.standardForm();
		return _kmerIndex.contains(kmer);
	}

	void outputProgress(bool set) 
	{
		_outputProgress = set;
	}

	const KmerDistribution& getKmerHist() const
	{
		return _kmerDistribution;
	}

private:
	const SequenceContainer& _seqContainer;

	void addFastaSequence(const FastaRecord& fastaRecord);

	cuckoohash_map<Kmer, CountOrReadVector> _kmerIndex;
	KmerDistribution 			 	_kmerDistribution;
	cuckoohash_map<Kmer, size_t>  	_kmerCounts;
	bool _outputProgress;
};
