//(c) 2016 by Authors
//This file is a part of ABruijn program.
//Released under the BSD license (see LICENSE file)

#include <stdexcept>
#include <iostream>
#include <unordered_set>
#include <algorithm>

#include "vertex_index.h"
#include "../common/logger.h"
#include "../common/parallel.h"
#include "../common/config.h"


void VertexIndex::countKmers(size_t hardThreshold)
{
	Logger::get().debug() << "Hard threshold set to " << hardThreshold;
	if (hardThreshold == 0 || hardThreshold > 100) 
	{
		throw std::runtime_error("Wrong hard threshold value: " + 
								 std::to_string(hardThreshold));
	}
	Logger::get().debug() << "Started kmer counting";

	//const size_t PRE_COUNT_SIZE = 1024 * 1024 * 1024;
	const size_t PRE_COUNT_SIZE = pow(4, Parameters::get().kmerSize);
	std::vector<unsigned char> preCounters(PRE_COUNT_SIZE, 0);

	//filling up bloom filter
	if (_outputProgress) Logger::get().info() << "Counting kmers (1/2):";
	ProgressPercent bloomProg(_seqContainer.getIndex().size());
	for (auto& seqPair : _seqContainer.getIndex())
	{
		if (_outputProgress) bloomProg.advance();
		for (auto kmerPos : IterKmers(_seqContainer.getSeq(seqPair.first)))
		{
			kmerPos.kmer.standardForm();
			if (preCounters[kmerPos.kmer.hash() % PRE_COUNT_SIZE] != 
				std::numeric_limits<unsigned char>::max())
				++preCounters[kmerPos.kmer.hash() % PRE_COUNT_SIZE];
		}
	}

	//counting only kmers that have passed the filter
	if (_outputProgress) Logger::get().info() << "Counting kmers (2/2):";

	std::function<void(const FastaRecord::Id&)> countUpdate = 
	[&preCounters, hardThreshold, this, PRE_COUNT_SIZE] 
		(const FastaRecord::Id& readId)
	{
		for (auto kmerPos : IterKmers(_seqContainer.getSeq(readId)))
		{
			kmerPos.kmer.standardForm();
			size_t count = preCounters[kmerPos.kmer.hash() % PRE_COUNT_SIZE];
			if (count >= hardThreshold * 2)		//fwd and rev
			{
				_kmerCounts.upsert(kmerPos.kmer, [](size_t& num){++num;}, 1);
			}
		}
	};
	std::vector<FastaRecord::Id> allReads;
	for (auto& hashPair : _seqContainer.getIndex())
	{
		allReads.push_back(hashPair.first);
	}
	processInParallel(allReads, countUpdate, 
					  Parameters::get().numThreads, _outputProgress);
	
	for (auto kmer : _kmerCounts.lock_table())
	{
		_kmerDistribution[kmer.second / 2] += 1;	//fwd and rev
	}
}


void VertexIndex::buildIndex(int minCoverage, int maxCoverage, int filterRatio)
{
	if (_outputProgress) Logger::get().info() << "Filling index table";

	size_t kmerEntries = 0;
	for (auto kmer : _kmerCounts.lock_table())
	{
		if ((size_t)minCoverage <= kmer.second && 
			kmer.second <= (size_t)maxCoverage)
		{
			kmerEntries += kmer.second;
		}
	}
	Logger::get().debug() << "Kmer index size: " << kmerEntries;

	std::function<void(const FastaRecord::Id&)> indexUpdate = 
	[minCoverage, maxCoverage, filterRatio, this] 
	(const FastaRecord::Id& readId)
	{
		for (auto kmerPos : IterKmers(_seqContainer.getSeq(readId)))
		{
			//keep only kmers in standard form in index
			bool revCmp = kmerPos.kmer.standardForm();
			if (revCmp) continue;

			//subsampling
			int32_t samplePos = kmerPos.position;
			if (!readId.strand())	//keeping strands synchonized
			{
				samplePos = _seqContainer.seqLen(readId) - samplePos -
							Parameters::get().kmerSize;
			}
			if (samplePos % filterRatio != 0) continue;

			size_t count = 0;
			_kmerCounts.find(kmerPos.kmer, count);
			count /= 2;	//fwd and rev

			if ((int)count < minCoverage) continue;
			if ((int)count > maxCoverage)
			{
				//downsampling, so as to have approximately
				//maxCoverage k-mer instances
				if (rand() % (int)count > maxCoverage) continue;
			}

			//all good
			_kmerIndex.insert(kmerPos.kmer, nullptr);
			_kmerIndex.update_fn(kmerPos.kmer, 
				[readId, &kmerPos, count](ReadVector*& vec)
				{
					if (vec == nullptr)
					{
						vec = new ReadVector;
						vec->reserve(count);
					}
					vec->emplace_back(readId, kmerPos.position);
				});
		}
	};
	std::vector<FastaRecord::Id> allReads;
	for (auto& hashPair : _seqContainer.getIndex())
	{
		allReads.push_back(hashPair.first);
	}
	processInParallel(allReads, indexUpdate, 
					  Parameters::get().numThreads, _outputProgress);

	_kmerCounts.clear();
	_kmerCounts.reserve(0);
}


void VertexIndex::clear()
{
	for (auto kmerHash : _kmerIndex.lock_table())
	{
		delete kmerHash.second;
	}
	_kmerIndex.clear();
	_kmerIndex.reserve(0);

	_kmerCounts.clear();
	_kmerCounts.reserve(0);

	//_repetitiveKmers.clear();
	//_repetitiveKmers.reserve(0);
}
