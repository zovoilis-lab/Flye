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

	size_t preCountSize = 1024 * 1024 * 1024;	//1G by default
	if (Parameters::get().kmerSize > 15)		//4 ^ 15 = 1G
	{
		preCountSize *= 4 * 4;					//16G in case of larger k-mers
	}
	auto preCounters = new std::atomic<unsigned char>[preCountSize];
	for (size_t i = 0; i < preCountSize; ++i) preCounters[i] = 0;

	std::vector<FastaRecord::Id> allReads;
	for (auto& hashPair : _seqContainer.getIndex())
	{
		allReads.push_back(hashPair.first);
	}

	//first pass: filling up naive hash counting filter
	if (_outputProgress) Logger::get().info() << "Counting kmers (1/2):";
	std::function<void(const FastaRecord::Id&)> preCountUpdate = 
	[&preCounters, hardThreshold, this, preCountSize] 
		(const FastaRecord::Id& readId)
	{
		if (!readId.strand()) return;
		
		for (auto kmerPos : IterKmers(_seqContainer.getSeq(readId)))
		{
			kmerPos.kmer.standardForm();
			size_t kmerBucket = kmerPos.kmer.hash() % preCountSize;

			unsigned char expected = 0;
			while (true)
			{
				expected = preCounters[kmerBucket]; 
				if (expected == std::numeric_limits<unsigned char>::max()) 
				{
					break;
				}
				if (preCounters[kmerBucket]
						.compare_exchange_weak(expected, expected + 1))
				{
					break;
				}
			}
		}
	};
	processInParallel(allReads, preCountUpdate, 
					  Parameters::get().numThreads, _outputProgress);

	//second pass: counting kmers that have passed the filter
	if (_outputProgress) Logger::get().info() << "Counting kmers (2/2):";

	std::function<void(const FastaRecord::Id&)> countUpdate = 
	[&preCounters, hardThreshold, this, preCountSize] 
		(const FastaRecord::Id& readId)
	{
		if (!readId.strand()) return;

		for (auto kmerPos : IterKmers(_seqContainer.getSeq(readId)))
		{
			kmerPos.kmer.standardForm();
			size_t count = preCounters[kmerPos.kmer.hash() % preCountSize];
			if (count >= hardThreshold)
			{
				_kmerCounts.upsert(kmerPos.kmer, [](size_t& num){++num;}, 1);
			}
		}
	};
	processInParallel(allReads, countUpdate, 
					  Parameters::get().numThreads, _outputProgress);
	
	for (auto kmer : _kmerCounts.lock_table())
	{
		_kmerDistribution[kmer.second] += 1;
	}
	delete[] preCounters;
}


void VertexIndex::buildIndex(int minCoverage, int maxCoverage, int filterRatio)
{
	if (_outputProgress) Logger::get().info() << "Filling index table";
	
	//"Replacing" k-mer couns with k-mer index. We need multiple passes
	//to avoid peaks in memory usage during the has table extensions +
	//prevent memory fragmentation
	
	size_t kmerEntries = 0;
	size_t solidKmers = 0;
	for (auto& kmer : _kmerCounts.lock_table())
	{
		if ((size_t)minCoverage <= kmer.second && 
			kmer.second <= (size_t)maxCoverage)
		{
			kmerEntries += kmer.second;
			++solidKmers;
		}
	}
	Logger::get().debug() << "Solid kmers: " << solidKmers;
	Logger::get().debug() << "Kmer index size: " << kmerEntries;

	_kmerIndex.reserve(solidKmers);
	for (auto& kmer : _kmerCounts.lock_table())
	{
		if ((size_t)minCoverage <= kmer.second && 
			kmer.second <= (size_t)maxCoverage)
		{
			CountOrReadVector corv;
			corv.count = kmer.second;
			_kmerIndex.insert(kmer.first, corv);
		}
	}
	_kmerCounts.clear();
	_kmerCounts.reserve(0);
	
	//replacing counts with vectors
	for (auto& kmer : _kmerIndex.lock_table())
	{
		size_t count = kmer.second.count;
		kmer.second.rv = new ReadVector;
		kmer.second.rv->reserve(count);
	}

	std::function<void(const FastaRecord::Id&)> indexUpdate = 
	[minCoverage, maxCoverage, filterRatio, this] 
	(const FastaRecord::Id& readId)
	{
		if (!readId.strand()) return;

		for (auto kmerPos : IterKmers(_seqContainer.getSeq(readId)))
		{
			FastaRecord::Id targetRead = readId;
			bool revCmp = kmerPos.kmer.standardForm();
			if (revCmp)
			{
				kmerPos.position = _seqContainer.seqLen(readId) - 
										kmerPos.position -
										Parameters::get().kmerSize;
				targetRead = targetRead.rc();
			}

			//subsampling
			if (kmerPos.position % filterRatio != 0) continue;

			_kmerIndex.update_fn(kmerPos.kmer, 
				[targetRead, &kmerPos](CountOrReadVector& corv)
				{
					corv.rv->emplace_back(targetRead, kmerPos.position);
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
}


void VertexIndex::clear()
{
	for (auto kmerHash : _kmerIndex.lock_table())
	{
		delete kmerHash.second.rv;
	}
	_kmerIndex.clear();
	_kmerIndex.reserve(0);

	_kmerCounts.clear();
	_kmerCounts.reserve(0);
}
