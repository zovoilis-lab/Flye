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
		preCountSize *= 4;						//4G in case of larger k-mers
	}
	//std::vector<std::atomic<unsigned char>> preCounters(preCountSize, 0);
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

			size_t count = 0;
			_kmerCounts.find(kmerPos.kmer, count);

			if ((int)count < minCoverage) continue;
			if ((int)count > maxCoverage)
			{
				//downsampling, so as to have approximately
				//maxCoverage k-mer instances
				//if (rand() % (int)count > maxCoverage / 10) continue;
				continue;
			}

			//all good
			_kmerIndex.insert(kmerPos.kmer, nullptr);
			_kmerIndex.update_fn(kmerPos.kmer, 
				[targetRead, &kmerPos, count](ReadVector*& vec)
				{
					if (vec == nullptr)
					{
						vec = new ReadVector;
						vec->reserve(count);
					}
					vec->emplace_back(targetRead, kmerPos.position);
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
}
