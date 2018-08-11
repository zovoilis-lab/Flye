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


void VertexIndex::countKmers(size_t hardThreshold, int genomeSize)
{
	if (Parameters::get().kmerSize > 31)
	{
		throw std::runtime_error("Maximum kmer size is 31");
	}

	Logger::get().debug() << "Hard threshold set to " << hardThreshold;
	if (hardThreshold == 0)
	{
		throw std::runtime_error("Wrong hard threshold value: " + 
								 std::to_string(hardThreshold));
	}
	Logger::get().debug() << "Started kmer counting";

	size_t preCountSize = 1024 * 1024 * 1024;	//1G by default
	if (genomeSize > (int)Config::get("big_genome_threshold"))
	{
		preCountSize *= 4 * 4;					//16G in case of larger genomes
	}
	auto preCounters = new std::atomic<unsigned char>[preCountSize];
	for (size_t i = 0; i < preCountSize; ++i) preCounters[i] = 0;

	std::vector<FastaRecord::Id> allReads;
	for (auto& seq : _seqContainer.iterSeqs())
	{
		allReads.push_back(seq.id);
	}

	//first pass: filling up naive hash counting filter
	if (_outputProgress) Logger::get().info() << "Counting kmers (1/2):";
	std::function<void(const FastaRecord::Id&)> preCountUpdate = 
	[&preCounters, hardThreshold, this, preCountSize] 
		(const FastaRecord::Id& readId)
	{
		if (!readId.strand()) return;
		
		int32_t nextKmerPos = _sampleRate;
		for (auto kmerPos : IterKmers(_seqContainer.getSeq(readId)))
		{
			if (_sampleRate > 1) //subsampling
			{
				if (--nextKmerPos > 0) continue;
				nextKmerPos = _sampleRate + 
					(int32_t)((kmerPos.kmer.hash() ^ readId.hash()) % 3) - 1;
			}

			bool revCmp = kmerPos.kmer.standardForm();
			if (revCmp)
			{
				kmerPos.position = _seqContainer.seqLen(readId) - 
										kmerPos.position -
										Parameters::get().kmerSize;
			}
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

		int32_t nextKmerPos = _sampleRate;
		for (auto kmerPos : IterKmers(_seqContainer.getSeq(readId)))
		{
			if (_sampleRate > 1) //subsampling
			{
				if (--nextKmerPos > 0) continue;
				nextKmerPos = _sampleRate + 
					(int32_t)((kmerPos.kmer.hash() ^ readId.hash()) % 3) - 1;
			}

			bool revCmp = kmerPos.kmer.standardForm();
			if (revCmp)
			{
				kmerPos.position = _seqContainer.seqLen(readId) - 
										kmerPos.position -
										Parameters::get().kmerSize;
			}

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
		_repetitiveFrequency = std::max(_repetitiveFrequency, kmer.second);
	}
	delete[] preCounters;
}

namespace
{
	template <class T>
	size_t getFreq(T& histIter)
		{return histIter->first;};

	template <class T>
	size_t getCount(T& histIter)
		{return histIter->second;};

}

void VertexIndex::setRepeatCutoff(int minCoverage)
{
	size_t totalKmers = 0;
	size_t uniqueKmers = 0;
	for (auto mapPair = this->getKmerHist().begin();
		 mapPair != this->getKmerHist().end(); ++mapPair)
	{
		if (minCoverage <= (int)getFreq(mapPair))
		{
			totalKmers += getCount(mapPair) * getFreq(mapPair);
			uniqueKmers += getCount(mapPair);
		}
	}
	float meanFrequency = (float)totalKmers / (uniqueKmers + 1);
	_repetitiveFrequency = (float)Config::get("repeat_kmer_rate") * meanFrequency;
	
	size_t repetitiveKmers = 0;
	for (auto mapPair = this->getKmerHist().rbegin();
		 mapPair != this->getKmerHist().rend(); ++mapPair)
	{
		if (getFreq(mapPair) > _repetitiveFrequency)
		{
			repetitiveKmers += getCount(mapPair);
		}
	}
	float filteredRate = (float)repetitiveKmers / uniqueKmers;
	Logger::get().debug() << "Repetetive k-mer frequency: " 
						  << _repetitiveFrequency;
	Logger::get().debug() << "Filtered " << repetitiveKmers 
						  << " repetitive kmers (" <<
						  filteredRate << ")";
}


void VertexIndex::buildIndex(int minCoverage)
{
	if (_outputProgress) Logger::get().info() << "Filling index table";
	
	//"Replacing" k-mer couns with k-mer index. We need multiple passes
	//to avoid peaks in memory usage during the hash table extensions +
	//prevent memory fragmentation
	
	size_t kmerEntries = 0;
	size_t solidKmers = 0;
	for (auto& kmer : _kmerCounts.lock_table())
	{
		if ((size_t)minCoverage <= kmer.second)
		{
			kmerEntries += kmer.second;
			++solidKmers;
		}
	}
	Logger::get().debug() << "Sampling rate: " << _sampleRate;
	Logger::get().debug() << "Solid kmers: " << solidKmers;
	Logger::get().debug() << "Kmer index size: " << kmerEntries;
	Logger::get().debug() << "Mean k-mer frequency: " 
		<< (float)kmerEntries / solidKmers;

	_kmerIndex.reserve(solidKmers);
	for (auto& kmer : _kmerCounts.lock_table())
	{
		if ((size_t)minCoverage <= kmer.second)
		{
			ReadVector rv{(uint32_t)kmer.second, 0, nullptr};
			_kmerIndex.insert(kmer.first, rv);
		}
		if (kmer.second > _repetitiveFrequency)
		{
			_repetitiveKmers.insert(kmer.first, true);
		}
	}
	_kmerCounts.clear();
	_kmerCounts.reserve(0);

	_memoryChunks.push_back(new ReadPosition[INDEX_CHUNK]);
	size_t chunkOffset = 0;
	for (auto& kmer : _kmerIndex.lock_table())
	{
		if (INDEX_CHUNK < kmer.second.capacity) 
		{
			throw std::runtime_error("k-mer is too frequent");
		}
		if (INDEX_CHUNK - chunkOffset < kmer.second.capacity)
		{
			_memoryChunks.push_back(new ReadPosition[INDEX_CHUNK]);
			chunkOffset = 0;
		}
		kmer.second.data = _memoryChunks.back() + chunkOffset;
		chunkOffset += kmer.second.capacity;
	}
	//Logger::get().debug() << "Total chunks " << _memoryChunks.size()
	//	<< " wasted space: " << wasted;

	std::function<void(const FastaRecord::Id&)> indexUpdate = 
	[this] (const FastaRecord::Id& readId)
	{
		if (!readId.strand()) return;

		int32_t nextKmerPos = _sampleRate;
		int32_t seqLen = _seqContainer.seqLen(readId);
		for (auto kmerPos : IterKmers(_seqContainer.getSeq(readId)))
		{
			if (_sampleRate > 1) //subsampling
			{
				if (--nextKmerPos > 0) continue;
				nextKmerPos = _sampleRate + 
					(int32_t)((kmerPos.kmer.hash() ^ readId.hash()) % 3) - 1;
			}

			FastaRecord::Id targetRead = readId;
			bool revCmp = kmerPos.kmer.standardForm();
			if (revCmp)
			{
				kmerPos.position = _seqContainer.seqLen(readId) - 
										kmerPos.position -
										Parameters::get().kmerSize;
				targetRead = targetRead.rc();
			}

			//filter out repetitive kmers (but allow some on sequence ends)
			/*if (_repetitiveKmers.contains(kmerPos.kmer) &&
				kmerPos.position >= _flankRepeatSize &&
				kmerPos.position < seqLen - _flankRepeatSize) continue;*/

			_kmerIndex.update_fn(kmerPos.kmer, 
				[targetRead, &kmerPos](ReadVector& rv)
				{

					rv.data[rv.size] = ReadPosition(targetRead, 
													kmerPos.position);
					++rv.size;
				});
		}
	};
	std::vector<FastaRecord::Id> allReads;
	for (auto& seq : _seqContainer.iterSeqs())
	{
		allReads.push_back(seq.id);
	}
	processInParallel(allReads, indexUpdate, 
					  Parameters::get().numThreads, _outputProgress);

	Logger::get().debug() << "Sorting k-mer index";
	for (auto& kmerVec : _kmerIndex.lock_table())
	{
		std::sort(kmerVec.second.data, kmerVec.second.data + kmerVec.second.size,
				  [](const ReadPosition& p1, const ReadPosition& p2)
				  	{return p1.readId < p2.readId;});
	}
}


void VertexIndex::clear()
{
	for (auto& chunk : _memoryChunks) delete[] chunk;
	_memoryChunks.clear();

	_kmerIndex.clear();
	_kmerIndex.reserve(0);

	_kmerCounts.clear();
	_kmerCounts.reserve(0);
}
