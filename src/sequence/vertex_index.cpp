//(c) 2016 by Authors
//This file is a part of ABruijn program.
//Released under the BSD license (see LICENSE file)

#include <stdexcept>
#include <iostream>
#include <unordered_set>
#include <algorithm>
#include <queue>

#include "vertex_index.h"
#include "../common/logger.h"
#include "../common/parallel.h"
#include "../common/config.h"


void VertexIndex::countKmers(size_t hardThreshold, int genomeSize)
{
	if (Parameters::get().kmerSize > 31)
	{
		throw std::runtime_error("Maximum k-mer size is 31");
	}

	Logger::get().debug() << "Hard threshold set to " << hardThreshold;
	if (hardThreshold == 0)
	{
		throw std::runtime_error("Wrong hard threshold value: " + 
								 std::to_string(hardThreshold));
	}
	Logger::get().debug() << "Started k-mer counting";

	size_t preCountSize = 1024 * 1024 * 1024;	//1G by default
	if (genomeSize > (int)Config::get("big_genome_threshold"))
	{
		preCountSize *= 4 * 4;					//16G in case of larger genomes
	}
	auto preCounters = new std::atomic<unsigned char>[preCountSize];
	for (size_t i = 0; i < preCountSize; ++i) preCounters[i] = 0;

	std::vector<FastaRecord::Id> allReads;
	for (const auto& seq : _seqContainer.iterSeqs())
	{
		allReads.push_back(seq.id);
	}

	//first pass: filling up naive hash counting filter
	if (_outputProgress) Logger::get().info() << "Counting k-mers (1/2):";
	std::function<void(const FastaRecord::Id&)> preCountUpdate = 
	[&preCounters, this, preCountSize] 
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
	if (_outputProgress) Logger::get().info() << "Counting k-mers (2/2):";

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
	
	for (const auto& kmer : _kmerCounts.lock_table())
	{
		_kmerDistribution[kmer.second] += 1;
		_repetitiveFrequency = std::max(_repetitiveFrequency, kmer.second);
	}
	delete[] preCounters;
}

void VertexIndex::buildIndexUnevenCoverage(int globalMinFreq, float selectRate,
										   int tandemFreq)
{
	this->setRepeatCutoff(globalMinFreq);

	//_solidMultiplier = 1;

	std::vector<FastaRecord::Id> allReads;
	for (const auto& seq : _seqContainer.iterSeqs())
	{
		allReads.push_back(seq.id);
	}

	//first, count the number of k-mers that will be actually stored in the index
	_kmerIndex.reserve(_kmerCounts.size() / 10);
	if (_outputProgress) Logger::get().info() << "Filling index table (1/2)";
	std::function<void(const FastaRecord::Id&)> initializeIndex = 
	[this, globalMinFreq, selectRate, tandemFreq] (const FastaRecord::Id& readId)
	{
		if (!readId.strand()) return;

		auto topKmers = this->yieldFrequentKmers(readId, selectRate, tandemFreq);
		for (auto kmerFreq : topKmers)
		{
			kmerFreq.kmer.standardForm();

			if (kmerFreq.freq < (size_t)globalMinFreq ||
				kmerFreq.freq > _repetitiveFrequency) continue;

			ReadVector defVec((uint32_t)1, (uint32_t)0);
			_kmerIndex.upsert(kmerFreq.kmer, 
							  [](ReadVector& rv){++rv.capacity;}, defVec);
		}
	};
	processInParallel(allReads, initializeIndex, 
					  Parameters::get().numThreads, _outputProgress);
	
	this->allocateIndexMemory();

	if (_outputProgress) Logger::get().info() << "Filling index table (2/2)";
	std::function<void(const FastaRecord::Id&)> indexUpdate = 
	[this, globalMinFreq, selectRate, tandemFreq] (const FastaRecord::Id& readId)
	{
		if (!readId.strand()) return;

		auto topKmers = this->yieldFrequentKmers(readId, selectRate, tandemFreq);
		for (auto kmerFreq : topKmers)
		{
			if (kmerFreq.freq < (size_t)globalMinFreq ||
				kmerFreq.freq > _repetitiveFrequency) continue;

			KmerPosition kmerPos(kmerFreq.kmer, kmerFreq.position);
			FastaRecord::Id targetRead = readId;
			bool revCmp = kmerPos.kmer.standardForm();
			if (revCmp)
			{
				kmerPos.position = _seqContainer.seqLen(readId) - 
										kmerPos.position -
										Parameters::get().kmerSize;
				targetRead = targetRead.rc();
			}

			//will not trigger update for k-mer not in the index
			_kmerIndex.update_fn(kmerPos.kmer, 
				[targetRead, &kmerPos, this](ReadVector& rv)
				{
					if (rv.size == rv.capacity) 
					{
						Logger::get().warning() << "Index size mismatch " << rv.capacity;
						return;
					}
					size_t globPos = _seqContainer
							.globalPosition(targetRead, kmerPos.position);
					rv.data[rv.size].set(globPos);
					++rv.size;
				});
		}
	};
	processInParallel(allReads, indexUpdate, 
					  Parameters::get().numThreads, _outputProgress);

	_kmerCounts.clear();
	_kmerCounts.reserve(0);

	Logger::get().debug() << "Sorting k-mer index";
	for (const auto& kmerVec : _kmerIndex.lock_table())
	{
		std::sort(kmerVec.second.data, kmerVec.second.data + kmerVec.second.size,
				  [](const IndexChunk& p1, const IndexChunk& p2)
				  	{return p1.get() < p2.get();});
	}
	
	size_t totalEntries = 0;
	for (const auto& kmerRec : _kmerIndex.lock_table())
	{
		totalEntries += kmerRec.second.size;
	}
	Logger::get().debug() << "Selected k-mers: " << _kmerIndex.size();
	Logger::get().debug() << "Index size: " << totalEntries;
	Logger::get().debug() << "Mean k-mer frequency: " 
		<< (float)totalEntries / _kmerIndex.size();
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

//TODO: switch to the markRepeats function below
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
	Logger::get().debug() << "Repetitive k-mer frequency: " 
						  << _repetitiveFrequency;
	Logger::get().debug() << "Filtered " << repetitiveKmers 
						  << " repetitive k-mers (" <<
						  filteredRate << ")";
}

void VertexIndex::filterFrequentKmers(float rate)
{
	size_t totalKmers = 0;
	size_t uniqueKmers = _kmerIndex.size();
	for (const auto& kmer : _kmerIndex.lock_table())
	{
		totalKmers += kmer.second.capacity;
	}
	float meanFrequency = (float)totalKmers / (uniqueKmers + 1);
	_repetitiveFrequency = rate * meanFrequency;
	
	size_t repetitiveKmers = 0;
	for (const auto& kmer : _kmerIndex.lock_table())
	{
		if (kmer.second.capacity > _repetitiveFrequency)
		{
			++repetitiveKmers;
			//repetitiveKmers += kmer.second.capacity;
			_repetitiveKmers.insert(kmer.first, true);
		}
	}

	for (const auto& kmer : _repetitiveKmers.lock_table())
	{
		_kmerIndex.erase(kmer.first);
	}

	float filteredRate = (float)repetitiveKmers / uniqueKmers;
	Logger::get().debug() << "Repetitive k-mer frequency: " 
						  << _repetitiveFrequency;
	Logger::get().debug() << "Filtered " << repetitiveKmers 
						  << " repetitive k-mers (" <<
						  filteredRate << ")";
}

void VertexIndex::buildIndex(int minCoverage)
{
	this->setRepeatCutoff(minCoverage);

	if (_outputProgress) Logger::get().info() << "Filling index table";
	//_solidMultiplier = 1;
	
	//"Replacing" k-mer couns with k-mer index. We need multiple passes
	//to avoid peaks in memory usage during the hash table extensions +
	//prevent memory fragmentation
	
	size_t kmerEntries = 0;
	size_t solidKmers = 0;
	for (const auto& kmer : _kmerCounts.lock_table())
	{
		if ((size_t)minCoverage <= kmer.second &&
			kmer.second < _repetitiveFrequency)
		{
			kmerEntries += kmer.second;
			++solidKmers;
		}
		if (kmer.second > _repetitiveFrequency)
		{
			_repetitiveKmers.insert(kmer.first, true);
		}
	}
	
	_kmerIndex.reserve(solidKmers);
	for (const auto& kmer : _kmerCounts.lock_table())
	{
		if ((size_t)minCoverage <= kmer.second &&
			kmer.second < _repetitiveFrequency)
		{
			ReadVector rv((uint32_t)kmer.second, 0);
			_kmerIndex.insert(kmer.first, rv);
		}
	}
	_kmerCounts.clear();
	_kmerCounts.reserve(0);

	this->allocateIndexMemory();

	std::function<void(const FastaRecord::Id&)> indexUpdate = 
	[this] (const FastaRecord::Id& readId)
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

			FastaRecord::Id targetRead = readId;
			bool revCmp = kmerPos.kmer.standardForm();
			if (revCmp)
			{
				kmerPos.position = _seqContainer.seqLen(readId) - 
										kmerPos.position -
										Parameters::get().kmerSize;
				targetRead = targetRead.rc();
			}
			
			//will not trigger update if the k-mer is not already in index
			_kmerIndex.update_fn(kmerPos.kmer, 
				[targetRead, &kmerPos, this](ReadVector& rv)
				{
					size_t globPos = _seqContainer
							.globalPosition(targetRead, kmerPos.position);
					//if (globPos > MAX_INDEX) throw std::runtime_error("Too much!");
					rv.data[rv.size].set(globPos);
					++rv.size;
				});
		}
	};
	std::vector<FastaRecord::Id> allReads;
	for (const auto& seq : _seqContainer.iterSeqs())
	{
		allReads.push_back(seq.id);
	}
	processInParallel(allReads, indexUpdate, 
					  Parameters::get().numThreads, _outputProgress);

	Logger::get().debug() << "Sorting k-mer index";
	for (const auto& kmerVec : _kmerIndex.lock_table())
	{
		std::sort(kmerVec.second.data, kmerVec.second.data + kmerVec.second.size,
				  [](const IndexChunk& p1, const IndexChunk& p2)
				  	{return p1.get() < p2.get();});
	}

	//Logger::get().debug() << "Sampling rate: " << _sampleRate;
	Logger::get().debug() << "Selected k-mers: " << solidKmers;
	Logger::get().debug() << "K-mer index size: " << kmerEntries;
	Logger::get().debug() << "Mean k-mer frequency: " 
		<< (float)kmerEntries / solidKmers;
}

std::vector<VertexIndex::KmerFreq>
	VertexIndex::yieldFrequentKmers(const FastaRecord::Id& seqId,
									float selectRate, int tandemFreq)
{
	thread_local std::unordered_map<Kmer, size_t> localFreq;
	localFreq.clear();
	std::vector<KmerFreq> topKmers;
	topKmers.reserve(_seqContainer.seqLen(seqId));

	for (const auto& kmerPos : IterKmers(_seqContainer.getSeq(seqId)))
	{
		auto stdKmer = kmerPos.kmer;
		stdKmer.standardForm();
		size_t freq = 1;
		_kmerCounts.find(stdKmer, freq);

		++localFreq[stdKmer];
		topKmers.push_back({kmerPos.kmer, kmerPos.position, freq});
	}

	if (topKmers.empty()) return {};
	std::sort(topKmers.begin(), topKmers.end(),
			  [](const KmerFreq& k1, const KmerFreq& k2)
			   {return k1.freq > k2.freq;});
	const size_t maxKmers = selectRate * topKmers.size();
	const size_t minFreq = topKmers[maxKmers].freq;

	auto itVec = topKmers.begin();
	while(itVec != topKmers.end() && itVec->freq >= minFreq) ++itVec;
	topKmers.erase(itVec, topKmers.end());

	topKmers.erase(std::remove_if(topKmers.begin(), topKmers.end(),
				   		[tandemFreq](KmerFreq kf)
						{
							kf.kmer.standardForm();
							return localFreq[kf.kmer] > (size_t)tandemFreq;
						}), 
				   topKmers.end());

	/*std::vector<KmerPosition> result;
	result.reserve(topKmers.size());
	for (auto kmerFreq : topKmers)
	{
		if (kmerFreq.freq <= tandemFreq) 
			result.push_back({kmerFreq.kmer, kmerFreq.position});
	}*/

	return topKmers;
}

std::vector<KmerPosition> 
	VertexIndex::yieldMinimizers(const FastaRecord::Id& seqId, int window)
{
	struct KmerAndHash
	{
		KmerPosition kp;
		size_t hash;
	};
	thread_local std::deque<KmerAndHash> miniQueue;
	miniQueue.clear();

	std::vector<KmerPosition> minimizers;
	minimizers.reserve(_seqContainer.seqLen(seqId) / window * 2);

	for (auto kmerPos : IterKmers(_seqContainer.getSeq(seqId)))
	{
		auto stdKmer = kmerPos.kmer;
		stdKmer.standardForm();
		size_t curHash = stdKmer.hash();
		
		while (!miniQueue.empty() && miniQueue.back().hash > curHash)
		{
			miniQueue.pop_back();
		}
		miniQueue.push_back({kmerPos, curHash});
		if (miniQueue.front().kp.position <= kmerPos.position - window)
		{
			while (miniQueue.front().kp.position <= kmerPos.position - window)
			{
				miniQueue.pop_front();
			}
			while (miniQueue.size() >= 2 && miniQueue[0].hash == miniQueue[1].hash)
			{
				miniQueue.pop_front();
			}
		}
		if (minimizers.empty() || minimizers.back().position != 
								  miniQueue.front().kp.position)
		{
			minimizers.push_back(miniQueue.front().kp);
		}
	}

	//Logger::get().debug() << _seqContainer.seqLen(seqId) << " " << minimizers.size();
	return minimizers;
}

void VertexIndex::allocateIndexMemory()
{
	_memoryChunks.push_back(new IndexChunk[MEM_CHUNK]);
	size_t chunkOffset = 0;
	//Important: since packed structures are apparently not thread-safe,
	//make sure that adjacent k-mer index arrays (that are accessed in parallel)
	//do not overlap within 8-byte window
	const size_t PADDING = 1;
	for (auto& kmer : _kmerIndex.lock_table())
	{
		if (MEM_CHUNK < kmer.second.capacity + PADDING) 
		{
			throw std::runtime_error("k-mer is too frequent");
		}
		if (MEM_CHUNK - chunkOffset < kmer.second.capacity + PADDING)
		{
			_memoryChunks.push_back(new IndexChunk[MEM_CHUNK]);
			chunkOffset = 0;
		}
		kmer.second.data = _memoryChunks.back() + chunkOffset;
		chunkOffset += kmer.second.capacity + PADDING;
	}

	//Logger::get().debug() << "Total chunks " << _memoryChunks.size()
	//	<< " wasted space: " << wasted;
}

void VertexIndex::buildIndexMinimizers(int minCoverage, int wndLen)
{
	if (_outputProgress) Logger::get().info() << "Building minimizer index";

	std::vector<FastaRecord::Id> allReads;
	size_t totalLen = 0;
	for (const auto& seq : _seqContainer.iterSeqs())
	{
		allReads.push_back(seq.id);
		if (seq.id.strand()) totalLen += seq.sequence.length();
	}

	_kmerIndex.reserve(1000000);
	if (_outputProgress) Logger::get().info() << "Pre-calculating index storage";
	std::function<void(const FastaRecord::Id&)> initializeIndex = 
	[this, minCoverage, wndLen] (const FastaRecord::Id& readId)
	{
		if (!readId.strand()) return;

		auto minimizers = this->yieldMinimizers(readId, wndLen);
		for (auto kmerPos : minimizers)
		{
			auto stdKmer = kmerPos.kmer;
			stdKmer.standardForm();
			ReadVector defVec((uint32_t)1, (uint32_t)0);
			_kmerIndex.upsert(stdKmer, 
							  [](ReadVector& rv){++rv.capacity;}, defVec);
		}
	};
	processInParallel(allReads, initializeIndex, 
					  Parameters::get().numThreads, _outputProgress);

	this->filterFrequentKmers((float)Config::get("repeat_kmer_rate"));
	this->allocateIndexMemory();
	
	if (_outputProgress) Logger::get().info() << "Filling index";
	std::function<void(const FastaRecord::Id&)> indexUpdate = 
	[this, minCoverage, wndLen] (const FastaRecord::Id& readId)
	{
		if (!readId.strand()) return;
		auto minimizers = this->yieldMinimizers(readId, wndLen);
		for (auto kmerPos : minimizers)
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

			if (_repetitiveKmers.contains(kmerPos.kmer)) continue;

			_kmerIndex.update_fn(kmerPos.kmer, 
				[targetRead, &kmerPos, this](ReadVector& rv)
				{
					if (rv.size == rv.capacity) 
					{
						Logger::get().warning() << "Index size mismatch " << rv.capacity;
						return;
					}
					size_t globPos = _seqContainer
							.globalPosition(targetRead, kmerPos.position);
					rv.data[rv.size].set(globPos);
					++rv.size;
				});
		}
	};
	processInParallel(allReads, indexUpdate, 
					  Parameters::get().numThreads, _outputProgress);

	Logger::get().debug() << "Sorting k-mer index";
	for (const auto& kmerVec : _kmerIndex.lock_table())
	{
		std::sort(kmerVec.second.data, kmerVec.second.data + kmerVec.second.size,
				  [](const IndexChunk& p1, const IndexChunk& p2)
				  	{return p1.get() < p2.get();});
	}
	
	size_t totalEntries = 0;
	for (const auto& kmerRec : _kmerIndex.lock_table())
	{
		totalEntries += kmerRec.second.size;
	}
	Logger::get().debug() << "Selected k-mers: " << _kmerIndex.size();
	Logger::get().debug() << "K-mer index size: " << totalEntries;
	Logger::get().debug() << "Mean k-mer frequency: " 
		<< (float)totalEntries / _kmerIndex.size();

	float minimizerRate = (float)totalLen / totalEntries;
	Logger::get().debug() << "Minimizer rate: " << minimizerRate;
	_sampleRate = minimizerRate;
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
