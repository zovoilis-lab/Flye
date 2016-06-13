//(c) 2016 by Authors
//This file is a part of ABruijn program.
//Released under the BSD license (see LICENSE file)

#include <stdexcept>
#include <iostream>
#include <unordered_set>
#include <algorithm>

#include "vertex_index.h"
#include "logger.h"

VertexIndex::VertexIndex():
	_kmerSize(0)
{
}

void VertexIndex::setKmerSize(unsigned int size)
{
	_kmerSize = size;
}

void VertexIndex::countKmers(const SequenceContainer& seqContainer,
							 size_t hardThreshold)
{
	Logger::get().debug() << "Hard threshold set to " << hardThreshold;
	if (hardThreshold == 0 || hardThreshold > 100) 
	{
		throw std::runtime_error("Wrong hard threshold value: " + 
								 std::to_string(hardThreshold));
	}

	//TODO: precount size should correlate with kmer size
	const size_t PRE_COUNT_SIZE = 1024 * 1024 * 1024;
	std::vector<unsigned char> preCounters(PRE_COUNT_SIZE, 0);

	//filling up bloom filter
	Logger::get().info() << "Counting kmers (1/2):";
	ProgressPercent bloomProg(seqContainer.getIndex().size());
	for (auto& seqPair : seqContainer.getIndex())
	{
		bloomProg.advance();
		for (auto kmerPos : IterKmers(seqPair.first))
		{
			if (preCounters[kmerPos.kmer.hash() % PRE_COUNT_SIZE] != 
				std::numeric_limits<unsigned char>::max())
				++preCounters[kmerPos.kmer.hash() % PRE_COUNT_SIZE];
		}
	}

	//counting only kmers that have passed the filter
	Logger::get().info() << "Counting kmers (2/2):";
	ProgressPercent indexProg(seqContainer.getIndex().size());

	for (auto& seqPair : seqContainer.getIndex())
	{
		indexProg.advance();
		for (auto kmerPos : IterKmers(seqPair.first))
		{
			size_t count = preCounters[kmerPos.kmer.hash() % PRE_COUNT_SIZE];
			if (count >= hardThreshold)
			{
				_kmerCounts.upsert(kmerPos.kmer, [](size_t& num){++num;}, 1);
			}
			else
			{
				_kmerDistribution[count] += 1;
			}
		}
	}
	
	for (auto kmer : _kmerCounts.lock_table())
	{
		_kmerDistribution[kmer.second] += 1;
	}
}


void VertexIndex::buildIndex(int minCoverage, int maxCoverage)
{
	const SequenceContainer& seqContainer = SequenceContainer::get();

	Logger::get().info() << "Building kmer index";
	ProgressPercent indexProg(seqContainer.getIndex().size());

	for (auto& seqPair : seqContainer.getIndex())
	{
		indexProg.advance();

		for (auto kmerPos : IterKmers(seqPair.first))
		{
			size_t count = 0;
			_kmerCounts.find(kmerPos.kmer, count);
			if ((size_t)minCoverage <= count && count <= (size_t)maxCoverage)
			{
				if (!_kmerIndex.contains(kmerPos.kmer))
				{
					auto ptr = new ReadVector;
					ptr->reserve(count);
					_kmerIndex.insert(kmerPos.kmer, ptr);
				}
				ReadVector* vec = _kmerIndex[kmerPos.kmer];
				vec->push_back(ReadPosition(seqPair.second.id, 
											kmerPos.position));
			}
		}
	}
	_kmerCounts.clear();
	_kmerCounts.reserve(0);

	//read indexing
	/*
	for (auto& kmerHash: _kmerIndex.lock_table())
	{
		for (auto& kmerPosPair : *kmerHash.second)
		{
			if (!_readIndex.contains(kmerPosPair.readId))
			{
				_readIndex.insert(kmerPosPair.readId, new KmerVector);
			}

			KmerVector* vec = _readIndex[kmerPosPair.readId];
			vec->push_back(KmerPosition(kmerHash.first, kmerPosPair.position));
		}
	}
	for (auto& readHash : _readIndex.lock_table())
	{
		std::sort(readHash.second->begin(), readHash.second->end(),
					[](const KmerPosition& p1, const KmerPosition& p2)
						{return p1.position < p2.position;});
		readHash.second->shrink_to_fit();
	}*/

	//distance spectrum
	/*
	const int POS_DIFF = 500;
	std::unordered_map<Kmer, std::vector<int>> distances;
	for (auto& readHash : _readIndex)
	{
		for (size_t i = 0; i < readHash.second.size() - 2; ++i)
		{
			//going right
			for (size_t j = i + 1; j < readHash.second.size(); ++j)
			{
				if (readHash.second[j].position - 
					readHash.second[i].position > POS_DIFF)
				{
					int dRight = _kmerIndex[readHash.second[j].kmer].size();
					distances[readHash.second[i].kmer].push_back(dRight);
					break;
				}
			}
			//going left
			for (int j = i; j >= 0; --j)
			{
				if (readHash.second[i].position - 
					readHash.second[j].position > POS_DIFF)
				{
					int dLeft = _kmerIndex[readHash.second[j].kmer].size();
					distances[readHash.second[i].kmer].push_back(dLeft);
					break;
				}
			}
		}
	}

	std::map<int, int> spectrum;
	for (auto& hash : distances)
	{
		double sum = 0.0f;
		for(size_t i = 0; i < hash.second.size(); i++)
		{
   			sum += hash.second[i];
		}
		double dist = !hash.second.empty() ? sum / hash.second.size() : 0.0f;
		int bin = dist;
		spectrum[bin] += 1;
	}

	for (auto& pair : spectrum)
	{
		std::cerr << pair.first << "\t" << pair.second << std::endl;
	}*/
}
