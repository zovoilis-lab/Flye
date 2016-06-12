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

	const size_t PRE_COUNT_SIZE = 1024 * 1024 * 1024;
	std::vector<unsigned char> preCounters(PRE_COUNT_SIZE, 0);

	//filling up bloom filter
	Logger::get().info() << "Counting kmers (1/2):";
	ProgressPercent bloomProg(seqContainer.getIndex().size());
	for (auto& seqPair : seqContainer.getIndex())
	{
		bloomProg.advance();
		if (seqPair.second.sequence.length() < _kmerSize) 
			continue;

		Kmer curKmer(seqPair.second.sequence.substr(0, _kmerSize));
		size_t pos = _kmerSize;
		while (true)
		{
			if (preCounters[curKmer.hash() % PRE_COUNT_SIZE] != 
				std::numeric_limits<unsigned char>::max())
				++preCounters[curKmer.hash() % PRE_COUNT_SIZE];

			if (pos == seqPair.second.sequence.length()) break;

			curKmer.appendRight(seqPair.second.sequence[pos++]);
		}
	}

	//adding kmers that have passed the filter
	Logger::get().info() << "Counting kmers (2/2):";
	ProgressPercent indexProg(seqContainer.getIndex().size());

	auto increaseFn = [](size_t& num) {++num;}; 
	for (auto& seqPair : seqContainer.getIndex())
	{
		indexProg.advance();
		if (seqPair.second.sequence.length() < _kmerSize) 
			continue;

		Kmer curKmer(seqPair.second.sequence.substr(0, _kmerSize));
		size_t pos = _kmerSize;
		while (true)
		{
			size_t count = preCounters[curKmer.hash() % PRE_COUNT_SIZE];
			if (count >= hardThreshold)
			{
				_kmerCounts.upsert(curKmer, increaseFn, 1);
			}
			else
			{
				_kmerDistribution[count] += 1;
			}

			if (pos == seqPair.second.sequence.length()) break;
			curKmer.appendRight(seqPair.second.sequence[pos++]);
		}
	}
	
	for (auto kmer : _kmerCounts.lock_table())
	{
		_kmerDistribution[kmer.second] += 1;
	}
}


void VertexIndex::buildIndex(const SequenceContainer& seqContainer,
							 int minCoverage, int maxCoverage)
{
	Logger::get().info() << "Building kmer index";
	ProgressPercent indexProg(seqContainer.getIndex().size());

	for (auto& seqPair : seqContainer.getIndex())
	{
		indexProg.advance();
		if (seqPair.second.sequence.length() < _kmerSize) 
			continue;

		Kmer curKmer(seqPair.second.sequence.substr(0, _kmerSize));
		size_t pos = _kmerSize;
		while (true)
		{
			size_t count = 0;
			_kmerCounts.find(curKmer, count);
			if ((size_t)minCoverage <= count && count <= (size_t)maxCoverage)
			{
				if (!_kmerIndex.contains(curKmer))
				{
					auto ptr = new ReadVector;
					ptr->reserve(count);
					_kmerIndex.insert(curKmer, ptr);
				}
				ReadVector* vec = _kmerIndex[curKmer];
				vec->push_back(ReadPosition(seqPair.second.id, pos));
			}

			if (pos == seqPair.second.sequence.length()) break;
			curKmer.appendRight(seqPair.second.sequence[pos++]);
		}
	}
	_kmerCounts.clear();
	_kmerCounts.reserve(0);

	//read indexing
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
	}

	_kmerCounts.clear();

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
