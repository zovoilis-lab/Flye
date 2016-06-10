//(c) 2016 by Authors
//This file is a part of ABruijn program.
//Released under the BSD license (see LICENSE file)

#include <cassert>
#include <stdexcept>
#include <iostream>
#include <unordered_set>
#include <algorithm>

#include "kmer_index.h"
#include "logger.h"

namespace
{
	static std::vector<size_t> table;
	size_t dnaToId(char c)
	{
		return table[(size_t)c];
	}
	struct TableFiller
	{
		TableFiller()
		{
			static bool tableFilled = false;
			if (!tableFilled)
			{
				tableFilled = true;
				table.assign(256, -1);	//256 chars
				table[(size_t)'A'] = 0;
				table[(size_t)'a'] = 0;
				table[(size_t)'C'] = 1;
				table[(size_t)'c'] = 1;
				table[(size_t)'T'] = 2;
				table[(size_t)'t'] = 2;
				table[(size_t)'G'] = 3;
				table[(size_t)'g'] = 3;
				table[(size_t)'-'] = 4;
			}
		}
	};
	TableFiller filler;
	char idToDna(unsigned int num)
	{
		static const char LETTERS[] = "ACTG-";
		if (num > 4) throw std::runtime_error("Error converting number to DNA");
		return LETTERS[num];
	}
}

Kmer::Kmer(const std::string& dnaString):
	_representation(0)
{
	if (dnaString.length() != VertexIndex::getInstance().getKmerSize())
	{
		throw std::runtime_error("Kmer length inconsistency");
	}

	for (auto dnaChar : dnaString)	
	{
		_representation <<= 2;
		_representation += dnaToId(dnaChar);
	}
}

void Kmer::appendRight(char dnaSymbol)
{
	_representation <<= 2;
	_representation += dnaToId(dnaSymbol);

	KmerRepr kmerSize = VertexIndex::getInstance().getKmerSize();
	KmerRepr kmerMask = ((KmerRepr)1 << kmerSize * 2) - 1;
	_representation &= kmerMask;
}

void Kmer::appendLeft(char dnaSymbol)
{
	_representation >>= 2;

	KmerRepr kmerSize = VertexIndex::getInstance().getKmerSize();
	KmerRepr shift = kmerSize * 2 - 2;
	_representation += dnaToId(dnaSymbol) << shift;
}

void Kmer::reverseComplement()
{
	KmerRepr tmpRepr = _representation;
	_representation = 0;
	KmerRepr mask = 3;

	for (unsigned int i = 0; i < VertexIndex::getInstance().getKmerSize(); ++i)
	{
		_representation <<= 2;
		_representation += ~(mask & tmpRepr);
		tmpRepr >>= 2;
	}
}

std::string Kmer::dnaRepresentation() const
{
	std::string repr;
	KmerRepr mask = 3;
	KmerRepr tempRepr = _representation;
	for (unsigned int i = 0; i < VertexIndex::getInstance().getKmerSize(); ++i)
	{
		repr.push_back(idToDna(tempRepr & mask));
		tempRepr >>= 2;
	}
	repr.reserve();
	return repr;
}

VertexIndex::VertexIndex():
	_kmerSize(0)
{
}

void VertexIndex::setKmerSize(unsigned int size)
{
	_kmerSize = size;
}

void VertexIndex::buildKmerIndex(const SequenceContainer& seqContainer,
								 size_t hardThreshold)
{
	Logger::get().debug() << "Hard threshold set to " << hardThreshold;
	if (hardThreshold == 0 || hardThreshold > 100) 
	{
		throw std::runtime_error("Wrong hard threshold value: " + 
								 std::to_string(hardThreshold));
	}
	static const size_t BLOOM_CELLS = (size_t)1 << 31;
	static const size_t BLOOM_HASH = 4;
	CountingBloom bloomFilter(hardThreshold, BLOOM_HASH, BLOOM_CELLS);

	//filling up bloom filter
	Logger::get().info() << "Indexing kmers (1/2):";
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
			bloomFilter.add(curKmer.hash());
			if (pos == seqPair.second.sequence.length()) break;
			curKmer.appendRight(seqPair.second.sequence[pos++]);
		}
	}

	//adding kmers that have passed the filter
	Logger::get().info() << "Indexing kmers (2/2):";
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
			size_t count = bloomFilter.count(curKmer.hash(), hardThreshold - 1);
			if (count >= hardThreshold)
			{
				_kmerIndex[curKmer]
					.push_back(ReadPosition(seqPair.second.id, pos));
			}
			else
			{
				_kmerDistribution[count] += 1;
			}

			if (pos == seqPair.second.sequence.length()) break;
			curKmer.appendRight(seqPair.second.sequence[pos++]);
		}
	}
	
	for (auto kmer : _kmerIndex)
	{
		_kmerDistribution[kmer.second.size()] += 1;
	}
}


void VertexIndex::applyKmerThresholds(unsigned int minCoverage, 
									  unsigned int maxCoverage)
{
	int removedCount = 0;
	Logger::get().debug() << "Initial size: " << _kmerIndex.size();
	for (auto itKmers = _kmerIndex.begin(); itKmers != _kmerIndex.end();)
	{
		if (itKmers->second.size() < minCoverage || 
			itKmers->second.size() > maxCoverage)
		{
			itKmers = _kmerIndex.erase(itKmers);
			++removedCount;
		}
		else
		{
			itKmers->second.shrink_to_fit();
			++itKmers;
		}
	}
	Logger::get().debug() << "Removed " << removedCount << " entries";
}

void VertexIndex::buildReadIndex()
{
	Logger::get().info() << "Building read index";
	for (auto& kmerHash: _kmerIndex)
	{
		for (auto& kmerPosPair : kmerHash.second)
		{
			_readIndex[kmerPosPair.readId]
				.push_back(KmerPosition(kmerHash.first, kmerPosPair.position));
		}
	}
	for (auto& readHash : _readIndex)
	{
		std::sort(readHash.second.begin(), readHash.second.end(),
					[](const KmerPosition& p1, const KmerPosition& p2)
						{return p1.position < p2.position;});
		readHash.second.shrink_to_fit();
	}

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

void VertexIndex::outputCounts() const
{
	for (auto& hashPair : _kmerIndex)
	{
		{
			std::cout << hashPair.first.dnaRepresentation() << "\t" 
					  << hashPair.second.size() << std::endl;
		}
	}
}
