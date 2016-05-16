//(c) 2016 by Authors
//This file is a part of ABruijn program.
//Released under the BSD license (see LICENSE file)

#include <cassert>
#include <stdexcept>
#include <iostream>
#include <unordered_set>
#include <algorithm>

#include "kmer_index.h"
#include "utility.h"

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
	LOG_PRINT("Hard threshold set to " << hardThreshold);
	if (hardThreshold == 0 || hardThreshold > 100) 
	{
		throw std::runtime_error("Wrong hard threshold value: " + 
								 std::to_string(hardThreshold));
	}
	static const size_t BLOOM_CELLS = (size_t)1 << 31;
	static const size_t BLOOM_HASH = 5;
	CountingBloom bloomFilter(hardThreshold, BLOOM_HASH, BLOOM_CELLS);

	//filling up bloom filter
	DEBUG_PRINT("First pass:");
	size_t counterDone = 0;
	int prevPercent = -1;
	for (auto& seqPair : seqContainer.getIndex())
	{
		++counterDone;
		if (seqPair.second.sequence.length() < _kmerSize) 
			continue;

		int percent = 10 * counterDone / seqContainer.getIndex().size();
		if (percent > prevPercent)
		{
			DEBUG_PRINT(percent * 10 << "%");
			prevPercent = percent;
		}

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
	DEBUG_PRINT("Second pass:");
	counterDone = 0;
	prevPercent = -1;
	for (auto& seqPair : seqContainer.getIndex())
	{
		++counterDone;
		if (seqPair.second.sequence.length() < _kmerSize) 
			continue;

		int percent = 10 * counterDone / seqContainer.getIndex().size();
		if (percent > prevPercent)
		{
			DEBUG_PRINT(percent * 10 << "% ");
			prevPercent = percent;
		}

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
	DEBUG_PRINT("Initial size: " << _kmerIndex.size());
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
			++itKmers;
		}
	}
	DEBUG_PRINT("Removed " << removedCount << " entries");
}

void VertexIndex::buildReadIndex()
{
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
	}
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
