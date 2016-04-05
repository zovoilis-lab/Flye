#include <cassert>
#include <stdexcept>
#include <iostream>
#include <unordered_set>

#include <bf.h>

#include "kmer_index.h"
#include "utility.h"

namespace
{
	unsigned int dnaToNumber(char dnaChar)
	{
		switch(dnaChar)
		{
		case 'A':
			return 0;
		case 'C':
			return 1;
		case 'G':
			return 2;
		case 'T':
			return 3;
		default:
			throw std::runtime_error("Processing non-DNA character");
		}
	}

	char numberToDna(unsigned int num)
	{
		static const char LETTERS[] = "ACGT";
		if (num > 3) throw std::runtime_error("Error converting number to DNA");
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
		_representation += dnaToNumber(dnaChar);
	}
}

Kmer Kmer::appendRight(char dnaSymbol) const
{
	Kmer newKmer(*this);
	newKmer._representation <<= 2;
	newKmer._representation += dnaToNumber(dnaSymbol);

	KmerRepr kmerSize = VertexIndex::getInstance().getKmerSize();
	KmerRepr kmerMask = ((KmerRepr)1 << kmerSize * 2) - 1;
	newKmer._representation &= kmerMask;

	return newKmer;
}

Kmer Kmer::appendLeft(char dnaSymbol) const
{
	Kmer newKmer(*this);
	newKmer._representation >>= 2;

	KmerRepr kmerSize = VertexIndex::getInstance().getKmerSize();
	KmerRepr shift = kmerSize * 2 - 2;
	newKmer._representation += dnaToNumber(dnaSymbol) << shift;

	return newKmer;
}

Kmer Kmer::reverseComplement() const
{
	Kmer newKmer(*this);
	newKmer._representation = 0;
	Kmer::KmerRepr tmpRepr = _representation;
	KmerRepr mask = 3;

	for (unsigned int i = 0; i < VertexIndex::getInstance().getKmerSize(); ++i)
	{
		newKmer._representation <<= 2;
		newKmer._representation += ~(mask & tmpRepr);
		tmpRepr >>= 2;
	}
	return newKmer;
}

std::string Kmer::dnaRepresentation() const
{
	std::string repr;
	KmerRepr mask = 3;
	KmerRepr tempRepr = _representation;
	for (unsigned int i = 0; i < VertexIndex::getInstance().getKmerSize(); ++i)
	{
		repr.push_back(numberToDna(tempRepr & mask));
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

void VertexIndex::buildKmerIndex(const SequenceContainer& seqContainer)
{
	static const size_t BLOOM_CELLS = (size_t)1 << 31;
	static const size_t BLOOM_WIDTH = 2;
	static const size_t BLOOM_HASH = 7;

	bf::spectral_mi_bloom_filter bloomFilter(bf::make_hasher(BLOOM_HASH), 
											  BLOOM_CELLS, BLOOM_WIDTH);
	std::unordered_set<Kmer, Kmer::KmerHash> goodKmers;

	//filling up bloom filter
	LOG_PRINT("First pass:");
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
			std::cerr << percent * 10 << "% ";
			prevPercent = percent;
		}

		Kmer curKmer(seqPair.second.sequence.substr(0, _kmerSize));
		size_t pos = _kmerSize;
		while (true)
		{
			if (bloomFilter.lookup(curKmer.hash()) < 3)
			{
				bloomFilter.add(curKmer.hash());
			}
			else
			{
				goodKmers.insert(curKmer);
			}
				
			if (pos == seqPair.second.sequence.length())
				break;
			curKmer = curKmer.appendRight(seqPair.second.sequence[pos++]);
		}
	}
	std::cerr << std::endl;

	//adding only good kmers to index
	LOG_PRINT("Second pass:");
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
			std::cerr << percent * 10 << "% ";
			prevPercent = percent;
		}

		Kmer curKmer(seqPair.second.sequence.substr(0, _kmerSize));
		size_t pos = _kmerSize;
		while (true)
		{
			if (goodKmers.count(curKmer))
				_kmerIndex[curKmer].push_back(ReadPosition(seqPair.second.id, pos));

			if (pos == seqPair.second.sequence.length())
				break;
			curKmer = curKmer.appendRight(seqPair.second.sequence[pos++]);
		}
	}
	std::cerr << std::endl;
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
	//fighting memory fragmentation
	//for (auto& indexPair : _kmerIndex)
	//{
	//	decltype(indexPair.second) tmpVector;
	//	tmpVector.reserve(indexPair.second.size());
	//	std::copy(indexPair.second.begin(), indexPair.second.end(), 
	//			  std::back_inserter(tmpVector));
	//	indexPair.second.swap(tmpVector);
	//}
	DEBUG_PRINT(_kmerIndex.size());
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
