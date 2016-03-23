#include <cassert>
#include <stdexcept>
#include <iostream>

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
{}

void VertexIndex::setKmerSize(unsigned int size)
{
	_kmerSize = size;
}

void VertexIndex::addFastaSequence(const FastaRecord& fastaRecord)
{
	if (fastaRecord.sequence.length() < _kmerSize) return;

	int32_t position = 0;
	Kmer curKmer(fastaRecord.sequence.substr(0, _kmerSize));
	_kmerIndex[curKmer].push_back(ReadPosition(fastaRecord.id, position));
	for (size_t pos = _kmerSize; pos < fastaRecord.sequence.length(); ++pos)
	{
		position += 1;
		curKmer = curKmer.appendRight(fastaRecord.sequence[pos]);
		_kmerIndex[curKmer].push_back(ReadPosition(fastaRecord.id, position));
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
