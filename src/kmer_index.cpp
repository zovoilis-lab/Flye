#include "kmer_index.h"
#include <cassert>
#include <stdexcept>
#include <iostream>

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
	if (dnaString.length() != KmerIndex::getIndex().getKmerSize())
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

	KmerRepr kmerSize = KmerIndex::getIndex().getKmerSize();
	KmerRepr kmerMask = ((KmerRepr)1 << kmerSize * 2) - 1;
	newKmer._representation &= kmerMask;

	return newKmer;
}

Kmer Kmer::appendLeft(char dnaSymbol) const
{
	Kmer newKmer(*this);
	newKmer._representation >>= 2;

	KmerRepr kmerSize = KmerIndex::getIndex().getKmerSize();
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

	for (unsigned int i = 0; i < KmerIndex::getIndex().getKmerSize(); ++i)
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
	for (unsigned int i = 0; i < KmerIndex::getIndex().getKmerSize(); ++i)
	{
		repr.push_back(numberToDna(tempRepr & mask));
		tempRepr >>= 2;
	}
	repr.reserve();
	return repr;
}

KmerIndex::KmerIndex():
	_kmerSize(0)
{}

void KmerIndex::setKmerSize(unsigned int size)
{
	_kmerSize = size;
}

void KmerIndex::addFastaSequence(const FastaRecord& fastaRecord)
{
	if (fastaRecord.sequence_.length() < _kmerSize) return;

	Kmer curKmer(fastaRecord.sequence_.substr(0, _kmerSize));
	_kmerCount[curKmer] += 1;
	_kmerCount[curKmer.reverseComplement()] += 1;
	for (size_t pos = _kmerSize; pos < fastaRecord.sequence_.length(); ++pos)
	{
		curKmer = curKmer.appendRight(fastaRecord.sequence_[pos]);
		_kmerCount[curKmer] += 1;
		_kmerCount[curKmer.reverseComplement()] += 1;
	}
}

void KmerIndex::outputCounts()
{
	for (auto hashPair : _kmerCount)
	{
		//if (hashPair.second > 0)
		{
			std::cout << hashPair.first.dnaRepresentation() << "\t" 
					  << hashPair.second << std::endl;
		}
	}
}
