//(c) 2016 by Authors
//This file is a part of ABruijn program.
//Released under the BSD license (see LICENSE file)

#include <stdexcept>
#include <cassert>

#include "kmer.h"
#include "vertex_index.h"
#include "sequence_container.h"
#include "../common/config.h"

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

Kmer::Kmer(const FastaRecord::DnaRepr& dnaString, 
		   size_t start, size_t length):
	_representation(0)
{
	if (length != Parameters::get().kmerSize)
	{
		throw std::runtime_error("Kmer length inconsistency");
	}

	//for (auto dnaChar : dnaString)	
	for (size_t i = start; i < start + length; ++i)	
	{
		_representation <<= 2;
		_representation += dnaToId(dnaString.at(i));
	}
}

void Kmer::appendRight(char dnaSymbol)
{
	_representation <<= 2;
	_representation += dnaToId(dnaSymbol);

	KmerRepr kmerSize = Parameters::get().kmerSize;
	KmerRepr kmerMask = ((KmerRepr)1 << kmerSize * 2) - 1;
	_representation &= kmerMask;
}

void Kmer::appendLeft(char dnaSymbol)
{
	_representation >>= 2;

	KmerRepr kmerSize = Parameters::get().kmerSize;
	KmerRepr shift = kmerSize * 2 - 2;
	_representation += dnaToId(dnaSymbol) << shift;
}

void Kmer::reverseComplement()
{
	KmerRepr tmpRepr = _representation;
	_representation = 0;
	KmerRepr mask = 3;

	for (unsigned int i = 0; i < Parameters::get().kmerSize; ++i)
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
	for (unsigned int i = 0; i < Parameters::get().kmerSize; ++i)
	{
		repr.push_back(idToDna(tempRepr & mask));
		tempRepr >>= 2;
	}
	std::reverse(repr.begin(), repr.end());
	return repr;
}

KmerIterator::KmerIterator(const FastaRecord::DnaRepr* readSeq, 
						   size_t position):
	_readSeq(readSeq),
	_position(position)
{
	if (position != readSeq->length() - Parameters::get().kmerSize)
	{
		//_kmer = Kmer(readSeq->substr(0, Parameters::get().kmerSize));
		_kmer = Kmer(*readSeq, 0, Parameters::get().kmerSize);
	}
}

bool KmerIterator::operator==(const KmerIterator& other) const
{
	return _readSeq == other._readSeq && _position == other._position;
}

bool KmerIterator::operator!=(const KmerIterator& other) const
{
	return !(*this == other);
}

KmerIterator& KmerIterator::operator++()
{
	size_t appendPos = _position + Parameters::get().kmerSize;
	_kmer.appendRight(_readSeq->at(appendPos));
	++_position;
	return *this;
}

KmerPosition KmerIterator::operator*() const
{
	return KmerPosition(_kmer, _position);
}


KmerIterator IterKmers::begin()
{
	if (_sequence.length() < Parameters::get().kmerSize) 
		return this->end();

	return KmerIterator(&_sequence, 0);
}

KmerIterator IterKmers::end()
{
	return KmerIterator(&_sequence, _sequence.length() - 
									Parameters::get().kmerSize);
}
