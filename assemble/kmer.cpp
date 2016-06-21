//(c) 2016 by Authors
//This file is a part of ABruijn program.
//Released under the BSD license (see LICENSE file)

#include <stdexcept>
#include <cassert>

#include "kmer.h"
#include "vertex_index.h"
#include "sequence_container.h"

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
	if (dnaString.length() != VertexIndex::get().getKmerSize())
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

	KmerRepr kmerSize = VertexIndex::get().getKmerSize();
	KmerRepr kmerMask = ((KmerRepr)1 << kmerSize * 2) - 1;
	_representation &= kmerMask;
}

void Kmer::appendLeft(char dnaSymbol)
{
	_representation >>= 2;

	KmerRepr kmerSize = VertexIndex::get().getKmerSize();
	KmerRepr shift = kmerSize * 2 - 2;
	_representation += dnaToId(dnaSymbol) << shift;
}

void Kmer::reverseComplement()
{
	KmerRepr tmpRepr = _representation;
	_representation = 0;
	KmerRepr mask = 3;

	for (unsigned int i = 0; i < VertexIndex::get().getKmerSize(); ++i)
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
	for (unsigned int i = 0; i < VertexIndex::get().getKmerSize(); ++i)
	{
		repr.push_back(idToDna(tempRepr & mask));
		tempRepr >>= 2;
	}
	repr.reserve();
	return repr;
}

KmerIterator::KmerIterator(const std::string* readSeq, size_t position):
	_readSeq(readSeq),
	_position(position)
{
	if (position != readSeq->length() - VertexIndex::get().getKmerSize())
	{
		_kmer = Kmer(readSeq->substr(0, VertexIndex::get().getKmerSize()));
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
	size_t appendPos = _position + VertexIndex::get().getKmerSize();
	_kmer.appendRight((*_readSeq)[appendPos]);
	++_position;
	return *this;
}

KmerPosition KmerIterator::operator*() const
{
	return KmerPosition(_kmer, _position);
}

KmerPosition SolidKmerIterator::operator*() const
{
	assert(VertexIndex::get().isSolid(_kmer));
	return KmerPosition(_kmer, _position);
}

SolidKmerIterator::SolidKmerIterator(const std::string* readSeq, 
									 size_t position):
	KmerIterator(readSeq, position)
{
	if (!VertexIndex::get().isSolid(_kmer) && _position == 0) ++(*this);
}


SolidKmerIterator& SolidKmerIterator::operator++()
{
	size_t appendPos = _position + VertexIndex::get().getKmerSize();
	do
	{
		_kmer.appendRight((*_readSeq)[appendPos++]);
	}
	while(!VertexIndex::get().isSolid(_kmer) &&
		  appendPos < _readSeq->length());

	_position = appendPos - VertexIndex::get().getKmerSize();
	return *this;
}

bool SolidKmerIterator::operator==(const SolidKmerIterator& other) const
{
	return _readSeq == other._readSeq && 
		   _position == other._position;
}

bool SolidKmerIterator::operator!=(const SolidKmerIterator& other) const
{
	return !(*this == other);
}

KmerIterator IterKmers::begin()
{
	const std::string& seq = SequenceContainer::get().getIndex()
										   .at(_readId).sequence;
	if (seq.length() < VertexIndex::get().getKmerSize()) 
		return this->end();

	return KmerIterator(&seq, 0);
}

KmerIterator IterKmers::end()
{
	const std::string& seq = SequenceContainer::get().getIndex()
										   .at(_readId).sequence;
	return KmerIterator(&seq, seq.length() - 
						VertexIndex::get().getKmerSize());
}

SolidKmerIterator IterSolidKmers::begin()
{
	const std::string& seq = SequenceContainer::get().getIndex()
										   .at(_readId).sequence;
	if (seq.length() < VertexIndex::get().getKmerSize()) 
		return this->end();

	return SolidKmerIterator(&seq, 0);
}

SolidKmerIterator IterSolidKmers::end()
{

	const std::string& seq = SequenceContainer::get().getIndex()
										   .at(_readId).sequence;

	return SolidKmerIterator(&seq, seq.length() - 
							 VertexIndex::get().getKmerSize());
}
