//(c) 2016 by Authors
//This file is a part of ABruijn program.
//Released under the BSD license (see LICENSE file)

#include <stdexcept>

#include "kmer.h"
#include "vertex_index.h"

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

