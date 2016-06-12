//(c) 2016 by Authors
//This file is a part of ABruijn program.
//Released under the BSD license (see LICENSE file)

#pragma once

#include <unordered_map>

#include "sequence_container.h"

static_assert(sizeof(size_t) == 8, "32-bit architectures are not supported");

class Kmer
{
public:
	Kmer(const std::string& dnaString);

	void reverseComplement();
	void appendRight(char dnaSymbol);
	void appendLeft(char dnaSymbol);
	std::string dnaRepresentation() const;
	typedef size_t KmerRepr;

	bool operator == (const Kmer& other) const
		{return this->_representation == other._representation;}
	bool operator != (const Kmer& other) const
		{return !(*this == other);}
	size_t hash() const
		{return _representation * 0x9ddfea08eb382d69ULL;}

private:
	KmerRepr _representation;
};

namespace std
{
	template <>
	struct hash<Kmer>
	{
		std::size_t operator()(const Kmer& kmer) const
		{
			return kmer.hash();
		}
	};
}

