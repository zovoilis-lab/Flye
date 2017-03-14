//(c) 2016 by Authors
//This file is a part of ABruijn program.
//Released under the BSD license (see LICENSE file)

#pragma once

#include <unordered_map>
#include <memory>

#include "sequence_container.h"

static_assert(sizeof(size_t) == 8, "32-bit architectures are not supported");

class Kmer
{
public:
	Kmer(): _representation(0) {}
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
	{
		size_t x = _representation;
		size_t z = (x += 0x9E3779B97F4A7C15ULL);
		z = (z ^ (z >> 30)) * 0xBF58476D1CE4E5B9ULL;
		z = (z ^ (z >> 27)) * 0x94D049BB133111EBULL;
		return z ^ (z >> 31);
	}

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

struct KmerPosition
{
	KmerPosition(Kmer kmer, int32_t position):
		kmer(kmer), position(position) {}
	Kmer kmer;
	int32_t position;
};

class KmerIterator
{
public:
    typedef std::forward_iterator_tag iterator_category;

	KmerIterator(const std::string* readSeq, size_t position);

    bool operator==(const KmerIterator&) const;
    bool operator!=(const KmerIterator&) const;

	KmerPosition operator*() const;
	KmerIterator& operator++();

protected:
	const std::string* 	_readSeq;
	size_t 				_position;
	Kmer 				_kmer;
};

class SolidKmerIterator : public KmerIterator
{
public:
	SolidKmerIterator(const std::string* readSeq, size_t position);

    bool operator==(const SolidKmerIterator&) const;
    bool operator!=(const SolidKmerIterator&) const;

	KmerPosition operator*() const;
	SolidKmerIterator& operator++();
};

class IterKmers
{
public:
	IterKmers(FastaRecord::Id readId):
		_readId(readId)
	{}

	KmerIterator begin();
	KmerIterator end();

private:
	FastaRecord::Id _readId;
};

class IterSolidKmers
{
public:
	IterSolidKmers(FastaRecord::Id readId):
		_readId(readId)
	{}

	SolidKmerIterator begin();
	SolidKmerIterator end();

private:
	FastaRecord::Id _readId;
};
