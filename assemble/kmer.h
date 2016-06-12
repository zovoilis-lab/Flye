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

class KmerIterator
{
public:
	typedef typename std::allocator<Kmer>::difference_type difference_type;
	typedef typename std::allocator<Kmer>::value_type value_type;
	typedef typename std::allocator<Kmer>::reference reference;
    typedef typename std::allocator<Kmer>::pointer pointer;
    typedef std::forward_iterator_tag iterator_category;

	KmerIterator(const std::string* readSeq, size_t position = 0);
	KmerIterator(const KmerIterator& other);
	KmerIterator& operator=(const KmerIterator& other);

    bool operator==(const KmerIterator&) const;
    bool operator!=(const KmerIterator&) const;

	value_type operator*() const;
	KmerIterator& operator++();

private:
	const std::string* 	_readSeq;
	size_t 				_position;
	Kmer 				_kmer;
};

class ReadKmersWrapper
{
public:
	ReadKmersWrapper(FastaRecord::Id readId):
		_readId(readId)
	{}

	KmerIterator begin();
	KmerIterator end();

private:
	FastaRecord::Id _readId;
};
