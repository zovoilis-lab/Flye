#pragma once

#include <string>
#include <unordered_map>

#include "fasta.h"

class Kmer
{
public:
	Kmer(const std::string& dnaString);

	Kmer reverseComplement() const;
	Kmer appendRight(char dnaSymbol) const;
	Kmer appendLeft(char dnaSymbol) const;
	std::string dnaRepresentation() const;
	typedef uint64_t KmerRepr;

	bool operator==(const Kmer& other) const
		{return this->_representation == other._representation;}

	struct KmerHash
	{
		std::size_t operator()(const Kmer& kmer) const
		{
			return std::hash<KmerRepr>()(kmer._representation);
		}
	};

private:
	KmerRepr _representation;
};

class KmerIndex
{
public:
	static KmerIndex& getIndex()
	{
		static KmerIndex instance;
		return instance;
	}
	KmerIndex(const KmerIndex&) = delete;
	void operator=(const KmerIndex&) = delete;

	void 		 setKmerSize(unsigned int size);
	//void 		 applyLimits(unsigned int minCoverage, unsigned int maxCoverage);
	unsigned int getKmerSize() const {return _kmerSize;}
	void 		 addFastaSequence(const FastaRecord& fastaRecord);
	void 		 outputCounts();

private:
	KmerIndex();

	unsigned int _kmerSize;
	std::unordered_map<Kmer, int, Kmer::KmerHash> _kmerCount;
};
