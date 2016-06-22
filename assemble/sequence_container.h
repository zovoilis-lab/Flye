//(c) 2016 by Authors
//This file is a part of ABruijn program.
//Released under the BSD license (see LICENSE file)

#pragma once

#include <vector>
#include <unordered_map>
#include <string>
#include <limits>

struct FastaRecord
{
	struct Id
	{
		explicit Id(uint32_t id): _id(id) {}

		bool operator==(const Id& other) const
			{return _id == other._id;}
		bool operator!=(const Id& other) const
			{return !(*this == other);}

		Id rc() const		//reverse complement 
			{return Id(_id + 1 - (_id % 2) * 2);}
		bool strand() const		//true = positive, false = negative
			{return _id % 2;}
		size_t hash() const 
			{return 0x9ddfea08eb382d69ULL * (size_t)_id;}
	private:
		uint32_t _id;
	};
	static const Id ID_NONE; 
	typedef std::tuple<Id, Id> IdPair;

	FastaRecord(): id(ID_NONE) {}
	FastaRecord(const std::string& sequence, const std::string& description, 
				Id id):
		id(id), sequence(sequence), description(description)
	{
	}
	
	Id id;
	std::string sequence;
	std::string description;		
};

namespace std
{
	template <>
	struct hash<FastaRecord::Id> 
	{
		size_t operator() (const FastaRecord::Id& h) const throw() 
		{
			 return h.hash();
		}
	};

	template <>
	struct hash<FastaRecord::IdPair> 
	{
		 size_t operator()(const FastaRecord::IdPair& k) const
		 {
			size_t lhs = std::get<0>(k).hash();
			size_t rhs = std::get<1>(k).hash();
			lhs ^= rhs + 0x9ddfea08eb382d69ULL + (lhs << 6) + (lhs >> 2);
			return lhs;
		 }
	};
}


class SequenceContainer
{
public:
	typedef std::unordered_map<FastaRecord::Id, 
							   FastaRecord> SequenceIndex;

	static SequenceContainer& get()
	{
		static SequenceContainer container;
		return container;
	}
	static void writeFasta(const std::vector<FastaRecord>& records,
						   const std::string& fileName);

	const SequenceIndex& getIndex() const
		{return _seqIndex;}
	size_t seqLen(FastaRecord::Id readId) const
		{return _seqIndex.at(readId).sequence.length();}
	std::string seqName(FastaRecord::Id readId) const
		{return _seqIndex.at(readId).description;}
	void 	readFasta(const std::string& filename);

private:
	SequenceContainer() {}

	size_t 	getSequences(std::vector<FastaRecord>& record, 
						 const std::string& fileName);
	size_t 	getSequencesWithComplements(std::vector<FastaRecord>& record, 
										const std::string& fileName);
	void 	validateSequence(const std::string& sequence);
	void 	validateHeader(std::string& header);

	SequenceIndex _seqIndex;
};
