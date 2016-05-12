//****************************************************************************
//* Copyright (c) 2012-2016 Saint-Petersburg Academic University
//* All Rights Reserved
//* See file LICENSE for details.
//****************************************************************************

#pragma once

#include <vector>
#include <unordered_map>
#include <string>
#include <limits>

struct FastaRecord
{
	struct Id
	{
		Id(uint32_t id): _id(id) {}
		operator uint32_t () const {return _id;}
		uint32_t rc() const		//reverse complement 
			{return _id + 1 - (_id % 2) * 2;}
		size_t hash() const {return _id;}
	private:
		uint32_t _id;
	};
	static const Id ID_NONE; 

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
			 return std::hash<uint32_t>()(h.hash());
		}
	};
}

class SequenceContainer
{
public:
	typedef std::unordered_map<FastaRecord::Id, 
							   FastaRecord> SequenceIndex;

	static SequenceContainer& getInstance()
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
