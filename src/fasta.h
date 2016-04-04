//****************************************************************************
//* Copyright (c) 2012 Saint-Petersburg Academic University
//* All Rights Reserved
//* See file LICENSE for details.
//****************************************************************************

#pragma once

#include <vector>
#include <unordered_map>

struct FastaRecord
{
	typedef int32_t ReadIdType;
	static const ReadIdType ID_NONE = 0;

	FastaRecord() {}
	FastaRecord(const std::string& sequence, const std::string& description, 
				ReadIdType id):
		id(id), sequence(sequence), description(description)
	{
	}
	
	ReadIdType id;
	std::string sequence;
	std::string description;		
};

class SequenceContainer
{
public:
	typedef std::unordered_map<FastaRecord::ReadIdType, 
							   FastaRecord> SequenceIndex;

	static SequenceContainer& getInstance()
	{
		static SequenceContainer container;
		return container;
	}
	const SequenceIndex& getIndex() const
		{return _seqIndex;}
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
