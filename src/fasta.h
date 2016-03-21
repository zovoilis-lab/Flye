//****************************************************************************
//* Copyright (c) 2012 Saint-Petersburg Academic University
//* All Rights Reserved
//* See file LICENSE for details.
//****************************************************************************

#pragma once

#include <sstream>
#include <fstream>
#include <iterator>
#include <vector>

struct FastaRecord
{
	typedef int32_t ReadIdType;

	FastaRecord() {}
	FastaRecord(const std::string& sequence, const std::string& description, 
				ReadIdType id):
		id_(id), sequence_(sequence), description_(description)
	{
	}
	
	ReadIdType id_;
	std::string sequence_;
	std::string description_;		
};

class FastaReader
{
public:
	explicit FastaReader(const std::string& fileName): 
		inputStream_(fileName.c_str()),
		fileName_(fileName) {}
	size_t 	GetSequences(std::vector<FastaRecord>& record);
	size_t 	GetSequencesWithComplements(std::vector<FastaRecord>& record);
	bool 	IsOk() const;
private:
	void ValidateSequence(std::string& sequence);
	void ValidateHeader(std::string& header);

	std::ifstream inputStream_;
	std::string fileName_;
};
