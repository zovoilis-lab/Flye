//****************************************************************************
//* Copyright (c) 2012 Saint-Petersburg Academic University
//* All Rights Reserved
//* See file LICENSE for details.
//****************************************************************************

#include "fasta.h"
#include <stdexcept>
#include <fstream>
#include <sstream>

namespace
{
	class ParseException : public std::runtime_error 
	{
	public:
		ParseException(const std::string & what):
			std::runtime_error(what)
		{}
	};

	char complementSymbol(char c)
	{
		switch (c)
		{
			case 'A': return 'T';
			case 'C': return 'G';
			case 'G': return 'C';
			case 'T': return 'A';
			case 'a': return 't';
			case 'c': return 'g';
			case 'g': return 'c';
			case 't': return 'a';
			default: return c;
		}
	}

	void reverseComplement(const std::string& dnaIn, std::string& dnaOut)
	{
		dnaOut.clear();
		for (size_t i = dnaIn.length(); i > 0; --i)
		{
			dnaOut += complementSymbol(dnaIn[i - 1]);
		}
	}
}

void SequenceContainer::readFasta(const std::string& fileName)
{
	std::vector<FastaRecord> records;
	this->getSequencesWithComplements(records, fileName);
	for (auto& rec : records)
	{
		_seqIndex[rec.id] = rec;
	}
}

size_t SequenceContainer::getSequences(std::vector<FastaRecord>& record, 
									   const std::string& fileName)
{
	std::string buffer;
	std::string sequence;
	std::string header; 
	std::ifstream inputStream(fileName);
	int line = 1;

	record.clear();
	FastaRecord::ReadIdType seqId = 1;

	try
	{
		while(!inputStream.eof())
		{
			std::getline(inputStream, buffer, '\n');
			if (*buffer.rbegin() == '\r') buffer.erase(buffer.size() - 1);
			if (buffer.empty()) continue;

			if (buffer[0] == '>')
			{
				if (!header.empty())
				{
					if (sequence.empty()) throw ParseException("empty sequence");

					record.push_back(FastaRecord(sequence, header, seqId));
					++seqId;
					sequence.clear();
					header.clear();
				}
				this->validateHeader(buffer);
				header = buffer;
			}
			else
			{
				this->validateSequence(buffer);
				sequence += buffer;
			}

			++line;
		}
		
		if (sequence.empty()) throw ParseException("empty sequence");
		record.push_back(FastaRecord(sequence, header, seqId));
	}
	catch (ParseException & e)
	{
		std::stringstream ss;
		ss << "parse error in " << fileName << " on line " << line << ": " << e.what();
		throw std::runtime_error(ss.str());
	}

	return record.size();
}

size_t SequenceContainer::getSequencesWithComplements(std::vector<FastaRecord>& records,
													  const std::string& fileName)
{
	this->getSequences(records, fileName);
	std::vector<FastaRecord> complements;
	for (auto &record : records)
	{
		std::string header = "-" + record.description;
		record.description = "+" + record.description;
		std::string revComplement;
		reverseComplement(record.sequence, revComplement);
		complements.push_back(FastaRecord(revComplement, header, -record.id));
	}
	std::copy(complements.begin(), complements.end(), std::back_inserter(records));
	return records.size();
}

void SequenceContainer::validateHeader(std::string& header)
{
	size_t delim = header.find(' ');
	if (delim == std::string::npos)
	{
		delim = header.length() - 1;
	}
	else
	{
		--delim;
	}

	header = header.substr(1, delim);
	if (header.empty()) throw ParseException("empty header");
}

void SequenceContainer::validateSequence(const std::string& sequence)
{
	const std::string VALID_CHARS = "ACGTURYKMSWBDHWNX-";
	for (size_t i = 0; i < sequence.length(); ++i)
	{
		if (VALID_CHARS.find(toupper(sequence[i])) == std::string::npos) 
		{
			throw ParseException((std::string("illegal character: ") + 
								  sequence[i]).c_str());
		}
	}
}

void SequenceContainer::writeFasta(const std::vector<FastaRecord>& records, 
								   const std::string& fileName)
{
	static const size_t FASTA_SLICE = 80;
	std::ofstream fout(fileName);
	for (const FastaRecord& rec : records)
	{
		fout << ">" << rec.description << std::endl;
		for (size_t c = 0; c < rec.sequence.length(); c += FASTA_SLICE)
			fout << rec.sequence.substr(c, FASTA_SLICE) << std::endl;
	}
}
