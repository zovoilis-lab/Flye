//(c) 2016 by Authors
//This file is a part of ABruijn program.
//Released under the BSD license (see LICENSE file)

#include <stdexcept>
#include <fstream>
#include <sstream>
#include <iostream>
#include <random>
#include <algorithm>

#include "sequence_container.h"

size_t SequenceContainer::g_nextSeqId = 0;

namespace
{
	class ParseException : public std::runtime_error 
	{
	public:
		ParseException(const std::string & what):
			std::runtime_error(what)
		{}
	};
}

const FastaRecord::Id FastaRecord::ID_NONE = 
			Id(std::numeric_limits<uint32_t>::max());


void SequenceContainer::readFasta(const std::string& fileName)
{
	std::vector<FastaRecord> records;
	this->getSequencesWithComplements(records, fileName);
	_seqIndex.reserve(records.size());

	//shuffling input reads
	std::vector<size_t> indicesPerm(records.size());
	for (size_t i = 0; i < indicesPerm.size(); ++i) indicesPerm[i] = i;
	std::random_shuffle(indicesPerm.begin(), indicesPerm.end());

	for (size_t i : indicesPerm)
	{
		_seqIndex[records[i].id] = std::move(records[i]);
	}
}


const FastaRecord& 
	SequenceContainer::addSequence(const DnaSequence& sequence, 
								   const std::string& description)
{
	FastaRecord::Id newId = FastaRecord::Id(g_nextSeqId);
	FastaRecord fwdRecord(sequence, "+" + description, newId);
	_seqIndex[fwdRecord.id] = std::move(fwdRecord);

	FastaRecord revRecord(sequence.complement(), "-" + description, 
						  newId.rc());
	_seqIndex[revRecord.id] = std::move(revRecord);

	g_nextSeqId += 2;
	return _seqIndex[newId];
}

size_t SequenceContainer::getSequences(std::vector<FastaRecord>& record, 
									   const std::string& fileName)
{
	std::string buffer;
	std::string header; 
	std::string sequence;
	std::ifstream inputStream(fileName);
	if (!inputStream.is_open())
	{
		throw std::runtime_error("Can't open reads file");
	}

	int line = 1;
	record.clear();

	try
	{
		while(!inputStream.eof())
		{
			std::getline(inputStream, buffer, '\n');
			if (buffer.empty()) continue;
			if (*buffer.rbegin() == '\r') buffer.erase(buffer.size() - 1);

			if (buffer[0] == '>')
			{
				if (!header.empty())
				{
					if (sequence.empty()) throw ParseException("empty sequence");

					sequence.shrink_to_fit();
					record.emplace_back(DnaSequence(sequence), header, 
										FastaRecord::Id(g_nextSeqId));
					g_nextSeqId += 2;
					sequence.clear();
					header.clear();
				}
				this->validateHeader(buffer);
				header = buffer;
			}
			else
			{
				this->validateSequence(buffer);
				std::copy(buffer.begin(), buffer.end(), 
						  std::back_inserter(sequence));
			}

			++line;
		}
		
		if (sequence.empty()) throw ParseException("empty sequence");
		record.emplace_back(DnaSequence(sequence), header, 
							FastaRecord::Id(g_nextSeqId));
		g_nextSeqId += 2;

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

	records.reserve(records.size() * 2);
	std::vector<FastaRecord> complements;
	for (auto &record : records)
	{
		std::string header = "-" + record.description;
		record.description = "+" + record.description;

		DnaSequence revComplement = record.sequence.complement();
		complements.emplace_back(revComplement, header, record.id.rc());
	}

	for (auto& rec : complements)
	{
		records.push_back(std::move(rec));
	}

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

void SequenceContainer::validateSequence(std::string& sequence)
{
	const std::string VALID_CHARS = "ACGT";
	for (size_t i = 0; i < sequence.length(); ++i)
	{
		if (VALID_CHARS.find(toupper(sequence[i])) == std::string::npos) 
		{
			sequence[i] = VALID_CHARS[rand() % 4];
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
			fout << rec.sequence.substr(c, FASTA_SLICE).str() << std::endl;
	}
}
