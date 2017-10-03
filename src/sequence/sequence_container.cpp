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

const FastaRecord::Id FastaRecord::ID_NONE = 
			Id(std::numeric_limits<uint32_t>::max());


bool SequenceContainer::isFasta(const std::string& fileName)
{
	size_t dotPos = fileName.rfind(".");
	if (dotPos == std::string::npos)
	{
		throw ParseException("Incorrect file suffix");
	}
	std::string suffix = fileName.substr(dotPos + 1);
	if (suffix == "fasta" || suffix == "fa")
	{
		return true;
	}
	else if (suffix == "fastq" || suffix == "fq")
	{
		return false;
	}
	throw ParseException("Incorrect file suffix");
}

void SequenceContainer::loadFromFile(const std::string& fileName)
{
	std::vector<FastaRecord> records;

	if (this->isFasta(fileName))
	{
		this->readFasta(records, fileName);
	}
	else
	{
		this->readFastq(records, fileName);
	}

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
	
	//shuffling input reads
	std::vector<size_t> indicesPerm(records.size());
	for (size_t i = 0; i < indicesPerm.size(); ++i) indicesPerm[i] = i;
	std::random_shuffle(indicesPerm.begin(), indicesPerm.end());
	//

	_seqIndex.reserve(records.size());
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

size_t SequenceContainer::readFasta(std::vector<FastaRecord>& record, 
									const std::string& fileName)
{
	std::string buffer;
	std::string header; 
	std::string sequence;
	std::ifstream inputStream(fileName);
	if (!inputStream.is_open())
	{
		throw ParseException("Can't open reads file");
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
		if (header.empty())
		{
			throw ParseException("Fasta fromat error");
		}
		record.emplace_back(DnaSequence(sequence), header, 
							FastaRecord::Id(g_nextSeqId));
		g_nextSeqId += 2;

	}
	catch (ParseException& e)
	{
		std::stringstream ss;
		ss << "parse error in " << fileName << " on line " << line << ": " << e.what();
		throw ParseException(ss.str());
	}

	return record.size();
}

size_t SequenceContainer::readFastq(std::vector<FastaRecord>& record, 
									const std::string& fileName)
{
	std::string buffer;
	std::string header; 
	std::ifstream inputStream(fileName);
	if (!inputStream.is_open())
	{
		throw ParseException("Can't open reads file");
	}

	int line = 1;
	int stateCounter = 0;
	record.clear();
	try
	{
		while(!inputStream.eof())
		{
			std::getline(inputStream, buffer, '\n');
			if (buffer.empty()) continue;
			if (*buffer.rbegin() == '\r') buffer.erase(buffer.size() - 1);

			if (stateCounter == 0)
			{
				if (buffer[0] != '@') throw ParseException("Fastq format error");
				header = buffer.substr(1);
				this->validateHeader(header);
			}
			else if (stateCounter == 1)
			{
				this->validateSequence(buffer);
				record.emplace_back(DnaSequence(buffer), header, 
									FastaRecord::Id(g_nextSeqId));
				g_nextSeqId += 2;
			}
			else if (stateCounter == 2)
			{
				if (buffer[0] != '+') throw ParseException("Fastq fromat error");
			}
			stateCounter = (stateCounter + 1) % 4;
			++line;
		}
	}
	catch (ParseException& e)
	{
		std::stringstream ss;
		ss << "parse error in " << fileName << " on line " << line << ": " << e.what();
		throw ParseException(ss.str());
	}

	return record.size();
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
