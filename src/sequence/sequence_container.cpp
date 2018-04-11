//(c) 2016 by Authors
//This file is a part of ABruijn program.
//Released under the BSD license (see LICENSE file)

#include <stdexcept>
#include <fstream>
#include <sstream>
#include <iostream>
#include <random>
#include <algorithm>
#include <zlib.h>

#include "sequence_container.h"
#include "../common/logger.h"

size_t SequenceContainer::g_nextSeqId = 0;

const FastaRecord::Id FastaRecord::ID_NONE = 
			Id(std::numeric_limits<uint32_t>::max());


bool SequenceContainer::isFasta(const std::string& fileName)
{
	std::string withoutGz = fileName;
	if (fileName.substr(fileName.size() - 3) == ".gz")
	{
		withoutGz = fileName.substr(0, fileName.size() - 3);
	}

	size_t dotPos = withoutGz.rfind(".");
	if (dotPos == std::string::npos)
	{
		throw ParseException("Can't identify input file type");
	}
	std::string suffix = withoutGz.substr(dotPos + 1);

	if (suffix == "fasta" || suffix == "fa")
	{
		return true;
	}
	else if (suffix == "fastq" || suffix == "fq")
	{
		return false;
	}
	throw ParseException("Can't identify input file type");
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

	_seqIndex.reserve(_seqIndex.size() + records.size());
	for (size_t i : indicesPerm)
	{
		_seqIndex[records[i].id] = std::move(records[i]);
	}
}


int SequenceContainer::computeNxStat(float fraction) const
{
	std::vector<int32_t> readLengths;
	int64_t totalLengh = 0;
	for (auto& read : _seqIndex) 
	{
		readLengths.push_back(read.second.sequence.length());
		totalLengh += read.second.sequence.length();
	}
	std::sort(readLengths.begin(), readLengths.end(),
			  [](int32_t a, int32_t b) {return a > b;});

	int32_t nx = 0;
	int64_t cummulativeLen = 0;
	for (auto l : readLengths)
	{
		cummulativeLen += l;
		if (cummulativeLen > fraction * totalLengh)
		{
			nx = l;
			break;
		}
	}
	return nx;
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
	size_t BUF_SIZE = 32 * 1024 * 1024;
	char* rawBuffer = new char[BUF_SIZE];
	auto* fd = gzopen(fileName.c_str(), "rb");
	if (!fd)
	{
		throw ParseException("Can't open reads file");
	}

	record.clear();
	int lineNo = 1;
	std::string header; 
	std::string sequence;
	try
	{
		while(!gzeof(fd))
		{
			//get a new line
			std::string nextLine;
			for (;;)
			{
				char* read = gzgets(fd, rawBuffer, BUF_SIZE);
				if (!read) break;
				nextLine += read;
				if (nextLine.empty()) break;
				if (nextLine.back() == '\n')
				{
					nextLine.pop_back();
					break;
				}
			}

			if (nextLine.empty()) continue;
			if (nextLine.back() == '\r') nextLine.pop_back();

			if (nextLine[0] == '>')
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
				this->validateHeader(nextLine);
				header = nextLine;
			}
			else
			{
				this->validateSequence(nextLine);
				std::copy(nextLine.begin(), nextLine.end(), 
						  std::back_inserter(sequence));
			}

			++lineNo;
		}
		
		if (sequence.empty()) throw ParseException("empty sequence");
		if (header.empty())
		{
			throw ParseException("Fasta format error");
		}
		record.emplace_back(DnaSequence(sequence), header, 
							FastaRecord::Id(g_nextSeqId));
		g_nextSeqId += 2;

	}
	catch (ParseException& e)
	{
		std::stringstream ss;
		ss << "parse error in " << fileName << " on line " << lineNo << ": " << e.what();
		gzclose(fd);
		throw ParseException(ss.str());
	}

	delete[] rawBuffer;
	gzclose(fd);
	return record.size();
}

size_t SequenceContainer::readFastq(std::vector<FastaRecord>& record, 
									const std::string& fileName)
{

	size_t BUF_SIZE = 32 * 1024 * 1024;
	char* rawBuffer = new char[BUF_SIZE];
	auto* fd = gzopen(fileName.c_str(), "rb");
	if (!fd)
	{
		throw ParseException("Can't open reads file");
	}

	record.clear();
	int lineNo = 1;
	int stateCounter = 0;
	std::string header; 
	try
	{
		while (!gzeof(fd))
		{
			//get a new line
			std::string nextLine;
			for (;;)
			{
				char* read = gzgets(fd, rawBuffer, BUF_SIZE);
				if (!read) break;
				nextLine += read;
				if (nextLine.empty()) break;
				if (nextLine.back() == '\n')
				{
					nextLine.pop_back();
					break;
				}
			}

			if (nextLine.empty()) 
			{
				stateCounter = (stateCounter + 1) % 4;
				continue;
			}
			if (nextLine.back() == '\r') nextLine.pop_back();

			if (stateCounter == 0)
			{
				if (nextLine[0] != '@') throw ParseException("Fastq format error");
				header = nextLine;
				this->validateHeader(header);
			}
			else if (stateCounter == 1)
			{
				this->validateSequence(nextLine);
				record.emplace_back(DnaSequence(nextLine), header, 
									FastaRecord::Id(g_nextSeqId));
				g_nextSeqId += 2;
			}
			else if (stateCounter == 2)
			{
				if (nextLine[0] != '+') throw ParseException("Fastq fromat error");
			}
			stateCounter = (stateCounter + 1) % 4;
			++lineNo;
		}
	}
	catch (ParseException& e)
	{
		std::stringstream ss;
		ss << "parse error in " << fileName << " on line " << lineNo << ": " << e.what();
		gzclose(fd);
		throw ParseException(ss.str());
	}

	gzclose(fd);
	delete[] rawBuffer;
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
								   const std::string& filename)
{
	static const size_t FASTA_SLICE = 80;

	Logger::get().debug() << "Writing FASTA";
	FILE* fout = fopen(filename.c_str(), "w");
	if (!fout) throw std::runtime_error("Can't open " + filename);
	
	for (auto& rec : records)
	{
		std::string contigSeq;
		for (size_t c = 0; c < rec.sequence.length(); c += FASTA_SLICE)
		{
			contigSeq += rec.sequence.substr(c, FASTA_SLICE).str() + "\n";
		}
		std::string header = ">" + rec.description + "\n";
		fwrite(header.data(), sizeof(header.data()[0]), 
			   header.size(), fout);
		fwrite(contigSeq.data(), sizeof(contigSeq.data()[0]), 
			   contigSeq.size(), fout);
	}
}
