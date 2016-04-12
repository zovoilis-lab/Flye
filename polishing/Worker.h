#pragma once
#include <string>
#include <vector>
#include "ScoringMatrix.h"
#include "Alignment.h"
#include <math.h>

struct Record 
{
	std::string read;
	std::string methodUsed;
	double score;

	int  del_index;

	int  sub_index;
	char sub_letter;

	int  ins_index;
	char ins_letter;

	Record():
		score(0.0f), del_index(-1), sub_index(-1), sub_letter('*'),
		ins_index(-1), ins_letter('*')	
	{}
};

class Worker {
public:
	Worker(const std::string& scoreMatPath);
	void runOneToAll(const std::string& candidate, Record& rec);
	void run(const std::string& dataPath, const std::string& format);
	void run(size_t start, size_t stop, 
			 const std::string& dataPath, const std::string& format);

public:
	ScoringMatrix  _scoreMat;
	std::vector<std::string> _reads;
	int _filePos;
	//Alignment align;

	void readFasta(std::vector<std::string>& reads, 
				   const std::string& path);
	std::string readFastaSpecial(std::vector<std::string>& reads, 
							     const std::string& path);
	//std::vector<std::string> split(const std::string& str, char delim);
	//std::vector<std::string> split(const std::string& str, char delim, 
	//							   std::vector<std::string>& e);
	void progressUpdate(int start, int stop, int initial, 
						int interval, double& prev);
	void outputRecord(Record rec);
	void outputSeparator();
};

