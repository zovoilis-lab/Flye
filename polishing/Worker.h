#pragma once
#include <string>
#include <vector>
#include "ScoringMatrix.h"
#include "Alignment.h"
#include <math.h>



class Worker 
{
public:
	Worker(const std::string& scoreMatPath);
	void run(const std::string& dataPath, const std::string& format);

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

	struct Bubble
	{
		std::string header;
		int position;

		std::string candidate;
		std::vector<std::string> branches;
	};
private:
	ScoringMatrix  _scoreMat;
	std::vector<Bubble> _bubbles;

	void runOneToAll(const std::string& candidate, 
					 const std::vector<std::string>& branches,
					 Record& rec);
	void readBubbles(const std::string& fileName);
	void outputRecord(const Record& rec);
	void outputSeparator();
};

