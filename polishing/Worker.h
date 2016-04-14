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
	void run(const std::string& dataPath);
	void writeConsensuses(const std::string& fileName);
	void writeLog(const std::string& fileName);

	struct StepInfo 
	{
		std::string read;
		std::string methodUsed;
		double score;

		int  del_index;
		int  sub_index;
		char sub_letter;
		int  ins_index;
		char ins_letter;

		StepInfo():
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
		std::vector<StepInfo> polishSteps;
	};
private:
	ScoringMatrix  _scoreMat;
	std::vector<Bubble> _bubbles;

	void processCandidate(const std::string& candidate, 
					  	 const std::vector<std::string>& branches,
					  	 StepInfo& rec);
	void readBubbles(const std::string& fileName);
};

