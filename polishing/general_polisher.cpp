//(c) 2016 by Authors
//This file is a part of ABruijn program.
//Released under the BSD license (see LICENSE file)

#include "general_polisher.h"
#include "alignment.h"

void GeneralPolisher::polishBubbles(std::vector<Bubble> bubbles)
{
	int prevPercent = -1;
	int counterDone = 0;
	for (auto& bubble : bubbles)
	{
		++counterDone;
		int percent = 10 * counterDone / bubbles.size();
		if (percent > prevPercent)
		{
			std::cerr << percent * 10 << "% ";
			prevPercent = percent;
		}

		std::string prevCandidate;
		std::string curCandidate = bubble.candidate;

		while (curCandidate != prevCandidate)
		{
			prevCandidate = curCandidate;
			StepInfo rec = this->makeStep(curCandidate, bubble.branches);
			curCandidate = rec.sequence;
			bubble.polishSteps.push_back(rec);
		}
		bubble.candidate = curCandidate;
	}
	std::cerr << std::endl;
}

StepInfo GeneralPolisher::makeStep(const std::string& candidate, 
				   				   const std::vector<std::string>& branches) 
{
	Alignment align(branches.size(), _subsMatrix);
	StepInfo stepResult;
	
	//Global
	double score = 0;
	for (size_t i = 0; i < branches.size(); ++i) 
	{
		score += align.globalAlignment(candidate, branches[i], i);
	}

	stepResult.score = score;
	stepResult.sequence = candidate;

	//Deletion
	for (size_t del_index = 0; del_index < candidate.size(); del_index++) 
	{
		double score = 0;
		for (size_t i = 0; i < branches.size(); i++) 
		{
			score += align.addDeletion(i, del_index + 1);
		}

		if (score > stepResult.score) 
		{
			std::string str = candidate;
			stepResult.methodUsed = StepDel;
			stepResult.score = score;
			stepResult.sequence = str.erase(del_index, 1);
			stepResult.changedIndex = del_index;
		}
	}

	//Substitution
	char alphabet[4] = {'A', 'C', 'G', 'T'};
	for (size_t sub_index = 0; sub_index < candidate.size(); sub_index++) 
	{
		for (char letter : alphabet)
		{
			if (letter == candidate[sub_index]) continue;

			double score = 0;
			for (size_t i = 0; i < branches.size(); i++) 
			{
				score += align.addSubstitution(i, sub_index + 1, letter, 
											   branches[i]);
			}

			if (score > stepResult.score) 
			{
				std::string str = candidate;
				stepResult.methodUsed = StepSub;
				stepResult.score = score;
				str.erase(sub_index, 1);
				stepResult.sequence = str.insert(sub_index, 1, letter);
				stepResult.changedIndex = sub_index;
				stepResult.changedLetter = letter;
			}
		}
	}

	//Insertion
	for (size_t ins_index = 0; ins_index < candidate.size()+1; ins_index++) 
	{
		for (char letter : alphabet)
		{
			double score = 0;
			for (size_t i = 0; i < branches.size(); i++) 
			{
				score += align.addInsertion(i, ins_index + 1, letter, 
											branches[i]);		
			}

			if (score > stepResult.score) 
			{
				std::string str = candidate;
				stepResult.methodUsed = StepIns;
				stepResult.score = score;
				stepResult.sequence = str.insert(ins_index, 1, letter);
				stepResult.changedIndex = ins_index;
				stepResult.changedLetter = letter;
			}
		}
	}	

	return stepResult;
}


