//(c) 2016 by Authors
//This file is a part of ABruijn program.
//Released under the BSD license (see LICENSE file)

#pragma once

#include <string>
#include <vector>

enum StepType {StepNone, StepIns, StepSub, StepDel, StepHopo};

struct StepInfo 
{
	std::string sequence;
	StepType methodUsed;
	double 	 score;
	int 	 changedIndex;
	char	 changedLetter;

	StepInfo():
		methodUsed(StepNone), score(0.0f), 
		changedIndex(0), changedLetter('*')
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
