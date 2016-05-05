//(c) 2016 by Authors
//This file is a part of ABruijn program.
//Released under the BSD license (see LICENSE file)

#pragma once

#include "bubble.h"
#include "subs_matrix.h"

class GeneralPolisher
{
public:
	GeneralPolisher(const SubstitutionMatrix& subsMatrix):
		_subsMatrix(subsMatrix)
	{}
	void polishBubble(Bubble& bubble);

private:
	StepInfo makeStep(const std::string& candidate, 
					  const std::vector<std::string>& branches);

	const SubstitutionMatrix& _subsMatrix;
};
