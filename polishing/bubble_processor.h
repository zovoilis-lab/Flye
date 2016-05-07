//(c) 2016 by Authors
//This file is a part of ABruijn program.
//Released under the BSD license (see LICENSE file)

#pragma once
#include <string>

#include <vector>
#include <math.h>

#include "subs_matrix.h"
#include "alignment.h"
#include "bubble.h"


class BubbleProcessor 
{
public:
	BubbleProcessor(const std::string& subsMatPath,
					const std::string& hopoMatrixPath);

	void polishAll(const std::string& dataPath);
	void writeConsensuses(const std::string& fileName);
	void writeLog(const std::string& fileName);

private:
	SubstitutionMatrix  _subsMatrix;
	HopoMatrix _hopoMatrix;
	std::vector<Bubble> _bubbles;

	void readBubbles(const std::string& fileName);
};
