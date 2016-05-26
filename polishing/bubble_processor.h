//(c) 2016 by Authors
//This file is a part of ABruijn program.
//Released under the BSD license (see LICENSE file)

#pragma once

#include <string>
#include <vector>
#include <cmath>
#include <mutex>

#include "subs_matrix.h"
#include "bubble.h"
#include "general_polisher.h"
#include "homo_polisher.h"
#include "utility.h"


class BubbleProcessor 
{
public:
	BubbleProcessor(const std::string& subsMatPath,
					const std::string& hopoMatrixPath);

	void polishAll(const std::string& dataPath, int numThreads);
	void writeConsensuses(const std::string& fileName);
	void writeLog(const std::string& fileName);

private:
	void parallelWorker();
	void readBubbles(const std::string& fileName);

	const SubstitutionMatrix  _subsMatrix;
	const HopoMatrix 		  _hopoMatrix;
	const GeneralPolisher 	  _generalPolisher;
	const HomoPolisher 		  _homoPolisher;

	std::vector<Bubble>		  _bubbles;
	ProgressPercent 		  _progress;
	size_t					  _nextJob;
	std::mutex				  _stateMutex;
};
