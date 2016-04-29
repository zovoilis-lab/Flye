//(c) 2016 by Authors
//This file is a part of ABruijn program.
//Released under the BSD license (see LICENSE file)

#include <vector>
#include <iostream>

#include "chimera.h"
#include "utility.h"

void ChimeraDetector::detectChimeras(const OverlapDetector& ovlpDetector,
									 const SequenceContainer& seqContainer)
{
	LOG_PRINT("Detecting chimeric sequences");

	static const int WINDOW = 100;
	static const float COV_THRESHOLD = 0.1f;
	static const float MAGIC_NUMBER = 2.5f;
	const int FLANK = (_maximumJump + _maximumOverhang) / WINDOW;

	std::unordered_map<FastaRecord::ReadIdType, 
					   std::vector<int>> localCoverage;
	for (auto& seqHash : seqContainer.getIndex())
	{
		int numWindows = seqHash.second.sequence.length() / WINDOW;
		if (numWindows - 2 * FLANK <= 0) continue;
		localCoverage[seqHash.first].assign(numWindows - 2 * FLANK, 0);

		for (auto& ovlp : ovlpDetector.getOverlapIndex().at(seqHash.first))
		{
			for (int pos = (ovlp.curBegin + _maximumJump) / WINDOW; 
				 pos < (ovlp.curEnd - _maximumJump) / WINDOW; ++pos)
			{
				if (pos - FLANK >= 0 && 
					pos - FLANK < (int)localCoverage[seqHash.first].size())
				{
					++localCoverage[seqHash.first][pos - FLANK];
				}
			}
		}
	}

	float covSum = 0;
	int numWindows = 0;
	for (auto& seqHash : seqContainer.getIndex())
	{
		for (int cov : localCoverage[seqHash.first])
		{
			covSum += cov;
			++numWindows;
		}
	}
	_coverage = (numWindows != 0) ? covSum / numWindows : 1;
	_coverage *= MAGIC_NUMBER;
	DEBUG_PRINT("Estimated coverage: " << _coverage);

	for (auto& seqHash : seqContainer.getIndex())
	{
		bool chimeric = false;
		for (int cov : localCoverage[seqHash.first])
		{
			if (cov < std::max(COV_THRESHOLD * _coverage, 1.0f))
			{
				chimeric = true;
				break;
			}
		}
		if (chimeric)
		{
			//DEBUG_PRINT("Chimeric: " << seqHash.second.description);
			_chimeras.insert(seqHash.first);
		}
	}
	LOG_PRINT(_chimeras.size() << " sequences were marked as chimeric");
}
