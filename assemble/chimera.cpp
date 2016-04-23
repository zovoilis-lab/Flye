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
	const size_t FLANK = (_maximumJump + _maximumOverhang) / WINDOW;

	std::unordered_map<FastaRecord::ReadIdType, 
					   std::vector<int>> localCoverage;
	for (auto& seqHash : seqContainer.getIndex())
	{
		auto& overlaps = ovlpDetector.getOverlapIndex().at(seqHash.first);
		int numWindows = seqHash.second.sequence.length() / WINDOW;
		localCoverage[seqHash.first].assign(numWindows, 0);

		for (auto& ovlp : overlaps)
		{
			for (int pos = (ovlp.curBegin + _maximumJump) / WINDOW; 
				 pos < (ovlp.curEnd - _maximumJump) / WINDOW; ++pos)
			{
				++localCoverage[seqHash.first][pos];
			}
		}
	}

	float covSum = 0;
	int numWindows = 0;
	for (auto& seqHash : seqContainer.getIndex())
	{
		for (size_t i = FLANK; 
			 i < localCoverage[seqHash.first].size() - FLANK; ++i)
		{
			covSum += localCoverage[seqHash.first][i];
			++numWindows;
		}
	}
	_coverage = covSum / numWindows;
	//float estCov = covSum / numWindows;
	//LOG_PRINT("Estimated coverage: " << estCov);
	//_coverage = estCov;

	for (auto& seqHash : seqContainer.getIndex())
	{
		bool chimeric = false;
		for (size_t i = FLANK; 
			 i < localCoverage[seqHash.first].size() - FLANK; ++i)
		{
			if (localCoverage[seqHash.first][i] < 
				int(COV_THRESHOLD * _coverage))
			{
				chimeric = true;
				break;
			}
		}
		if (chimeric)
		{
			DEBUG_PRINT("Chimeric: " << seqHash.second.description);
			_chimeras.insert(seqHash.first);
		}
	}
}
