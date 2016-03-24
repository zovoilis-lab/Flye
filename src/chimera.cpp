#include "chimera.h"
#include <vector>

void ChimeraDetector::detectChimeras(const OverlapDetector& ovlpDetector,
									 const SequenceContainer& seqContainer)
{
	static const int WINDOW = 100;
	for (auto& seqHash : seqContainer.getIndex())
	{
		auto& overlaps = ovlpDetector.getOverlapIndex().at(seqHash.first);
		int numWindows = seqHash.second.sequence.length() / WINDOW;
		std::vector<int> localCoverage(numWindows, 0);

		int total = 0;
		for (auto& ovlp : overlaps)
		{
			for (int pos = (ovlp.curBegin + _maximumJump) / WINDOW; 
				 pos < (ovlp.curEnd - _maximumJump) / WINDOW; ++pos)
				++localCoverage[pos];
				++total;
		}

		float avgCoverage = float(total) / localCoverage.size();
		bool chimeric = false;
		for (int cov : localCoverage)
		{
			if (cov < 0.1 * avgCoverage)
			{
				chimeric = true;
				break;
			}
		}
		if (chimeric)
		{
			std::cout << "Chimeric: " << seqHash.second.description 
					  << std::endl;
		}
	}
}
