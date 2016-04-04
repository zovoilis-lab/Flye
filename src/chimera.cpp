#include "chimera.h"
#include <vector>

void ChimeraDetector::detectChimeras(const OverlapDetector& ovlpDetector,
									 const SequenceContainer& seqContainer)
{
	static const int WINDOW = 100;
	static const float COV_THRESHOLD = 0.1f;
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

		//float avgCoverage = float(total) / localCoverage.size();
		bool chimeric = false;
		size_t flank = (_maximumJump + _maximumOverhang) / WINDOW;
		for (size_t i = flank; i < localCoverage.size() - flank; ++i)
		{
			if (localCoverage[i] < int(COV_THRESHOLD * _coverage))
			//if (localCoverage[i] <= 2)
			{
				chimeric = true;
				break;
			}
		}
		if (chimeric)
		{
			std::cout << "Chimeric: " << seqHash.second.description 
					  << std::endl;
			_chimeras.insert(seqHash.first);
		}
	}
}
