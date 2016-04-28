//(c) 2016 by Authors
//This file is a part of ABruijn program.
//Released under the BSD license (see LICENSE file)

#include <vector>
#include <iostream>
#include <unordered_set>
#include <unordered_map>

#include "chimera.h"
#include "utility.h"
#include "disjoint_set.h"


void ChimeraDetector::detectChimeras()
{
	LOG_PRINT("Detecting chimeric sequences");

	static const int WINDOW = 100;
	static const float COV_THRESHOLD = 0.1f;
	static const float MAGIC_NUMBER = 2.5f;
	const size_t FLANK = (_maximumJump + _maximumOverhang) / WINDOW;

	std::unordered_map<FastaRecord::ReadIdType, 
					   std::vector<int>> localCoverage;
	for (auto& seqHash : _seqContainer.getIndex())
	{
		auto& overlaps = _ovlpDetector.getOverlapIndex().at(seqHash.first);
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
	for (auto& seqHash : _seqContainer.getIndex())
	{
		for (size_t i = FLANK; 
			 i < localCoverage[seqHash.first].size() - FLANK; ++i)
		{
			covSum += localCoverage[seqHash.first][i];
			++numWindows;
		}
	}
	_coverage = (numWindows != 0) ? covSum / numWindows : 1;
	_coverage *= MAGIC_NUMBER;
	DEBUG_PRINT("Estimated coverage: " << _coverage);

	for (auto& seqHash : _seqContainer.getIndex())
	{
		bool chimeric = false;
		for (size_t i = FLANK; 
			 i < localCoverage[seqHash.first].size() - FLANK; ++i)
		{
			if (localCoverage[seqHash.first][i] < 
				std::max(COV_THRESHOLD * _coverage, 2.0f))
			{
				chimeric = true;
				break;
			}
		}
		if (chimeric)
		//if (this->testRead(seqHash.first))
		{
			//DEBUG_PRINT("Chimeric: " << seqHash.second.description);
			_chimeras.insert(seqHash.first);
		}
	}
	LOG_PRINT(_chimeras.size() << " sequences were marked as chimeric");
}

bool ChimeraDetector::testRead(FastaRecord::ReadIdType readId)
{
	auto& overlaps = _ovlpDetector.getOverlapIndex().at(readId);
	std::unordered_map<FastaRecord::ReadIdType,
					   SetNode<FastaRecord::ReadIdType>*> extensions;
    for (auto& ovlp : overlaps) 
	{
		extensions[ovlp.extId] = new SetNode<FastaRecord::ReadIdType>(ovlp.extId);
	}

	for (auto& extHash : extensions)
	{
		auto& extOverlaps = _ovlpDetector.getOverlapIndex().at(extHash.first);
		for (auto& extOvlp : extOverlaps)
		{
			if (extensions.count(extOvlp.extId))
			{
				unionSet(extHash.second, extensions[extOvlp.extId]);
			}
		}
	}

	std::unordered_set<SetNode<FastaRecord::ReadIdType>*> clusters;
	for (auto& extHash : extensions) clusters.insert(findSet(extHash.second));

	//DEBUG_PRINT("Clusters: " << clusters.size());
	return clusters.size() != 1;
}
