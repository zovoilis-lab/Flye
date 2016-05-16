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

	DEBUG_PRINT("Overlap-estimated coverage: " << 
				this->estimateOverlapCoverage());
	for (auto& seqHash : _seqContainer.getIndex())
	{
		//if (this->testReadByClusters(seqHash.first))
		//if (this->testReadByCoverage(seqHash.first))
		if (this->testSelfOverlap(seqHash.first))
		{
			//DEBUG_PRINT("Chimeric: " << seqHash.second.description);
			_chimeras.insert(seqHash.first);
		}
	}

	LOG_PRINT(_chimeras.size() / 2 << " sequences were marked as unreliable");
}

int ChimeraDetector::estimateOverlapCoverage()
{
	static const float MAGIC_25 = 2.5f;
	static const int WINDOW = 100;
	const int FLANK = (_maximumJump + _maximumOverhang) / WINDOW;

	std::unordered_map<FastaRecord::Id, 
					   std::vector<int>> localCoverage;
	for (auto& seqHash : _seqContainer.getIndex())
	{
		int numWindows = seqHash.second.sequence.length() / WINDOW;
		if (numWindows - 2 * FLANK <= 0) continue;
		localCoverage[seqHash.first].assign(numWindows - 2 * FLANK, 0);

		for (auto& ovlp : _ovlpDetector.getOverlapIndex().at(seqHash.first))
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
	for (auto& seqHash : _seqContainer.getIndex())
	{
		for (int cov : localCoverage[seqHash.first])
		{
			covSum += cov;
			++numWindows;
		}
	}
	int estCoverage = (numWindows != 0) ? covSum / numWindows : 1;
	return estCoverage * MAGIC_25;
}

bool ChimeraDetector::testReadByCoverage(FastaRecord::Id readId)
{
	static const int WINDOW = 100;
	static const float COV_THRESHOLD = 0.1f;
	const int FLANK = (_maximumJump + _maximumOverhang) / WINDOW;

	std::vector<int> localCoverage;
	int numWindows = _seqContainer.seqLen(readId) / WINDOW;
	if (numWindows - 2 * FLANK <= 0) return false;

	localCoverage.assign(numWindows - 2 * FLANK, 0);

	for (auto& ovlp : _ovlpDetector.getOverlapIndex().at(readId))
	{
		for (int pos = (ovlp.curBegin + _maximumJump) / WINDOW; 
			 pos < (ovlp.curEnd - _maximumJump) / WINDOW; ++pos)
		{
			if (pos - FLANK >= 0 && 
				pos - FLANK < (int)localCoverage.size())
			{
				++localCoverage[pos - FLANK];
			}
		}
	}

	for (int cov : localCoverage)
	{
		if (cov < std::max(COV_THRESHOLD * _coverage, 1.0f))
		{
			return true;
		}
	}
	return false;
}

bool ChimeraDetector::testSelfOverlap(FastaRecord::Id readId)
{
	auto& overlaps = _ovlpDetector.getOverlapIndex().at(readId);
	for (auto& ovlp : overlaps)
	{
		if (ovlp.extId == readId.rc()) return true;
	}
	return false;
}

bool ChimeraDetector::testReadByClusters(FastaRecord::Id readId)
{
	auto& overlaps = _ovlpDetector.getOverlapIndex().at(readId);
	std::unordered_map<FastaRecord::Id,
					   SetNode<FastaRecord::Id>*> extensions;
    for (auto& ovlp : overlaps) 
	{
		extensions[ovlp.extId] = new SetNode<FastaRecord::Id>(ovlp.extId);
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

	std::unordered_set<SetNode<FastaRecord::Id>*> clusters;
	for (auto& extHash : extensions) clusters.insert(findSet(extHash.second));

	//DEBUG_PRINT("Clusters: " << clusters.size());
	return clusters.size() != 1;
}
