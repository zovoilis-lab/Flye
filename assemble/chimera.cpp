//(c) 2016 by Authors
//This file is a part of ABruijn program.
//Released under the BSD license (see LICENSE file)

#include <vector>
#include <iostream>
#include <unordered_set>
#include <unordered_map>

#include "chimera.h"
#include "logger.h"


void ChimeraDetector::detectChimeras()
{
	Logger::get().debug() << "Overlap-estimated coverage: "
						  << this->estimateOverlapCoverage();
	for (auto& seqHash : _seqContainer.getIndex())
	{
		//if (this->testReadByCoverage(seqHash.first))
		if (this->testSelfOverlap(seqHash.first))
		{
			_chimeras.insert(seqHash.first);
		}

		//Logger::get().debug() << "Chimeric: " << seqHash.second.description;
	}

	Logger::get().debug() << _chimeras.size() / 2 
						  << " sequences were marked as chimeric";
}

int ChimeraDetector::estimateOverlapCoverage()
{
	static const int WINDOW = 100;
	const int FLANK = _minimumOverlap / WINDOW;

	std::unordered_map<FastaRecord::Id, 
					   std::vector<int>> localCoverage;
	for (auto& seqHash : _seqContainer.getIndex())
	{
		int numWindows = seqHash.second.sequence.length() / WINDOW;
		if (numWindows - 2 * FLANK <= 0) continue;
		localCoverage[seqHash.first].assign(numWindows - 2 * FLANK, 0);

		for (auto& ovlp : _ovlpDetector.getOverlapIndex().at(seqHash.first))
		{
			for (int pos = ovlp.curBegin / WINDOW; 
				 pos < ovlp.curEnd / WINDOW; ++pos)
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
	return estCoverage;
}

bool ChimeraDetector::testReadByCoverage(FastaRecord::Id readId)
{
	static const int WINDOW = 100;
	static const float COV_THRESHOLD = 0.1f;
	const int FLANK = (_maximumJump + _maximumOverhang) / WINDOW;

	std::vector<int> coverage;
	int numWindows = _seqContainer.seqLen(readId) / WINDOW;
	if (numWindows - 2 * FLANK <= 0) return false;

	coverage.assign(numWindows - 2 * FLANK, 0);
	for (auto& ovlp : _ovlpDetector.getOverlapIndex().at(readId))
	{
		if (ovlp.curId == ovlp.extId.rc()) continue;

		for (int pos = (ovlp.curBegin + _maximumJump) / WINDOW; 
			 pos < (ovlp.curEnd - _maximumJump) / WINDOW; ++pos)
		{
			if (pos - FLANK >= 0 && 
				pos - FLANK < (int)coverage.size())
			{
				++coverage[pos - FLANK];
			}
		}
	}

	int sum = 0;
	for (int cov : coverage) sum += cov;
	int minCoverage = (float)sum / coverage.size() * COV_THRESHOLD;

	for (int cov : coverage)
	{
		if (cov < minCoverage)
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
