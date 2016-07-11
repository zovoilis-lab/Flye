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
		if (this->testReadByCoverage(seqHash.first))
		{
			_chimeras.insert(seqHash.first);
		}
	}

	Logger::get().debug() << _chimeras.size() / 2 
						  << " sequences were marked as chimeric";
}

int ChimeraDetector::estimateOverlapCoverage()
{
	size_t total = 0;
	for (auto& seqHash : _seqContainer.getIndex())
	{
		total += _ovlpDetector.getOverlapIndex().at(seqHash.first).size();
	}
	return total / _seqContainer.getIndex().size();
}

bool ChimeraDetector::testReadByCoverage(FastaRecord::Id readId)
{
	static const int WINDOW = _maximumJump;
	const int FLANK = (_maximumJump + _maximumOverhang) / WINDOW;

	std::vector<int> coverage;
	int numWindows = _seqContainer.seqLen(readId) / WINDOW;
	if (numWindows - 2 * FLANK <= 0) return false;

	coverage.assign(numWindows - 2 * FLANK, 0);
	for (auto& ovlp : _ovlpDetector.getOverlapIndex().at(readId))
	{
		if (ovlp.curId == ovlp.extId.rc()) return true;

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

	Logger::get().debug() << _seqContainer.seqName(readId);
	std::string covStr;
	for (int cov : coverage)
	{
		covStr += std::to_string(cov) + " ";
	}
	Logger::get().debug() << covStr;

	bool zeroStrip = false;
	for (size_t i = 0; i < coverage.size() - 1; ++i)
	{
		if (!zeroStrip && coverage[i] != 0 && coverage[i + 1] == 0)
		{
			zeroStrip = true;
		}
		if (zeroStrip && coverage[i + 1] != 0)
		{
			Logger::get().debug() << "Chimeric!";
			return true;
		}
	}
	return false;
}
