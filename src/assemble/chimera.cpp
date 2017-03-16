//(c) 2016 by Authors
//This file is a part of ABruijn program.
//Released under the BSD license (see LICENSE file)

#include <vector>
#include <iostream>
#include <unordered_set>
#include <unordered_map>

#include "../sequence/config.h"
#include "chimera.h"
#include "logger.h"


/*
void ChimeraDetector::detectChimeras()
{
	Logger::get().debug() << "Overlap-estimated coverage: "
						  << this->estimateOverlapCoverage();
	for (auto& seqHash : _seqContainer.getIndex())
	{
		if (this->testReadByCoverage(seqHash.first))
		{
			_chimeras.insert(seqHash.first);
			_chimeras.insert(seqHash.first.rc());
			//Logger::get().debug() << "Chimeric!";
		}
	}

	Logger::get().debug() << _chimeras.size() / 2 
						  << " sequences were marked as chimeric";
}*/

bool ChimeraDetector::isChimeric(FastaRecord::Id readId)
{
	if (!_chimeras.count(readId))
	{
		bool result = this->testReadByCoverage(readId);
		_chimeras[readId] = result;
		_chimeras[readId.rc()] = result;
	}
	return _chimeras[readId];
}

/*
int ChimeraDetector::estimateOverlapCoverage()
{
	size_t total = 0;
	for (auto& seqHash : _seqContainer.getIndex())
	{
		total += _ovlpDetector.getOverlapIndex().at(seqHash.first).size();
	}
	return total / _seqContainer.getIndex().size();
}*/

bool ChimeraDetector::testReadByCoverage(FastaRecord::Id readId)
{
	//Logger::get().debug() << _seqContainer.seqName(readId);
	//self overlaps
	for (auto& ovlp : _ovlpContainer.getSeqOverlaps(readId))
	{
		if (ovlp.curId == ovlp.extId.rc()) 
		{
			Logger::get().debug() << "yes!";
			return true;
		}
	}

	static const int WINDOW = Constants::chimeraWindow;
	const int FLANK = (Constants::maximumJump + 
					   Constants::maximumOverhang) / WINDOW;

	std::vector<int> coverage;
	int numWindows = _seqContainer.seqLen(readId) / WINDOW;
	if (numWindows - 2 * FLANK <= 0) return false;

	coverage.assign(numWindows - 2 * FLANK, 0);
	for (auto& ovlp : _ovlpContainer.getSeqOverlaps(readId))
	{
		if (ovlp.curId == ovlp.extId.rc()) continue;

		for (int pos = ovlp.curBegin / WINDOW + 1; 
			 pos < ovlp.curEnd / WINDOW; ++pos)
		{
			if (pos - FLANK >= 0 && 
				pos - FLANK < (int)coverage.size())
			{
				++coverage[pos - FLANK];
			}
		}
	}

	/*
	std::string covStr;
	for (int cov : coverage)
	{
		covStr += std::to_string(cov) + " ";
	}
	Logger::get().debug() << covStr;
	*/

	//exclude low-coverage reads
	int maxCoverage = *std::max_element(coverage.begin(), coverage.end());
	if (maxCoverage < Constants::maxCoverageDropRate) return false;

	for (size_t i = 0; i < coverage.size() - 1; ++i)
	{
		int low = std::min(coverage[i], coverage[i + 1]);
		int high = std::max(coverage[i], coverage[i + 1]);
		if (low == 0 && high > 0) return true;
		if (low != 0 && high / low > Constants::maxCoverageDropRate) return true;
	}
	return false;
}
