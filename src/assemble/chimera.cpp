//(c) 2016 by Authors
//This file is a part of ABruijn program.
//Released under the BSD license (see LICENSE file)

#include <vector>
#include <iostream>
#include <unordered_set>
#include <unordered_map>

#include "../common/config.h"
#include "../common/logger.h"
#include "chimera.h"



bool ChimeraDetector::isChimeric(FastaRecord::Id readId)
{
	if (!_chimeras.contains(readId))
	{
		bool result = this->testReadByCoverage(readId);
		_chimeras[readId] = result;
		_chimeras[readId.rc()] = result;
	}
	return _chimeras[readId];
}

void ChimeraDetector::estimateGlobalCoverage()
{
	const int NUM_SAMPLES = 1000;
	const int MAX_COVERAGE = _inputCoverage * 
					(int)Config::get("max_coverage_drop_rate");
	int sampleNum = 0;
	std::vector<int32_t> readMeans;

	for (auto& seq : _seqContainer.getIndex())
	{
		auto coverage = this->getReadCoverage(seq.first);

		bool isZero = false;
		int64_t sum = 0;
		for (auto c : coverage)
		{
			if (c < MAX_COVERAGE)
			{
				isZero |= (c == 0);
				sum += c;
			}
		}
		
		if (!isZero)
		{
			readMeans.push_back(sum / coverage.size());
			++sampleNum;
		
			/*std::string covStr;
			for (int cov : coverage)
			{
				covStr += std::to_string(cov) + " ";
			}
			Logger::get().debug() << "\t" << covStr;*/
		}


		if (sampleNum >= NUM_SAMPLES) break;
	}

	int64_t sum = 0;
	for (auto m : readMeans) sum += m;

	if (!readMeans.empty())
	{
		_overlapCoverage = sum / readMeans.size();
		Logger::get().debug() << "Mean read coverage: " << _overlapCoverage;
	}
}

std::vector<int32_t> ChimeraDetector::getReadCoverage(FastaRecord::Id readId)
{
	static const int WINDOW = Config::get("chimera_window");
	const int FLANK = (int)Config::get("maximum_overhang") / WINDOW;

	std::vector<int> coverage;
	int numWindows = _seqContainer.seqLen(readId) / WINDOW;
	if (numWindows - 2 * FLANK <= 0) return {0};

	coverage.assign(numWindows - 2 * FLANK, 0);
	for (auto& ovlp : _ovlpContainer.lazySeqOverlaps(readId))
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

	return coverage;
}


bool ChimeraDetector::testReadByCoverage(FastaRecord::Id readId)
{
	auto coverage = this->getReadCoverage(readId);
	const float MAX_DROP_RATE = Config::get("max_coverage_drop_rate");
	/*std::string covStr;
	for (int cov : coverage)
	{
		covStr += std::to_string(cov) + " ";
	}*/
	//Logger::get().debug() << "\t" << _seqContainer.seqName(readId) << covStr;

	//int LOW_COV_THRESHOLD = 10;
	//int maxCoverage = *std::max_element(coverage.begin(), coverage.end());
	for (auto cov : coverage)
	{
		if (cov == 0) return true;

		if ((float)_overlapCoverage / cov > MAX_DROP_RATE) 
		{
			//Logger::get().debug() << "Chimeric: " 
			//	<< _seqContainer.seqName(readId) << covStr;
			return true;
		}
	}

	//self overlaps
	for (auto& ovlp : _ovlpContainer.lazySeqOverlaps(readId))
	{
		if (ovlp.curId == ovlp.extId.rc()) 
		{
			//Logger::get().debug() << "Self-ovlp: " 
			//	<< _seqContainer.seqName(readId) << covStr;
			return true;
		}
	}

	return false;
}
