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
	Logger::get().debug() << "Estimating overlap coverage";

	const int NUM_SAMPLES = 1000;
	const int sampleRate = _seqContainer.getIndex().size() / NUM_SAMPLES;
	static const int minCoverage = _inputCoverage / 
					(int)Config::get("max_coverage_drop_rate") + 1;
	static const int maxCoverage = _inputCoverage * 
					(int)Config::get("max_coverage_drop_rate");
	static const int flankSize = 0;

	std::unordered_map<int32_t, int32_t> readHist;

	for (auto& seq : _seqContainer.getIndex())
	{
		if (rand() % sampleRate) continue;
		auto coverage = this->getReadCoverage(seq.first);
		bool nonZero = false;
		for (auto c : coverage) nonZero |= (c != 0);
		if (!nonZero) continue;


		//int64_t sum = 0;
		//int64_t num = 0;
		for (size_t i = flankSize; i < coverage.size() - flankSize; ++i)
		{
			if (minCoverage < coverage[i] && coverage[i] < maxCoverage)
			{
				//fout << coverage[i] << std::endl;
				++readHist[coverage[i]];
				//sum += coverage[i];
				//++num;
			}
		}
		//fout << _ovlpContainer.lazySeqOverlaps(seq.first).size() << std::endl;

		/*std::string covStr;
		for (int cov : coverage)
		{
			covStr += std::to_string(cov) + " ";
		}
		Logger::get().debug() << "\t" << covStr;*/
		//++sampleNum;
		
		/*if (num)
		{
			++readHist[sum / num];
			++sampleNum;
		}*/

		//if (sampleNum >= NUM_SAMPLES) break;
	}

	if (readHist.empty())
	{
		Logger::get().warning() << "No overlaps found!";
		_overlapCoverage = 0;
	}

	int32_t maxCount = 0;
	int32_t peakCoverage = 0;
	for (auto& histIt : readHist)
	{
		if (histIt.second > maxCount)
		{
			maxCount = histIt.second;
			peakCoverage = histIt.first;
		}
	}
	_overlapCoverage = peakCoverage;
	Logger::get().info() << "Overlap-based coverage: " << _overlapCoverage;
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

	for (auto cov : coverage)
	{
		if (cov == 0) return true;

		if ((float)_overlapCoverage / cov > MAX_DROP_RATE) 
		{
			return true;
		}
	}

	//chimera detection based self-overlaps (typical PacBio pattern)
	if (_ovlpContainer.hasSelfOverlaps(readId))
	{
		//Logger::get().info() << "Self-ovlp!";
		return true;
	}

	return false;
}
