//(c) 2016 by Authors
//This file is a part of ABruijn program.
//Released under the BSD license (see LICENSE file)

#include <vector>
#include <iostream>
#include <unordered_set>
#include <unordered_map>
#include <iomanip>

#include "../common/config.h"
#include "../common/logger.h"
#include "chimera.h"



bool ChimeraDetector::isChimeric(FastaRecord::Id readId)
{
	if (!_chimeras.contains(readId))
	{
		bool result = this->testReadByCoverage(readId);
		_chimeras.insert(readId, result);
		_chimeras.insert(readId.rc(), result);
	}
	return _chimeras.find(readId);
}

void ChimeraDetector::estimateGlobalCoverage()
{
	Logger::get().debug() << "Estimating overlap coverage";

	int numSamples = std::min(1000, (int)_seqContainer.iterSeqs().size());
	int sampleRate = (int)_seqContainer.iterSeqs().size() / numSamples;
	//int minCoverage = _inputCoverage / 
	//				(int)Config::get("max_coverage_drop_rate") + 1;
	//int maxCoverage = _inputCoverage * 
	//				(int)Config::get("max_coverage_drop_rate");
	int flankSize = 0;

	std::unordered_map<int32_t, int32_t> readHist;
	std::vector<int32_t> covList;
	std::vector<float> ovlpDivergence;
	
	//std::ofstream fout("../cov_hist.txt");

	int64_t sum = 0;
	int64_t num = 0;
	for (auto& seq : _seqContainer.iterSeqs())
	{
		if (rand() % sampleRate) continue;
		auto coverage = this->getReadCoverage(seq.id);
		bool nonZero = false;
		for (auto c : coverage) nonZero |= (c != 0);
		if (!nonZero) continue;

		for (size_t i = flankSize; i < coverage.size() - flankSize; ++i)
		{
			{
				++readHist[coverage[i]];
				sum += coverage[i];
				++num;
				covList.push_back(coverage[i]);
			}
		}

		//getting divergence
		for (auto& ovlp : _ovlpContainer.lazySeqOverlaps(seq.id)) 
		{
			ovlpDivergence.push_back(ovlp.seqDivergence);
		}
	}

	if (readHist.empty())
	{
		Logger::get().warning() << "No overlaps found!";
		_overlapCoverage = 0;
	}
	else
	{
		_overlapCoverage = median(covList);
	}

	Logger::get().info() << "Overlap-based coverage: " << _overlapCoverage;
	Logger::get().info() << "Median read-read divergence: "
		<< std::setprecision(2)
		<< median(ovlpDivergence) << " (Q10 = " << quantile(ovlpDivergence, 10)
		<< ", Q90 = " << quantile(ovlpDivergence, 90) << ")";
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

	//int64_t sumCov = 0;
	int maxCov = 0;
	for (auto cov : coverage)
	{
		//sumCov += cov;
		maxCov = std::max(maxCov, cov);
	}
	//int meanCov = sumCov ? sumCov / coverage.size() : 0;
	int threshold = round((float)std::min(_overlapCoverage, maxCov) /
						  MAX_DROP_RATE);

	for (auto cov : coverage)
	{
		if (cov == 0) return true;

		if (cov < threshold)
		{
			return true;
		}
	}

	//chimera detection based self-overlaps (typical PacBio pattern)
	/*if (_ovlpContainer.hasSelfOverlaps(readId))
	{
		//Logger::get().info() << "Self-ovlp!";
		return true;
	}*/

	return false;
}
