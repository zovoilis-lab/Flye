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
		bool result = this->testReadByCoverage(readId) ||
					  _ovlpContainer.hasSelfOverlaps(readId); //typical PacBio pattern
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
			++readHist[coverage[i]];
			sum += coverage[i];
			++num;
			covList.push_back(coverage[i]);
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
}

std::vector<int32_t> ChimeraDetector::getReadCoverage(FastaRecord::Id readId) const
{
	static const int WINDOW = (int)Config::get("chimera_window");
	//const int FLANK = (int)Config::get("maximum_overhang") / WINDOW;
	const int FLANK = 0;

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

int ChimeraDetector::getLeftTrim(const std::vector<int32_t>& coverage) const
{
	static const int MAX_DROP_RATE = Config::get("max_coverage_drop_rate");

	int64_t sumCov = 0;
	for (auto c: coverage) sumCov += c;
	int32_t meanCov = !coverage.empty() ? sumCov / coverage.size() : 0;
	int32_t threshold = std::max(2, meanCov / MAX_DROP_RATE);
	
	size_t pos = 0;
	for (pos = 0; pos < coverage.size(); ++pos)
	{
		if (coverage[pos] > threshold) break;
	}

	return (pos + 1) * (int)Config::get("chimera_window");
}

int ChimeraDetector::getRightTrim(const std::vector<int32_t>& coverage) const
{
	static const int MAX_DROP_RATE = Config::get("max_coverage_drop_rate");

	int64_t sumCov = 0;
	for (auto c: coverage) sumCov += c;
	int32_t meanCov = !coverage.empty() ? sumCov / coverage.size() : 0;
	int32_t threshold = std::max(2, meanCov / MAX_DROP_RATE);
	
	size_t pos = 0;
	for (pos = 0; pos < coverage.size(); ++pos)
	{
		if (coverage[coverage.size() - 1 - pos] > threshold) break;
	}

	return (pos + 1) * (int)Config::get("chimera_window");
}

int ChimeraDetector::getRightTrim(FastaRecord::Id readId) const
{
	auto coverage = this->getReadCoverage(readId);
	return this->getRightTrim(coverage);
}

int ChimeraDetector::getLeftTrim(FastaRecord::Id readId) const
{
	auto coverage = this->getReadCoverage(readId);
	return this->getLeftTrim(coverage);
}

bool ChimeraDetector::testReadByCoverage(FastaRecord::Id readId)
{
	auto coverage = this->getReadCoverage(readId);
	static const int MAX_DROP_RATE = Config::get("max_coverage_drop_rate");

	/*std::string covStr;
	for (int cov : coverage)
	{
		covStr += std::to_string(cov) + " ";
	}
	Logger::get().debug() << "\t" << _seqContainer.seqName(readId) << " " << covStr;
	Logger::get().debug() << "Left trim: " << this->getLeftTrim(coverage);
	Logger::get().debug() << "Right trim: " << this->getRightTrim(coverage);*/
	
	//const int HARD_MIN_COV = 2;
	int32_t threshold = std::max(2, _overlapCoverage / MAX_DROP_RATE);

	size_t leftFlank = this->getLeftTrim(coverage) / 
						(int)Config::get("chimera_window");
	size_t rightFlank = this->getRightTrim(coverage) / 
						(int)Config::get("chimera_window");
	for (size_t i = leftFlank; i < coverage.size() - rightFlank; ++i)
	{
		if (coverage[i] < threshold) return true;
	}

	return false;
}
