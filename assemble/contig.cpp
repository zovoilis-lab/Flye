//(c) 2016 by Authors
//This file is a part of ABruijn program.
//Released under the BSD license (see LICENSE file)

#include <cassert>
#include <algorithm>
#include <fstream>

#include "contig.h"
#include "logger.h"

void ContigGenerator::generateContigs()
{
	Logger::get().info() << "Generating contig sequences";
	for (const ContigPath& path : _extender.getContigPaths())
	{
		if (path.reads.size() < 2) continue;

		std::vector<FastaRecord> contigParts;
		if (path.circular)
		{
			_contigs.push_back(this->generateCircular(path));
		}
		else
		{
			_contigs.push_back(this->generateLinear(path));
		}
	}
}


std::vector<FastaRecord> 
ContigGenerator::generateLinear(const ContigPath& path)
{
	std::vector<FastaRecord> contigParts;
	auto prevSwitch = std::make_pair(0, 0);
	int partNum = 1;
	for (size_t i = 0; i < path.reads.size(); ++i)
	{
		int32_t leftCut = prevSwitch.second;
		int32_t rightCut = _seqContainer.seqLen(path.reads[i]);
		if (i != path.reads.size() - 1)
		{
			auto curSwitch = this->getSwitchPositions(path.reads[i], 
													  path.reads[i + 1], 
													  prevSwitch.second);
			rightCut = curSwitch.first;
			prevSwitch = curSwitch;
		}

		std::string partSeq = _seqContainer.getIndex().at(path.reads[i])
							 .sequence.substr(leftCut, rightCut - leftCut);
		std::string partName = 
			(path.circular ? "circular_" : "linear_") + 
			std::to_string(_contigs.size()) +
			"_part_" + std::to_string(partNum) + "_" +
			_seqContainer.getIndex().at(path.reads[i]).description + 
			"[" + std::to_string(leftCut) + ":" + 
			std::to_string(rightCut) + "]";
		contigParts.push_back(FastaRecord(partSeq, partName, 
										  FastaRecord::Id(partNum)));
		++partNum;
	}
	return contigParts;
}

std::vector<FastaRecord> 
ContigGenerator::generateCircular(const ContigPath& path)
{
	std::vector<FastaRecord> contigParts;
	auto circPath = path.reads;
	circPath.push_back(path.reads[0]);
	circPath.push_back(path.reads[1]);

	int32_t initPivot = _seqContainer.seqLen(circPath[0]) / 2;
	auto firstSwitch = this->getSwitchPositions(circPath[0], circPath[1], 
												initPivot);
	auto prevSwitch = firstSwitch;
	int partNum = 1;
	for (size_t i = 1; i < circPath.size() - 1; ++i)
	{
		auto curSwitch = this->getSwitchPositions(circPath[i], 
												  circPath[i + 1], 
												  prevSwitch.second);
		int32_t leftCut = prevSwitch.second;
		int32_t rightCut = curSwitch.first;
		prevSwitch = curSwitch;

		if (i == circPath.size() - 2)	//finishing circle
		{
			rightCut = firstSwitch.first;
			if (rightCut - leftCut <= 0) rightCut = leftCut;
		}

		std::string partSeq = _seqContainer.getIndex().at(circPath[i])
							 .sequence.substr(leftCut, rightCut - leftCut);
		std::string partName = 
			(path.circular ? "circular_" : "linear_") + 
			std::to_string(_contigs.size()) +
			"_part_" + std::to_string(partNum) + "_" +
			_seqContainer.getIndex().at(circPath[i]).description + 
			"[" + std::to_string(leftCut) + ":" + 
			std::to_string(rightCut) + "]";
		contigParts.push_back(FastaRecord(partSeq, partName, 
										  FastaRecord::Id(partNum)));
		++partNum;
	}
	return contigParts;
}


void ContigGenerator::outputContigs(const std::string& fileName)
{
	std::vector<FastaRecord> allSeqs;
	for (auto& ctg : _contigs)
	{
		allSeqs.insert(allSeqs.end(), ctg.begin(), ctg.end());
	}
	SequenceContainer::writeFasta(allSeqs, fileName);
}


std::pair<int32_t, int32_t> 
ContigGenerator::getSwitchPositions(FastaRecord::Id leftRead, 
									FastaRecord::Id rightRead,
									int32_t prevSwitch)
{
	//find repeats
	std::unordered_map<Kmer, int> leftCount;
	for (auto& leftKmer : _vertexIndex.getIndexByRead().at(leftRead))
	{
		leftCount[leftKmer.kmer] += 1;
	}
	std::unordered_map<Kmer, int> rightCount;
	for (auto& rightKmer : _vertexIndex.getIndexByRead().at(rightRead))
	{
		rightCount[rightKmer.kmer] += 1;
	}

	std::vector<std::pair<int32_t, int32_t>> sharedKmers;
	for (auto& leftKmer : _vertexIndex.getIndexByRead().at(leftRead))
	{
		for (auto& rightKmer : _vertexIndex.getIndexByKmer().at(leftKmer.kmer))
		{
			if (rightKmer.readId == rightRead &&
				leftKmer.position > prevSwitch &&
				leftCount[leftKmer.kmer] == 1 &&
				rightCount[leftKmer.kmer] == 1)
			{
				sharedKmers.push_back({leftKmer.position, 
									   rightKmer.position});
			}
		}
	}

	//filter possible outliers
	std::sort(sharedKmers.begin(), sharedKmers.end(), 
			  [](const std::pair<int32_t, int32_t>& a, 
				 const std::pair<int32_t, int32_t>& b)
				 {return a.first - a.second < b.first - b.second;});
	size_t leftQ = sharedKmers.size() / 4;
	size_t rightQ = sharedKmers.size() * 3 / 4;

	if (leftQ >= rightQ)
	{
		Logger::get().warning() << "No jump found! " +
					  _seqContainer.getIndex().at(leftRead).description +
					  " : " + _seqContainer.getIndex().at(rightRead).description;
		return {prevSwitch + 1, 0};
	}

	std::pair<int32_t, int32_t> bestLeft = {0, 0};
	bool found = false;
	for (size_t i = leftQ; i < rightQ; ++i)
	{
		if (!found || bestLeft.first > sharedKmers[i].first)
		{
			found = true;
			bestLeft = sharedKmers[i];
		}
	}

	return bestLeft;
}
