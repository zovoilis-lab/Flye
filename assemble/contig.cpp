//(c) 2016 by Authors
//This file is a part of ABruijn program.
//Released under the BSD license (see LICENSE file)

#include <cassert>
#include <algorithm>
#include <fstream>

#include "contig.h"
#include "utility.h"

void ContigGenerator::generateContigs()
{
	LOG_PRINT("Generating contig sequences");
	const auto& readPaths = _extender.getContigPaths();
	for (const ContigPath& path : readPaths)
	{
		if (path.reads.size() < 2) continue;

		std::vector<FastaRecord> contigParts;
		int partNum = 1;
		auto circPath = path.reads;
		if (path.circular)
		{
			circPath.push_back(path.reads[0]);
			circPath.push_back(path.reads[1]);
		}

		int initPivot = _seqContainer.getIndex().at(circPath[0])
												.sequence.length() / 2;
		auto firstSwitch = this->getSwitchPositions(circPath[0], circPath[1], 
													initPivot);
		auto prevSwitch = firstSwitch;

		int readShift = 1;
		for (size_t i = 1; i < circPath.size() - 1; ++i)
		{
			auto curSwitch = this->getSwitchPositions(circPath[i], 
													  circPath[i + readShift], 
													  prevSwitch.second);
			//TODO: a better solution
			if (curSwitch.first == prevSwitch.second)
			{
				prevSwitch.second = 0;
				continue;
			}

			int32_t leftCut = prevSwitch.second;
			int32_t rightCut = curSwitch.first;
			prevSwitch = curSwitch;

			if (path.circular && i == circPath.size() - 2)	//finishing circle
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
		_contigs.push_back(contigParts);
	}
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
	/*int32_t ovlpShift = 0;
	for (auto& ovlp : _overlapDetector.getOverlapIndex().at(leftRead))
	{
		if (ovlp.extId == rightRead)
		{
			ovlpShift = (ovlp.curBegin - ovlp.extBegin + 
						 ovlp.curEnd - ovlp.extEnd) / 2;
			break;
		}
	}*/

	std::vector<std::pair<int32_t, int32_t>> acceptedKmers;
	for (auto& leftKmer : _vertexIndex.getIndexByRead().at(leftRead))
	{
		if (leftKmer.position <= prevSwitch)
			continue;

		for (auto& rightKmer : _vertexIndex.getIndexByKmer().at(leftKmer.kmer))
		{
			if (rightKmer.readId != rightRead)
				continue;

			//int32_t kmerShift = leftKmer.position - rightKmer.position;
			//if (abs(kmerShift - ovlpShift) < _maximumJump / 2)
			//{
			acceptedKmers.push_back(std::make_pair(leftKmer.position, 
												   rightKmer.position));
			//}
		}
	}
	if (acceptedKmers.empty())
	{
		WARNING_PRINT("Warning: no jump found! " +
					  _seqContainer.getIndex().at(leftRead).description +
					  " : " + _seqContainer.getIndex().at(rightRead).description);
		return std::make_pair(prevSwitch, prevSwitch);
	}

	//returna median among kmer shifts
	std::nth_element(acceptedKmers.begin(), 
					 acceptedKmers.begin() + acceptedKmers.size() / 2,
					 acceptedKmers.end(),
					 [](const std::pair<int32_t, int32_t>& a, 
					 	const std::pair<int32_t, int32_t>& b)
					 		{return a.first - a.second < b.first - b.second;});
	return acceptedKmers[acceptedKmers.size() / 2];
}
