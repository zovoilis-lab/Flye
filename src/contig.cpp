#include <cassert>
#include <algorithm>
#include <fstream>

#include "contig.h"

void ContigGenerator::generateContigs()
{
	const auto& readPaths = _extender.getContigPaths();
	for (const Extender::ReadPath& path : readPaths)
	{
		std::vector<FastaRecord> contigParts;
		FastaRecord::ReadIdType partId = 1;
		Extender::ReadPath circPath = path;
		circPath.push_back(path[0]);
		circPath.push_back(path[1]);
		auto prevSwitch = this->getSwitchPositions(circPath[0], circPath[1], 0);

		for (size_t i = 1; i < circPath.size() - 1; ++i)
		{
			auto curSwitch = this->getSwitchPositions(circPath[i - 1], circPath[i], 
													  prevSwitch.second);
			int32_t leftCut = prevSwitch.second;
			int32_t rightCut = curSwitch.first;

			std::string partSeq = _seqContainer.getIndex().at(circPath[i - 1])
								 .sequence.substr(leftCut, rightCut - leftCut);
			std::string partName = 
				"part_" + std::to_string(partId) + "_" +
				_seqContainer.getIndex().at(circPath[i - 1]).description + 
				"[" + std::to_string(leftCut) + ":" + 
				std::to_string(rightCut) + "]";
			contigParts.push_back(FastaRecord(partSeq, partName, partId));
			++partId;
		}
		_contigs.push_back(contigParts);
	}
}


void ContigGenerator::outputContigs(const std::string& fileName)
{
	std::ofstream fout(fileName);
	for (auto& part : _contigs.front())
	{
		fout << ">" << part.description << std::endl 
			 << part.sequence << std::endl;
	}
}


std::pair<int32_t, int32_t> 
ContigGenerator::getSwitchPositions(FastaRecord::ReadIdType leftRead, 
									FastaRecord::ReadIdType rightRead,
									int32_t prevSwitch)
{
	int32_t ovlpShift = -1;
	for (auto& ovlp : _overlapDetector.getOverlapIndex().at(leftRead))
	{
		if (ovlp.extId == rightRead)
		{
			ovlpShift = (ovlp.extBegin - ovlp.curBegin + 
						 ovlp.extEnd - ovlp.curEnd) / 2;
			break;
		}
	}
	assert(ovlpShift > 0);

	std::vector<std::pair<int32_t, int32_t>> acceptedKmers;
	for (auto& leftKmer : _vertexIndex.getIndexByRead().at(leftRead))
	{
		if (leftKmer.position <= prevSwitch)
			continue;

		for (auto& rightKmer : _vertexIndex.getIndexByKmer().at(leftKmer.kmer))
		{
			if (rightKmer.readId != rightRead)
				continue;

			int32_t kmerShift = rightKmer.position - leftKmer.position;
			if (abs(kmerShift - ovlpShift) > _maximumJump / 2)
			{
				acceptedKmers.push_back(std::make_pair(leftKmer.position, 
													   rightKmer.position));
			}
		}
	}
	assert(!shifts.empty());

	//return median by kmer shifts
	std::nth_element(acceptedKmers.begin(), 
					 acceptedKmers.begin() + acceptedKmers.size() / 2,
					 acceptedKmers.end(),
					 [](const std::pair<int32_t, int32_t>& a, 
					 	const std::pair<int32_t, int32_t>& b)
					 		{return a.second - a.first < b.second - b.first;});
	return acceptedKmers[acceptedKmers.size() / 2];
}
