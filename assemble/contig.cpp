//(c) 2016 by Authors
//This file is a part of ABruijn program.
//Released under the BSD license (see LICENSE file)

#include <cassert>
#include <algorithm>
#include <fstream>

#include "contig.h"
#include "logger.h"
#include "matrix.h"

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

namespace
{
	void pairwiseAlignment(const std::string& seqOne, const std::string& seqTwo,
						   std::string& outOne, std::string& outTwo, 
						   int bandWidth)
	{
		static const int32_t MATCH = 2;
		static const int32_t SUBST = -2;
		static const int32_t INDEL = -1;
		auto matchScore = [](char a, char b){return a == b ? MATCH : SUBST;};

		Matrix<int32_t> scoreMat(seqOne.length() + 1, seqTwo.length() + 1);
		Matrix<char> backtrackMat(seqOne.length() + 1, seqTwo.length() + 1);

		scoreMat.at(0, 0) = 0;
		backtrackMat.at(0, 0) = 0;
		for (size_t i = 0; i < seqOne.length(); ++i) 
		{
			scoreMat.at(i + 1, 0) = scoreMat.at(i, 0) + INDEL;
			backtrackMat.at(i + 1, 0) = 1;
		}
		for (size_t i = 0; i < seqTwo.length(); ++i) 
		{
			scoreMat.at(0, i + 1) = scoreMat.at(0, i) + INDEL;
			backtrackMat.at(0, i + 1) = 0;
		}

		//filling DP matrices
		//TODO: more accurate indices
		for (size_t i = 1; i < seqOne.length() + 1; ++i)
		{
			size_t diagLeft = std::max(1, (int)i - bandWidth);
			size_t diagRight = std::min((int)seqTwo.length() + 1, 
										(int)i + bandWidth);
			for (size_t j = diagLeft; j < diagRight; ++j)
			{
				int32_t left = scoreMat.at(i, j - 1) + INDEL;
				int32_t up = scoreMat.at(i - 1, j) + INDEL;
				int32_t cross = scoreMat.at(i - 1, j - 1) + 
								matchScore(seqOne[i - 1], seqTwo[j - 1]);

				int prev = 0;
				int32_t score = left;
				if (up > score)
				{
					prev = 1;
					score = up;
				}
				if (cross > score)
				{
					prev = 2;
					score = cross;
				}
				scoreMat.at(i, j) = score;
				backtrackMat.at(i, j) = prev;
			}
		}

		//backtrack
		int i = seqOne.length();
		int j = seqTwo.length();
		outOne.reserve(i * 1.5);
		outTwo.reserve(j * 1.5);

		while (i != 0 || j != 0) 
		{
			if(backtrackMat.at(i, j) == 1) 
			{
				outOne += seqOne[i - 1];
				outTwo += '-';
				i -= 1;
			}
			else if (backtrackMat.at(i, j) == 0) 
			{
				outOne += '-';
				outTwo += seqTwo[j - 1];
				j -= 1;
			}
			else
			{
				outOne += seqOne[i - 1];
				outTwo += seqTwo[j - 1];
				i -= 1;
				j -= 1;
			}
		}
		std::reverse(outOne.begin(), outOne.end());
		std::reverse(outTwo.begin(), outTwo.end());
	}
}


std::pair<int32_t, int32_t> 
ContigGenerator::getSwitchPositions(FastaRecord::Id leftRead, 
									FastaRecord::Id rightRead,
									int32_t prevSwitch)
{
	//find overlap
	const OverlapRange* readsOvlp = nullptr;
	for (auto& ovlp : _overlapDetector.getOverlapIndex().at(leftRead))
	{
		if (ovlp.extId == rightRead) readsOvlp = &ovlp;
	}

	//Alignment for a precise shift calculation
	std::string leftSeq = _seqContainer.getIndex().at(leftRead)
								.sequence.substr(readsOvlp->curBegin, 
											 	 readsOvlp->curRange());
	std::string rightSeq = _seqContainer.getIndex().at(rightRead)
								.sequence.substr(readsOvlp->extBegin, 
											 	 readsOvlp->extRange());
	std::string alignedLeft;
	std::string alignedRight;
	int width = std::max(abs((int)leftSeq.length() - (int)rightSeq.length()) * 2,
						 _maximumJump);
	pairwiseAlignment(leftSeq, rightSeq, alignedLeft, 
					  alignedRight, width);

	//_vertexIndex.getKmerSize();
	int leftPos = readsOvlp->curBegin;
	int rightPos = readsOvlp->extBegin;
	int matchRun = 0;
	for (size_t i = 0; i < alignedLeft.length(); ++i)
	{
		if (alignedLeft[i] != '-') ++leftPos;
		if (alignedRight[i] != '-') ++rightPos;

		if (alignedLeft[i] == alignedRight[i] &&
			leftPos > prevSwitch + _maximumJump)
		{
			++matchRun;
		}
		else
		{
			matchRun = 0;
		}
		if (matchRun == (int)_vertexIndex.getKmerSize())
		{
			return {leftPos, rightPos};
		}
	}
	
	Logger::get().warning() << "No jump found! "
				<< _seqContainer.seqName(leftRead) << " : "
				<< _seqContainer.seqName(rightRead);
	return {prevSwitch + 1, 0};
}
