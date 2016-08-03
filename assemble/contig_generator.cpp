//(c) 2016 by Authors
//This file is a part of ABruijn program.
//Released under the BSD license (see LICENSE file)

#include <cassert>
#include <algorithm>
#include <fstream>

#include "contig_generator.h"
#include "logger.h"
#include "config.h"

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

//banded glocal alignment
void ContigGenerator::pairwiseAlignment(const std::string& seqOne, 
										const std::string& seqTwo,
						   				std::string& outOne, 
										std::string& outTwo)
{
	static const int32_t MATCH = 5;
	static const int32_t SUBST = -5;
	static const int32_t INDEL = -3;

	static const int32_t matchScore[] = {SUBST, MATCH};
	static const int32_t indelScore[] = {INDEL, 0};

	const int bandWidth = abs((int)seqOne.length() - 
							  (int)seqTwo.length()) + Constants::maxumumJump;
	const size_t numRows = seqOne.length() + 1;
	const size_t numCols = 2 * bandWidth + 1;

	//reallocate
	if (_scoreMatrix.nrows() < numRows || _scoreMatrix.ncols() < numCols)
	{
		_scoreMatrix = Matrix<int32_t>(numRows, numCols);
		_backtrackMatrix = Matrix<char>(numRows, numCols);
	}


	for (size_t i = 0; i < numRows; ++i) 
	{
		size_t j = std::max(0, bandWidth - (int)i);
		_scoreMatrix.at(i, j) = 0;
		_backtrackMatrix.at(i, j) = 1;
	}
	for (size_t i = 0; i < numCols; ++i) 
	{
		_scoreMatrix.at(0, i) = 0;
		_backtrackMatrix.at(0, i) = 0;
	}

	//filling DP matrices
	for (size_t i = 1; i < numRows; ++i)
	{
		int leftOverhang = bandWidth - (int)i + 1;
		int rightOverhand = (int)i + bandWidth - (int)seqTwo.length();
		size_t colLeft = std::max(0, leftOverhang);
		size_t colRight = std::min((int)numCols, (int)numCols - rightOverhand);

		for (int j = colLeft; j < (int)colRight; ++j)
		{
			size_t twoCoord = j + i - bandWidth;
			int32_t cross = _scoreMatrix.at(i - 1, j) + 
							matchScore[seqOne[i - 1] == seqTwo[twoCoord - 1]];
			char maxStep = 2;
			int32_t maxScore = cross;

			if (j < (int)numCols - 1) //up
			{
				int32_t up = _scoreMatrix.at(i - 1, j + 1) +
								indelScore[twoCoord == seqTwo.length()];
 				if (up > maxScore)
				{
					maxStep = 1;
					maxScore = up;
				}
			}

			if (j > 0) //left
			{
				int32_t left = _scoreMatrix.at(i, j - 1) + 
								indelScore[i == seqOne.length()];
				if (left > maxScore)
				{
					maxStep = 0;
					maxScore = left;
				}
			}

			_scoreMatrix.at(i, j) = maxScore;
			_backtrackMatrix.at(i, j) = maxStep;
		}
	}

	//backtrack
	outOne.reserve(seqOne.length() * 9 / 8);
	outTwo.reserve(seqTwo.length() * 9 / 8);

	int i = numRows - 1;
	int j = bandWidth - (int)numRows + (int)seqTwo.length() + 1;
	while (i != 0 || j != bandWidth) 
	{
		if(_backtrackMatrix.at(i, j) == 1) //up
		{
			outOne += seqOne[i - 1];
			outTwo += '-';
			i -= 1;
			j += 1;
		}
		else if (_backtrackMatrix.at(i, j) == 0) //left
		{
			size_t twoCoord = j + i - bandWidth;
			outOne += '-';
			outTwo += seqTwo[twoCoord - 1];
			j -= 1;
		}
		else	//cross
		{
			size_t twoCoord = j + i - bandWidth;
			outOne += seqOne[i - 1];
			outTwo += seqTwo[twoCoord - 1];
			i -= 1;
		}
	}
	std::reverse(outOne.begin(), outOne.end());
	std::reverse(outTwo.begin(), outTwo.end());

	/*
	std::cerr << seqOne.substr(seqOne.size() - 100) << std::endl 
			  << seqTwo.substr(seqTwo.size() - 100) << std::endl;
	std::cerr << outOne.substr(outOne.size() - 100) << std::endl 
			  << outTwo.substr(outTwo.size() - 100) << std::endl << std::endl;
	*/
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
	pairwiseAlignment(leftSeq, rightSeq, alignedLeft, alignedRight);

	int leftPos = readsOvlp->curBegin;
	int rightPos = readsOvlp->extBegin;
	int matchRun = 0;
	for (size_t i = 0; i < alignedLeft.length(); ++i)
	{
		if (alignedLeft[i] != '-') ++leftPos;
		if (alignedRight[i] != '-') ++rightPos;

		if (alignedLeft[i] == alignedRight[i] &&
			leftPos > prevSwitch + Constants::maxumumJump)
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
