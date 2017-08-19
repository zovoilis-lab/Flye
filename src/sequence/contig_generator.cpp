//(c) 2016 by Authors
//This file is a part of ABruijn program.
//Released under the BSD license (see LICENSE file)

#include <cassert>
#include <algorithm>
#include <fstream>

#include "contig_generator.h"

#include "../common/config.h"
#include "../common/logger.h"
#include "../common/parallel.h"
#include "../common/matrix.h"


namespace
{
	//banded glocal alignment
	void pairwiseAlignment(const std::string& seqOne, const std::string& seqTwo,
						   std::string& outOne, std::string& outTwo, 
						   int bandWidth)
	{
		static const int32_t MATCH = 5;
		static const int32_t SUBST = -5;
		static const int32_t INDEL = -3;

		static const int32_t matchScore[] = {SUBST, MATCH};
		static const int32_t indelScore[] = {INDEL, 0};

		const size_t numRows = seqOne.length() + 1;
		const size_t numCols = 2 * bandWidth + 1;

		Matrix<char> backtrackMatrix(numRows, numCols);
		std::vector<int32_t> scoreRowOne(numCols, 0);
		std::vector<int32_t> scoreRowTwo(numCols, 0);


		for (size_t i = 0; i < numRows; ++i) 
		{
			size_t j = std::max(0, bandWidth - (int)i);
			backtrackMatrix.at(i, j) = 1;
		}
		for (size_t i = 0; i < numCols; ++i) 
		{
			backtrackMatrix.at(0, i) = 0;
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
				int32_t cross = scoreRowOne[j] + 
								matchScore[seqOne[i - 1] == seqTwo[twoCoord - 1]];
				char maxStep = 2;
				int32_t maxScore = cross;

				if (j < (int)numCols - 1) //up
				{
					int32_t up = scoreRowOne[j + 1] +
								 indelScore[twoCoord == seqTwo.length()];
					if (up > maxScore)
					{
						maxStep = 1;
						maxScore = up;
					}
				}

				if (j > 0) //left
				{
					int32_t left = scoreRowTwo[j - 1] + 
								   indelScore[i == seqOne.length()];
					if (left > maxScore)
					{
						maxStep = 0;
						maxScore = left;
					}
				}

				scoreRowTwo[j] = maxScore;
				backtrackMatrix.at(i, j) = maxStep;
			}
			scoreRowOne.swap(scoreRowTwo);
		}

		//backtrack
		outOne.reserve(seqOne.length() * 9 / 8);
		outTwo.reserve(seqTwo.length() * 9 / 8);

		int i = numRows - 1;
		int j = bandWidth - (int)numRows + (int)seqTwo.length() + 1;
		while (i != 0 || j != bandWidth) 
		{
			size_t twoCoord = j + i - bandWidth;
			if(backtrackMatrix.at(i, j) == 1) //up
			{
				outOne += seqOne[i - 1];
				outTwo += '-';
				i -= 1;
				j += 1;
			}
			else if (backtrackMatrix.at(i, j) == 0) //left
			{
				outOne += '-';
				outTwo += seqTwo[twoCoord - 1];
				j -= 1;
			}
			else	//cross
			{
				outOne += seqOne[i - 1];
				outTwo += seqTwo[twoCoord - 1];
				i -= 1;
			}
		}
		std::reverse(outOne.begin(), outOne.end());
		std::reverse(outTwo.begin(), outTwo.end());
	}
}


void ContigGenerator::generateContigs(const std::vector<ContigPath>& contigs)
{
	Logger::get().info() << "Generating contig sequences";

	for (const ContigPath& path : contigs)
	{
		if (path.sequences.size() < 2) 
		{
			if (!path.sequences.empty())
			{
				FastaRecord rec(path.sequences.front(), path.name, 
								FastaRecord::ID_NONE);
				_contigs.push_back({rec});
			}
			continue;
		}

		std::vector<FastaRecord> contigParts;
		_contigs.push_back(this->generateLinear(path));
	}
}

FastaRecord ContigGenerator::generateLinear(const ContigPath& path)
{
	auto alignments = this->generateAlignments(path);
	std::vector<FastaRecord> contigParts;

	auto prevSwitch = std::make_pair(0, 0);
	std::string contigSequence;
	for (size_t i = 0; i < path.sequences.size(); ++i)
	{
		auto& sequence = path.sequences[i];
		int32_t leftCut = prevSwitch.second;
		int32_t rightCut = sequence.length();
		if (i != path.sequences.size() - 1)
		{
			auto curSwitch = this->getSwitchPositions(alignments[i], 
													  prevSwitch.second);
			rightCut = curSwitch.first;
			prevSwitch = curSwitch;
		}

		contigSequence += sequence.substr(leftCut, rightCut - leftCut).str();
	}
	return FastaRecord(DnaSequence(contigSequence), path.name, 
					   FastaRecord::ID_NONE);
}


void ContigGenerator::outputContigs(const std::string& fileName)
{
	SequenceContainer::writeFasta(_contigs, fileName);
}


std::vector<ContigGenerator::AlignmentInfo> 
ContigGenerator::generateAlignments(const ContigPath& path)
{
	typedef size_t AlnTask;

	std::vector<AlignmentInfo> alnResults;
	std::function<void(const AlnTask&)> alnFunc =
	[this, &alnResults, &path](const AlnTask& i)
	{
		int32_t leftStart = path.overlaps[i].curBegin;
		int32_t leftLen = path.overlaps[i].curRange();
		std::string leftSeq = path.sequences[i].substr(leftStart, leftLen).str();

		int32_t rightStart = path.overlaps[i].extBegin;
		int32_t rightLen = path.overlaps[i].extRange();
		std::string rightSeq = path.sequences[i + 1].substr(rightStart, 
															rightLen).str();
		
		const int bandWidth = abs((int)leftSeq.length() - 
								  (int)rightSeq.length()) + 
								  		Constants::maximumJump;
		if (abs((int)leftSeq.length() - (int)rightSeq.length()) >
			std::min(leftSeq.length(), rightSeq.length()))
		{
			Logger::get().warning() << "Aligning sequence that are too "
				<< " different - something is terribly wrong!";
		}

		std::string alignedLeft;
		std::string alignedRight;
		pairwiseAlignment(leftSeq, rightSeq, alignedLeft, 
						  alignedRight, bandWidth);

		alnResults[i] = {alignedLeft, alignedRight, 
						 leftStart, rightStart};
	};

	std::vector<AlnTask> tasks;
	for (size_t i = 0; i < path.sequences.size() - 1; ++i)
	{
		tasks.push_back(i);
	}
	alnResults.resize(tasks.size());
	processInParallel(tasks, alnFunc, Parameters::get().numThreads, false);

	return alnResults;
}


std::pair<int32_t, int32_t> 
ContigGenerator::getSwitchPositions(const AlignmentInfo& aln,
									int32_t prevSwitch)
{
	int leftPos = aln.startOne;
	int rightPos = aln.startTwo;
	int matchRun = 0;
	for (size_t i = 0; i < aln.alnOne.length(); ++i)
	{
		if (aln.alnOne[i] != '-') ++leftPos;
		if (aln.alnTwo[i] != '-') ++rightPos;

		if (aln.alnOne[i] == aln.alnTwo[i] &&
			leftPos > prevSwitch + Constants::maximumJump)
		{
			++matchRun;
		}
		else
		{
			matchRun = 0;
		}
		if (matchRun == (int)Parameters::get().kmerSize)
		{
			return {leftPos, rightPos};
		}
	}

	Logger::get().debug() << "No jump found!";
	prevSwitch = std::max(prevSwitch + 1, aln.startOne);
	return {prevSwitch, aln.startTwo};
}
