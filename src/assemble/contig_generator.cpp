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


void ContigGenerator::generateContigs()
{
	Logger::get().info() << "Generating contig sequences";

	for (const ContigPath& path : _extender.getContigPaths())
	{
		if (path.reads.size() < 2) 
		{
			if (!path.reads.empty())
			{
				FastaRecord rec(_seqContainer.getSeq(path.reads.front()),
								"read_" + std::to_string(_contigs.size()) + 
									"_part_0_" + 
									_seqContainer.seqName(path.reads.front()),
								FastaRecord::Id(_contigs.size()));
				_contigs.push_back({rec});
			}
			continue;
		}
		

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
	auto alignments = this->generateAlignments(path);
	std::vector<FastaRecord> contigParts;

	auto prevSwitch = std::make_pair(0, 0);
	int partNum = 1;
	for (size_t i = 0; i < path.reads.size(); ++i)
	{
		int32_t leftCut = prevSwitch.second;
		int32_t rightCut = _seqContainer.seqLen(path.reads[i]);
		if (i != path.reads.size() - 1)
		{
			auto curSwitch = this->getSwitchPositions(alignments[i], 
													  prevSwitch.second);
			rightCut = curSwitch.first;
			prevSwitch = curSwitch;
		}

		auto partSeq = _seqContainer.getSeq(path.reads[i])
								.substr(leftCut, rightCut - leftCut);
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
	auto circPath = path;
	circPath.reads.push_back(path.reads[0]);
	circPath.reads.push_back(path.reads[1]);

	auto alignments = this->generateAlignments(circPath);
	int32_t initPivot = _seqContainer.seqLen(circPath.reads[0]) / 2;
	auto firstSwitch = this->getSwitchPositions(alignments[0], 
												initPivot);

	auto prevSwitch = firstSwitch;
	int partNum = 1;
	for (size_t i = 1; i < circPath.reads.size() - 1; ++i)
	{
		auto curSwitch = this->getSwitchPositions(alignments[i], 
												  prevSwitch.second);
		int32_t leftCut = prevSwitch.second;
		int32_t rightCut = curSwitch.first;
		prevSwitch = curSwitch;

		if (i == circPath.reads.size() - 2)	//finishing circle
		{
			rightCut = firstSwitch.first;
			if (rightCut - leftCut <= 0) rightCut = leftCut;
		}

		auto partSeq = _seqContainer.getSeq(circPath.reads[i])
							 	.substr(leftCut, rightCut - leftCut);
		std::string partName = 
			(path.circular ? "circular_" : "linear_") + 
			std::to_string(_contigs.size()) +
			"_part_" + std::to_string(partNum) + "_" +
			_seqContainer.getIndex().at(circPath.reads[i]).description + 
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


std::vector<ContigGenerator::AlignmentInfo> 
ContigGenerator::generateAlignments(const ContigPath& path)
{
	typedef std::tuple<FastaRecord::Id, FastaRecord::Id, size_t> AlnTask;

	std::vector<AlignmentInfo> alnResults;
	std::function<void(const AlnTask&)> alnFunc =
	[this, &alnResults](const AlnTask& aln)
	{
		OverlapRange readsOvlp;
		bool found = false;
		FastaRecord::Id idLeft = std::get<0>(aln);
		FastaRecord::Id idRight = std::get<1>(aln);
		for (auto& ovlp : _overlapContainer.lazySeqOverlaps(idLeft))
		{
			if (ovlp.extId == idRight) 
			{
				readsOvlp = ovlp;
				found = true;
				break;
			}
		}
		if (!found) 
		{
			Logger::get().debug() << _seqContainer.seqName(idLeft) << " " 
				<< _seqContainer.seqName(idRight);
			throw std::runtime_error("Ovlp not found!");
		}

		std::string leftSeq = _seqContainer.getIndex().at(idLeft)
									.sequence.substr(readsOvlp.curBegin,
													 readsOvlp.curRange()).str();
		std::string rightSeq = _seqContainer.getIndex().at(idRight)
									.sequence.substr(readsOvlp.extBegin, 
													 readsOvlp.extRange()).str();

		const int bandWidth = abs((int)leftSeq.length() - 
								  (int)rightSeq.length()) + 
								  		Constants::maximumJump;
		std::string alignedLeft;
		std::string alignedRight;
		pairwiseAlignment(leftSeq, rightSeq, alignedLeft, 
						  alignedRight, bandWidth);

		alnResults[std::get<2>(aln)] = {alignedLeft, alignedRight, 
							  		    readsOvlp.curBegin, 
										readsOvlp.extBegin};
	};

	std::vector<AlnTask> tasks;
	for (size_t i = 0; i < path.reads.size() - 1; ++i)
	{
		tasks.push_back(std::make_tuple(path.reads[i], path.reads[i + 1], i));
	}
	alnResults.resize(tasks.size());
	processInParallel(tasks, alnFunc, Parameters::get().numThreads, false);

	return alnResults;
}


std::pair<int32_t, int32_t> 
ContigGenerator::getSwitchPositions(AlignmentInfo aln,
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
	//			<< _seqContainer.seqName(leftRead) << " : "
	//			<< _seqContainer.seqName(rightRead);
	return {prevSwitch + 1, 0};
}
