//(c) 2016 by Authors
//This file is a part of ABruijn program.
//Released under the BSD license (see LICENSE file)

#include <algorithm>

#include "homo_polisher.h"
#include "matrix.h"


namespace
{
	//Computes global pairwise alignment with custom substitution matrix
	void pairwiseAlignment(const std::string& seqOne, const std::string& seqTwo,
						   const SubstitutionMatrix& subsMat,
						   std::string& outOne, std::string& outTwo)
	{
		Matrix<float> scoreMat(seqOne.length() + 1, seqTwo.length() + 1);
		Matrix<char> backtrackMat(seqOne.length() + 1, seqTwo.length() + 1);

		scoreMat.at(0, 0) = 0.0f;
		backtrackMat.at(0, 0) = 0;
		for (size_t i = 0; i < seqOne.length(); ++i) 
		{
			scoreMat.at(i + 1, 0) = scoreMat.at(i, 0) + 
									subsMat.getScore(seqOne[i], '-');
			backtrackMat.at(i + 1, 0) = 1;
		}
		for (size_t i = 0; i < seqTwo.length(); ++i) 
		{
			scoreMat.at(0, i + 1) = scoreMat.at(0, i) + 
									subsMat.getScore('-', seqTwo[i]);
			backtrackMat.at(0, i + 1) = 0;
		}

		//filling DP matrices
		for (size_t i = 1; i < seqOne.length() + 1; ++i)
		{
			for (size_t j = 1; j < seqTwo.length() + 1; ++j) 
			{
				double left = scoreMat.at(i, j - 1) + 
								subsMat.getScore('-', seqTwo[j - 1]);
				double up = scoreMat.at(i - 1, j) + 
								subsMat.getScore(seqOne[i - 1], '-');
				double cross = scoreMat.at(i - 1, j - 1) + 
								subsMat.getScore(seqOne[i - 1], seqTwo[j - 1]);

				int prev = 0;
				float score = left;
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
		outOne.clear();
		outTwo.clear();

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
		outOne += "$";
		outTwo += "$";
	}

	//Splits aligned strings into homopolymer runs (wrt to candAln)
	std::vector<std::pair<HopoMatrix::State, HopoMatrix::Observation>>
	splitBranchHopos(const std::string& candAln, const std::string& branchAln)
	{
		std::vector<std::pair<HopoMatrix::State, 
							  HopoMatrix::Observation>> result;
		size_t prevPos = 0;
		while (candAln[prevPos] == '-') ++prevPos;
		char prevNucl = candAln[prevPos];

		int runLength = 0;
		int gapLength = 0;
		for (size_t pos = prevPos + 1; pos < candAln.length(); ++pos)
		{
			if (candAln[pos] != '-') ++runLength;
			if (candAln[pos] == '-')
			{
				++gapLength;
			}
			else 
			{
				if (candAln[pos] != prevNucl)
				{
					//std::string branchSubseq(branchAln, prevPos, pos - prevPos);
					//std::string candSubseq(candAln, prevPos, pos - prevPos);
					//std::cout << candSubseq << " " << std::endl 
					//		  << branchSubseq << std::endl << std::endl;

					auto state = HopoMatrix::State(candAln, prevPos, pos);
					auto observ = HopoMatrix::strToObs(branchAln, prevPos, 
													   pos);
					result.push_back(std::make_pair(state, observ));

					prevNucl = candAln[pos];
					prevPos = pos - gapLength;
					gapLength = 0;
					runLength = 0;
				}
				else
				{
					gapLength = 0;
				}
			}
			
		}

		//std::cout << "----\n\n";
		return result;
	}
}

void HomoPolisher::polishBubble(Bubble& bubble)
{
	std::string prevCandidate;
	std::string curCandidate = bubble.candidate;

	std::vector<HopoMatrix::State> states;
	std::vector<std::vector<HopoMatrix::Observation>> observations;

	for (auto& branch : bubble.branches)
	{
		std::string alnCand;
		std::string alnBranch;
		pairwiseAlignment(bubble.candidate, branch, _subsMatrix,
						  alnCand, alnBranch);

		//std::cout << alnCand << std::endl << alnBranch << std::endl << std::endl;
		auto splitHopo = splitBranchHopos(alnCand, alnBranch);
		if (states.empty())
		{
			states.assign(splitHopo.size(), HopoMatrix::State());
			observations.assign(splitHopo.size(), 
								std::vector<HopoMatrix::Observation>());
		}
		assert(states.size() == splitHopo.size());

		for (size_t i = 0; i < splitHopo.size(); ++i)
		{
			states[i] = splitHopo[i].first;
			observations[i].push_back(splitHopo[i].second);
		}
	}

	std::string newConsensus;
	for (size_t i = 0; i < states.size(); ++i)
	{
		size_t length = states[i].length;
		if (length > 1)	//only homopolymers
		{
			length = this->mostLikelyLen(states[i], observations[i]);
		}
		newConsensus += std::string(length, states[i].nucl);
		/*if (len != (size_t)states[i].length)
		{

			std::cout << (int)states[i].length << states[i].nucl 
					  << " -> " << len << std::endl;
			for (auto obs : observations[i]) 
				std::cout << HopoMatrix::obsToStr(obs) << std::endl;
			std::cout << std::endl;
		}*/
	}

	if (newConsensus != bubble.candidate)
	{
		StepInfo info;
		info.methodUsed = StepHopo;
		info.sequence = newConsensus;
		bubble.polishSteps.push_back(info);
		bubble.candidate = newConsensus;
	}
}


size_t HomoPolisher::mostLikelyLen(HopoMatrix::State state,
								   const std::vector<HopoMatrix::Observation>&
								   								 observations)
{
	assert(!observations.empty());

	double maxScore = std::numeric_limits<double>::lowest();
	size_t maxRun = 0;
	for (size_t len = 1; len <= 10; ++len)
	{
		double likelihood = 0.0f;
		auto newState = HopoMatrix::State(state.nucl, len);
		for (auto obs : observations)
		{
			likelihood += _hopoMatrix.getScore(newState, obs);
		}
		if (likelihood > maxScore)
		{
			maxScore = likelihood;
			maxRun = len;
		}
		//std::cout << likelihood << " ";
	}
	//std::cout << std::endl << maxScore << std::endl;
	return maxRun;
}
