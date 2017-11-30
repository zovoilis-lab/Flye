//(c) 2016 by Authors
//This file is a part of ABruijn program.
//Released under the BSD license (see LICENSE file)

#pragma once

#include <string>
#include <vector>
#include <algorithm>
#include <iomanip>
#include <stdexcept>

#include "../common/matrix.h"
#include "subs_matrix.h"


class Alignment 
{

public:
	Alignment(size_t size, const SubstitutionMatrix& sm);

	typedef Matrix<AlnScoreType> ScoreMatrix;

	AlnScoreType globalAlignment(const std::string& v, const std::string& w, 
						   int index);
	AlnScoreType addDeletion(unsigned int index, unsigned int letterIndex);
	AlnScoreType addSubstitution(unsigned int wordIndex, 
						   unsigned int letterIndex,
						   char base, const std::string& read);

	AlnScoreType addInsertion(unsigned int wordIndex,
						   unsigned int positionIndex,
						   char base, const std::string& read);

private:
	std::vector<ScoreMatrix> _forwardScores;
	std::vector<ScoreMatrix> _reverseScores;
	const SubstitutionMatrix& _subsMatrix;

	AlnScoreType getBacktrackMatrix(const std::string& v, const std::string& w,
							     ScoreMatrix& scoreMat);

	void traceback(ScoreMatrix& backtrack, const std::string& v, 
				   const std::string& w, std::string& o_v, 
				   std::string& o_w);
};
