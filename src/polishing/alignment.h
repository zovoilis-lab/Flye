//(c) 2016 by Authors
//This file is a part of ABruijn program.
//Released under the BSD license (see LICENSE file)

#pragma once

#include <string>
#include <vector>
#include <algorithm>
#include <iomanip>
#include <stdexcept>

#include "matrix.h"
#include "subs_matrix.h"


class Alignment 
{

public:
	Alignment(size_t size, const SubstitutionMatrix& sm);

	double globalAlignment(const std::string& v, const std::string& w, 
						   int index);
	double addDeletion(unsigned int index, unsigned int letterIndex);
	double addSubstitution(unsigned int wordIndex, 
						   unsigned int letterIndex,
						   char base, const std::string& read);

	double addInsertion(unsigned int wordIndex,
						unsigned int positionIndex,
						char base, const std::string& read);

	typedef Matrix<double> FloatMatrix;
private:
	std::vector<FloatMatrix> _forwardScores;
	std::vector<FloatMatrix> _reverseScores;
	const SubstitutionMatrix& _subsMatrix;

	double getBacktrackMatrix(const std::string& v, const std::string& w,
							  FloatMatrix& scoreMat);

	void traceback(FloatMatrix& backtrack, const std::string& v, 
				   const std::string& w, std::string& o_v, 
				   std::string& o_w);
};
