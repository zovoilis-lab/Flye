//(c) 2016 by Authors
//This file is a part of ABruijn program.
//Released under the BSD license (see LICENSE file)

#pragma once

#include <string>
#include <vector>
#include <algorithm>
#include <iomanip>
#include <stdexcept>

#include "Matrix.h"
#include "ScoringMatrix.h"


class Alignment 
{

public:
	Alignment(size_t size);
	double globalAlignment(const std::string& v, const std::string& w, 
						   ScoringMatrix* sm, int index);
	double addDeletion(unsigned int index, unsigned int letterIndex);
	double addSubstitution(unsigned int wordIndex, 
		unsigned int letterIndex,
		char base,
		const std::string& read,
		ScoringMatrix* sm);

	double addInsertion(unsigned int wordIndex,
						unsigned int positionIndex,
						char base, const std::string& read,
						ScoringMatrix* sm);
	void clean();
	~Alignment();

	typedef Matrix<double> FloatMatrix;

private:
	std::vector<FloatMatrix> _forwardScores;
	std::vector<FloatMatrix> _reverseScores;

	//Global alignment help funcs:
	double getBacktrackMatrix(const std::string& v, const std::string& w, 
							  ScoringMatrix* sm,
							  FloatMatrix& scoreMat);

	void traceback(FloatMatrix& backtrack, const std::string& v, 
				   const std::string& w, std::string& o_v, 
				   std::string& o_w);

	//Auxiliary funcs:
	void writeMatToFile(const FloatMatrix& scoreMat);
	void writeStringsToFile(const std::string& v, const std::string& w, 
							float score);
};

