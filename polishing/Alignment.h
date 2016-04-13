#pragma once
#include <string>
#include <vector>
#include <algorithm>
#include <iomanip>
#include <stdexcept>
//#include <armadillo>

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

	typedef Matrix<double> IntMatrix;

private:
	std::vector<IntMatrix> _forwardScores;
	std::vector<IntMatrix> _reverseScores;

	//Global alignment help funcs:
	double getBacktrackMatrix(const std::string& v, const std::string& w, 
							  ScoringMatrix* sm, IntMatrix& backtrack,
							  IntMatrix& scoreMat);

	void traceback(IntMatrix& backtrack, const std::string& v, 
				   const std::string& w, std::string& o_v, 
				   std::string& o_w);

	//Auxiliary funcs:
	void writeMatToFile(const IntMatrix& scoreMat);
	void writeStringsToFile(const std::string& v, const std::string& w, 
							float score);
};

