#pragma once
#include <string>
#include "ScoringMatrix.h"
#include <vector>
#include <algorithm>
#include <iomanip>
#include <armadillo>

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

private:
	std::vector<arma::mat> forwardScores;
	std::vector<arma::mat> reverseScores;

	//Global alignment help funcs:
	double getBacktrackMatrix(const std::string& v, const std::string& w, 
							  ScoringMatrix* sm, arma::mat& backtrack,
							  arma::mat& scoreMat);

	void traceback(arma::mat& backtrack, const std::string& v, 
				   const std::string& w, std::string& o_v, 
				   std::string& o_w);

	//Auxiliary funcs:
	void writeMatToFile(arma::mat& scoreMat);
	void writeStringsToFile(const std::string& v, const std::string& w, 
							float score);
};

