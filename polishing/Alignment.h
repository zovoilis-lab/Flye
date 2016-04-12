#pragma once
#include <string>
#include "ScoringMatrix.h"
#include <vector>
#include <algorithm>
#include <iomanip>
#include <armadillo>

using namespace std;

class Alignment {
	vector<arma::mat> forwardScores;
	vector<arma::mat> reverseScores;

	//Global alignment help funcs:
	double getBacktrackMatrix(
		string v, 
		string w, 
		ScoringMatrix* sm,
		arma::mat& backtrack,
		arma::mat& scoreMat);

	void traceback(arma::mat& backtrack, string v, string w,
		string& o_v, string& o_w);

	//Auxiliary funcs:
	void writeMatToFile(arma::mat& scoreMat);
	void writeStringsToFile(string v, string w, float score);

public:
	Alignment(size_t size);
	double globalAlignment(string v, string w, ScoringMatrix* sm, int index);
	double addDeletion(unsigned int index, unsigned int letterIndex);
	double addSubstitution(unsigned int wordIndex, 
		unsigned int letterIndex,
		char base,
		std::string read,
		ScoringMatrix* sm);

	double addInsertion(unsigned int wordIndex,
		unsigned int positionIndex,
		char base,
		std::string read,
		ScoringMatrix* sm);

	void clean();
	~Alignment();
};

