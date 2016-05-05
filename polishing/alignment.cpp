//(c) 2016 by Authors
//This file is a part of ABruijn program.
//Released under the BSD license (see LICENSE file)

#include "alignment.h"
#include <chrono>


Alignment::Alignment(size_t size, const SubstitutionMatrix& sm):
	_forwardScores(size),
	_reverseScores(size),
	_subsMatrix(sm)
{ 
}

double Alignment::globalAlignment(const std::string& v, const std::string& w,
								  int index) 
{
	unsigned int x = v.size() + 1;
	unsigned int y = w.size() + 1;
	FloatMatrix scoreMat(x, y, 0);
		
	double score = this->getBacktrackMatrix(v, w, scoreMat);
	_forwardScores[index] = std::move(scoreMat);

	//std::string string1;
	//std::string string2;
	//this->traceback(backtrack, v, w, string1, string2);
	//writeStringsToFile(string1, string2, score);
	//writeMatToFile(scoreMat);
	//---------------------------------------------
	//The reverse alignment returns the same score
	//---------------------------------------------
	std::string rev_v(v.rbegin(), v.rend());
	std::string rev_w(w.rbegin(), w.rend());

	FloatMatrix scoreMatRev(x, y, 0);
	this->getBacktrackMatrix(rev_v, rev_w, scoreMatRev);
	_reverseScores[index] = std::move(scoreMatRev);

	return score;
}

double Alignment::addDeletion(unsigned int wordIndex, unsigned int letterIndex) 
{
	const FloatMatrix& forwardScore = _forwardScores[wordIndex];
	const FloatMatrix& reverseScore = _reverseScores[wordIndex];
	//writeMatToFile(forwardScore);
	//writeMatToFile(reverseScore);

	//Note: We subtract 2 because of zero indexing and an extra added row and column count
	//unsigned int index = (reverseScore.nrows() - 1) - letterIndex;
	size_t frontRow = letterIndex - 1;
	size_t revRow = reverseScore.nrows() - 1 - letterIndex;
	
	double maxVal = std::numeric_limits<double>::lowest();
	for (size_t col = 0; col < forwardScore.ncols(); ++col)
	{
		size_t backCol = forwardScore.ncols() - col - 1;
		double sum = forwardScore.at(frontRow, col) + 
				  reverseScore.at(revRow, backCol);
		maxVal = std::max(maxVal, sum);
	}

	return maxVal;
}

double Alignment::addSubstitution(unsigned int wordIndex, 
								  unsigned int letterIndex,
								  char base, const std::string& read) 
{
	//LetterIndex must start with 1 and go until (row.size - 1)
	const FloatMatrix& forwardScore = _forwardScores[wordIndex];
	const FloatMatrix& reverseScore = _reverseScores[wordIndex];

	size_t frontRow = letterIndex - 1;
	size_t revRow = reverseScore.nrows() - 1 - letterIndex;

	std::vector<double> sub(read.size() + 1);
	sub[0] = forwardScore.at(frontRow, 0) + _subsMatrix.getScore(base, '-');
	for (size_t i = 0; i < read.size(); ++i)
	{
		double match = forwardScore.at(frontRow, i) + 
						_subsMatrix.getScore(base, read[i]);
		double ins = forwardScore.at(frontRow, i + 1) + 
						_subsMatrix.getScore(base, '-');
		sub[i + 1] = std::max(match, ins);
	}

	double maxVal = std::numeric_limits<double>::lowest();
	for (size_t col = 0; col < forwardScore.ncols(); ++col)
	{
		size_t backCol = forwardScore.ncols() - col - 1;
		double sum = sub[col] + reverseScore.at(revRow, backCol);
		maxVal = std::max(maxVal, sum);
	}
	return maxVal;
}


double Alignment::addInsertion(unsigned int wordIndex,
							   unsigned int pos, 
							   char base, const std::string& read) 
{
	//LetterIndex must start with 1 and go until (row.size - 1)
	const FloatMatrix& forwardScore = _forwardScores[wordIndex];
	const FloatMatrix& reverseScore = _reverseScores[wordIndex];

	size_t frontRow = pos - 1;
	size_t revRow = reverseScore.nrows() - pos;

	std::vector<double> sub(read.size() + 1);
	sub[0] = forwardScore.at(frontRow, 0) + _subsMatrix.getScore(base, '-');
	for (size_t i = 0; i < read.size(); ++i)
	{
		double match = forwardScore.at(frontRow, i) + 
						_subsMatrix.getScore(base, read[i]);
		double ins = forwardScore.at(frontRow, i + 1) + 
						_subsMatrix.getScore(base, '-');
		sub[i + 1] = std::max(match, ins);
	}

	double maxVal = std::numeric_limits<double>::lowest();
	for (size_t col = 0; col < forwardScore.ncols(); ++col)
	{
		size_t backCol = forwardScore.ncols() - col - 1;
		double sum = sub[col] + reverseScore.at(revRow, backCol);
		maxVal = std::max(maxVal, sum);
	}
	return maxVal;
}


double Alignment::getBacktrackMatrix(const std::string& v, const std::string& w,
									 FloatMatrix& scoreMat) 
{
	double score = 0.0f;
	
	for (size_t i = 0; i < v.size(); i++) 
	{
		double score = _subsMatrix.getScore(v[i], '-');
		scoreMat.at(i + 1, 0) = scoreMat.at(i, 0) + score;
	}


	for (size_t i = 0; i < w.size(); i++) {
		double score = _subsMatrix.getScore('-', w[i]);
		scoreMat.at(0, i + 1) = scoreMat.at(0, i) + score;
	}


	for (size_t i = 1; i < v.size() + 1; i++)
	{
		char key1 = v[i - 1];
		for (size_t j = 1; j < w.size() + 1; j++) 
		{
			char key2 = w[j - 1];

			double left = scoreMat.at(i, j - 1) + 
							_subsMatrix.getScore('-', key2);
			double up = scoreMat.at(i - 1, j) + 
							_subsMatrix.getScore(key1, '-');
			score = std::max(left, up);

			double cross = scoreMat.at(i - 1, j - 1) + 
							_subsMatrix.getScore(key1, key2);
			score = std::max(score, cross);
			scoreMat.at(i, j) = score;
		}
	}

	return score;
}


void Alignment::traceback(FloatMatrix& backtrack, const std::string& v, 
						  const std::string& w, std::string& o_v, 
						  std::string& o_w) 
{
	int i = v.size();
	int j = w.size();
	o_v = "";
	o_w = "";

	while (i != 0 || j != 0) {
		if(backtrack.at(i, j) == 1){ //Del
			o_v += v[i - 1];
			o_w += '-';
			i -= 1;
		}
		else if (backtrack.at(i, j) == -1) { //Ins
			o_v += '-';
			o_w += w[j - 1];
			j -= 1;
		}
		else if (backtrack.at(i, j) == 7) { //Mat
			o_v += v[i - 1];
			o_w += w[j - 1];
			i -= 1;
			j -= 1;
		}
		else if (backtrack.at(i, j) == 0) { //Sub
			o_v += v[i - 1];
			o_w += w[j - 1];
			i -= 1;
			j -= 1;
		}
	}
	std::reverse(o_v.begin(), o_v.end());
	std::reverse(o_w.begin(), o_w.end());
}

void Alignment::writeMatToFile(const FloatMatrix& matrix) 
{
	std::ofstream file;
	//arma::mat scoreMat = matrix.t();
	const FloatMatrix& scoreMat = matrix;
	file.open("test.txt", std::ios::app);
	for (size_t i = 0; i < scoreMat.nrows(); i++) {
		for (size_t j = 0; j < scoreMat.ncols(); j++) {
			file << std::fixed;
			file << std::setw(4) << std::left << std::setprecision(2) 
				 << scoreMat.at(i, j) << "\t";
		}
		file << "\n";
	}
}

void Alignment::writeStringsToFile(const std::string& v, const std::string& w, 
								   float score) 
{
	std::ofstream file;
	file.open("test.txt", std::ios::app);
	file << "--------------------------------\n";
	file << "Score: " << score << "\n";
	file <<  v  << "\n";
	file << w;
}
