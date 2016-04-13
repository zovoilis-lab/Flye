#include "Alignment.h"
#include <chrono>


Alignment::Alignment(size_t size):
	_forwardScores(size),
	_reverseScores(size)
{ 
	//ofstream file;
	//file.open("test.txt");
	//std::chrono::time_point<std::chrono::system_clock> now;
	//now = std::chrono::system_clock::now();
	//std::time_t time = std::chrono::system_clock::to_time_t(now);
	//file << "File was produced at: " << std::ctime(&time);
	//file << "\n";
}

double Alignment::globalAlignment(const std::string& v, const std::string& w,
								  ScoringMatrix* sm, int index) 
{
	unsigned int x = v.size() + 1;
	unsigned int y = w.size() + 1;
	FloatMatrix backtrack(x, y, 0);
	FloatMatrix scoreMat(x, y, 0);
		
	double score = this->getBacktrackMatrix(v, w, sm, backtrack, scoreMat);
	_forwardScores[index] = scoreMat;

	//std::string string1;
	//std::string string2;
	//this->traceback(backtrack, v, w, string1, string2);
	//writeStringsToFile(string1, string2, score);
	//writeMatToFile(scoreMat);
	//---------------------------------------------
	//The reverse alignment returns the same score
	//---------------------------------------------
	std::string rev_v = v;
	std::string rev_w = w;
	std::reverse(rev_v.begin(), rev_v.end());
	std::reverse(rev_w.begin(), rev_w.end());

	FloatMatrix backtrackRev(x, y, 0);
	FloatMatrix scoreMatRev(x, y, 0);
	this->getBacktrackMatrix(rev_v, rev_w, sm, backtrackRev, scoreMatRev);
	_reverseScores[index] = scoreMatRev;

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
								  char base, const std::string& read,
								  ScoringMatrix* sm) 
{
	//LetterIndex must start with 1 and go until (row.size - 1)
	const FloatMatrix& forwardScore = _forwardScores[wordIndex];
	const FloatMatrix& reverseScore = _reverseScores[wordIndex];

	size_t frontRow = letterIndex - 1;
	size_t revRow = reverseScore.nrows() - 1 - letterIndex;

	std::vector<double> sub(read.size() + 1);
	sub[0] = forwardScore.at(frontRow, 0) + sm->getScore(base, '-');
	for (size_t i = 0; i < read.size(); ++i)
	{
		double match = forwardScore.at(frontRow, i) + sm->getScore(base, read[i]);
		double ins = forwardScore.at(frontRow, i + 1) + sm->getScore(base, '-');
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
							   char base, const std::string& read,
							   ScoringMatrix* sm) 
{
	//LetterIndex must start with 1 and go until (row.size - 1)
	//using namespace arma;
	const FloatMatrix& forwardScore = _forwardScores[wordIndex];
	const FloatMatrix& reverseScore = _reverseScores[wordIndex];

	size_t frontRow = pos - 1;
	size_t revRow = reverseScore.nrows() - pos;

	std::vector<double> sub(read.size() + 1);
	sub[0] = forwardScore.at(frontRow, 0) + sm->getScore(base, '-');
	for (size_t i = 0; i < read.size(); ++i)
	{
		double match = forwardScore.at(frontRow, i) + sm->getScore(base, read[i]);
		double ins = forwardScore.at(frontRow, i + 1) + sm->getScore(base, '-');
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
									 ScoringMatrix* sm, FloatMatrix& backtrack,
									 FloatMatrix& scoreMat) 
{
	//using namespace std;
	double score = 0.0f;
	
	for (size_t i = 0; i < v.size(); i++) 
	{
		double score = sm->getScore(v[i], '-');
		scoreMat.at(i + 1, 0) = scoreMat.at(i, 0) + score;
		backtrack.at(i + 1, 0) = 1; //Del
	}


	for (size_t i = 0; i < w.size(); i++) {
		double score = sm->getScore('-', w[i]);
		scoreMat.at(0, i + 1) = scoreMat.at(0, i) + score;
		backtrack.at(0, i + 1) = -1; //Ins
	}


	for (size_t i = 1; i < v.size() + 1; i++)
	{
		char key1 = v[i - 1];
		for (size_t j = 1; j < w.size() + 1; j++) 
		{
			char key2 = w[j - 1];

			double left = scoreMat.at(i, j - 1) + sm->getScore('-', key2);
			double up = scoreMat.at(i - 1, j) + sm->getScore(key1, '-');
			score = std::max(left, up);

			double cross = scoreMat.at(i - 1, j - 1) + sm->getScore(key1, key2);
			score = std::max(score, cross);
			scoreMat.at(i, j) = score;

			/*
			if (score == cross) {
				if (v[i - 1] == w[j - 1]) {
					backtrack(i, j) = 7; //Mat
				}
				else {
					backtrack(i, j) = 0; //Sub
				}
			}
			else if (score == up) {
				backtrack(i, j) = 1; //Del
			}
			else if (score == left) {
				backtrack(i, j) = -1; //ins
			}
			*/
		}
	}

	return score;
}


void Alignment::traceback(FloatMatrix& backtrack, const std::string& v, 
						  const std::string& w, std::string& o_v, 
						  std::string& o_w) {
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


Alignment::~Alignment() { }

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

	file << "\n"; 
	file << "\n";

	file.close();
}

void Alignment::writeStringsToFile(const std::string& v, const std::string& w, 
								   float score) 
{
	std::ofstream file;
	file.open("test.txt", std::ios::app);
	file << "--------------------------------\n";
	file << "Score: " << score << "\n";
	file <<  v  << "\n";
	file << w << "\n";
	file << "\n";
	file.close();
}

void Alignment::clean() 
{
	_forwardScores.clear();
	_reverseScores.clear();
}
