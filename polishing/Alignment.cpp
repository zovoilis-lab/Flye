#include "Alignment.h"
#include <chrono>

Alignment::Alignment(size_t size){ 
	//ofstream file;
	//file.open("test.txt");
	//std::chrono::time_point<std::chrono::system_clock> now;
	//now = std::chrono::system_clock::now();
	//std::time_t time = std::chrono::system_clock::to_time_t(now);
	//file << "File was produced at: " << std::ctime(&time);
	//file << "\n";

	forwardScores.resize(size);
	reverseScores.resize(size);
}

double Alignment::globalAlignment(std::string v, std::string w, ScoringMatrix* sm, int index) {
	unsigned int x = v.size() + 1;
	unsigned int y = w.size() + 1;
	arma::mat backtrack(x, y);
	backtrack.fill(0);

	arma::mat scoreMat(x, y);
	scoreMat.fill(0);
		
	double score = getBacktrackMatrix(v, w, sm, backtrack, scoreMat);
	

	forwardScores[index] = scoreMat;

	string string1;
	string string2;
	traceback(backtrack, v, w, string1, string2);
	//writeStringsToFile(string1, string2, score);
	//writeMatToFile(scoreMat);
	//---------------------------------------------
	//The reverse alignment returns the same score
	//---------------------------------------------
	string rev_v = v;
	string rev_w = w;
	std::reverse(rev_v.begin(), rev_v.end());
	std::reverse(rev_w.begin(), rev_w.end());

	arma::mat backtrackRev(x, y);
	backtrackRev.fill(0);

	arma::mat scoreMatRev(x, y);
	scoreMatRev.fill(0);

	getBacktrackMatrix(rev_v, rev_w, sm, backtrackRev, scoreMatRev);
	reverseScores[index] = scoreMatRev;

	return(score);
}

double Alignment::addDeletion(unsigned int wordIndex, unsigned int letterIndex) {
	using namespace arma;
	mat forwardScore = forwardScores[wordIndex];
	mat reverseScore = reverseScores[wordIndex];

	//writeMatToFile(forwardScore);
	//writeMatToFile(reverseScore);

	//Forward matrix part
	rowvec front = forwardScore.row(letterIndex-1);

	//Reverse matrix part
	//Note: We subtract 2 because of zero indexing and an extra added row and column count
	unsigned int index = (reverseScore.n_rows - 1) - letterIndex;
	rowvec rev_row = reverseScore.row(index);
	rowvec flipped = fliplr(rev_row);

	rowvec sum = front + flipped;
	return arma::max(sum);
}

double Alignment::addSubstitution(unsigned int wordIndex, 
	unsigned int letterIndex,
	char base,
	std::string read,
	ScoringMatrix* sm) {

	//LetterIndex must start with 1 and go until (row.size - 1)
	using namespace arma;
	mat forwardScore = forwardScores[wordIndex];
	mat reverseScore = reverseScores[wordIndex];

	//Forward matrix part (Grabs a column before where we make a substitution)
	rowvec front = forwardScore.row(letterIndex-1);
	
	//Fill sub
	rowvec sub(read.size()+1);
	sub[0] = front[0] + sm->getScore(base, '-');
	for (int i = 0; i < read.size(); i++) {
		double match = front[i] + sm->getScore(base, read[i]);
		double insertion = front[i+1] + sm->getScore(base, '-');
		sub[i+1] = std::max(match, insertion);
	}

	//Reverse matrix part
	//Note: We subtract 1 because of zero indexing and an extra added by size()
	unsigned int index = (reverseScore.n_rows - 1) - letterIndex;
	rowvec rev_row = reverseScore.row(index);
	rowvec flipped = fliplr(rev_row);

	rowvec sum = sub + flipped;
	return arma::max(sum);
}


double Alignment::addInsertion(unsigned int wordIndex,
	unsigned int pos, 
	char base,
	std::string read,
	ScoringMatrix* sm) {

	//LetterIndex must start with 1 and go until (row.size - 1)
	using namespace arma;
	mat forwardScore = forwardScores[wordIndex];
	mat reverseScore = reverseScores[wordIndex];

	//Forward matrix part (Grabs a column before where we make a substitution)
	rowvec front = forwardScore.row(pos - 1);

	//Fill sub
	rowvec sub(read.size() + 1);
	sub[0] = front[0] + sm->getScore(base, '-');
	for (int i = 0; i < read.size(); i++) {
		double match = front[i] + sm->getScore(base, read[i]);
		double insertion = front[i + 1] + sm->getScore(base, '-');
		sub[i + 1] = std::max(match, insertion);
	}

	//Reverse matrix part
	//Note: We subtract 1 because of zero indexing and an extra added by size()
	unsigned int index = reverseScore.n_rows - pos;
	rowvec rev_row = reverseScore.row(index);
	rowvec flipped = fliplr(rev_row);

	rowvec sum = sub + flipped;
	return arma::max(sum);
}



double Alignment::getBacktrackMatrix(
	std::string v, 
	std::string w, 
	ScoringMatrix* sm, 
	arma::mat& backtrack,
	arma::mat& s) {

	using namespace std;
	double score;
	
	for (int i = 0; i < v.size(); i++) {
		double score = sm->getScore(v[i], '-');
		s(i+1,0) = s(i, 0) + score;
		backtrack(i+1,0) = 1; //Del
	}


	for (int i = 0; i < w.size(); i++) {
		double score = sm->getScore('-', w[i]);
		s(0, i+1) = s(0, i) + score;
		backtrack(0, i+1) = -1; //Ins
	}


	for (int i = 1; i < v.size()+1; i++){
		char key1 = v[i - 1];
		for (int j = 1; j < w.size()+1; j++) {
			char key2 = w[j - 1];

			double left = s(i, j - 1) + sm->getScore('-', key2);
			double up = s(i - 1, j) + sm->getScore(key1, '-');
			score = max(left, up);

			double cross = s(i - 1, j - 1) + sm->getScore(key1, key2);
			score = max(score, cross);
			s(i, j) = score;

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
		}
	}

	return(score);
}


void Alignment::traceback(arma::mat& backtrack, string v, string w,
	string& o_v, string& o_w) {
	int i = v.size();
	int j = w.size();
	o_v = "";
	o_w = "";

	while (i != 0 || j != 0) {
		if(backtrack(i, j) == 1){ //Del
			o_v += v[i - 1];
			o_w += '-';
			i -= 1;
		}
		else if (backtrack(i, j) == -1) { //Ins
			o_v += '-';
			o_w += w[j - 1];
			j -= 1;
		}
		else if (backtrack(i, j) == 7) { //Mat
			o_v += v[i - 1];
			o_w += w[j - 1];
			i -= 1;
			j -= 1;
		}
		else if (backtrack(i, j) == 0) { //Sub
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

void Alignment::writeMatToFile(arma::mat& matrix) {
	ofstream file;
	arma::mat scoreMat = matrix.t();
	file.open("test.txt", ios::app);
	for (int i = 0; i < scoreMat.n_rows; i++) {
		for (int j = 0; j < scoreMat.n_cols; j++) {
			file << fixed;
			file << setw(4) << left << std::setprecision(2) <<scoreMat(i, j)  << "\t";
		}
		file << "\n";
	}

	file << "\n"; 
	file << "\n";

	file.close();
}

void Alignment::writeStringsToFile(string v, string w, float score) {
	ofstream file;
	file.open("test.txt", ios::app);
	file << "--------------------------------\n";
	file << "Score: " << score << "\n";
	file <<  v  << "\n";
	file << w << "\n";
	file << "\n";
	file.close();
}

void Alignment::clean() {
	forwardScores.clear();
	reverseScores.clear();
}
