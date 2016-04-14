#include "ScoringMatrix.h"
#include <sstream>
#include <vector>
#include <cmath>
#include <stdexcept>


ScoringMatrix::ScoringMatrix(int xrange, int yrange):
	_xrange(xrange), _yrange(yrange), _transTable(nullptr)
{	
	_matrix = new double*[xrange];
	for (int i = 0; i < xrange; i++) {
		_matrix[i] = new double[yrange];
	}
	_matrix[_xrange - 1][_yrange - 1] = 0;

	_transTable = new int[256];	//256 chars
	for (size_t i = 0; i < 256; ++i)
		_transTable[i] = -1;
	_transTable[(size_t)'A'] = 0;
	_transTable[(size_t)'a'] = 0;
	_transTable[(size_t)'C'] = 1;
	_transTable[(size_t)'c'] = 1;
	_transTable[(size_t)'T'] = 2;
	_transTable[(size_t)'t'] = 2;
	_transTable[(size_t)'G'] = 3;
	_transTable[(size_t)'g'] = 3;
	_transTable[(size_t)'-'] = 4;
}

void ScoringMatrix::loadMatrix(std::string path) {
	std::string line;
	std::ifstream file(path);

	if (file.is_open()) {
		while (getline(file, line)){
			if (line == "") {
				continue;
			}
			std::istringstream iss(line);
			std::string token;
			std::vector<std::string> items;
			while (std::getline(iss, token, '\t')) {
				items.push_back(token);
			}

			if (items[0] == "mat") {
				int index = dnaToNum(items[1][0]); 
				double probability = std::stod(items[2]);
				_matrix[index][index] = std::log(probability);
			}
			else if (items[0] == "mis") {
				int x = dnaToNum(items[1][0]);
				int y = dnaToNum(items[1][items[1].size() - 1]);
				double probability = std::stod(items[2]);
				_matrix[x][y] = std::log(probability);
			}
			else if (items[0] == "del") {
				int x = dnaToNum(items[1][0]);
				int y = _yrange - 1;
				double probability = std::stod(items[2]);
				_matrix[x][y] = std::log(probability);
			}
			else if (items[0] == "ins") {
				int y = dnaToNum(items[1][0]);
				int x = _xrange - 1;
				double probability = std::stod(items[2]);
				_matrix[x][y] = std::log(probability);
			}	
		}
		file.close();
	}
}

double ScoringMatrix::getScore(char v, char w) 
{
	return _matrix[dnaToNum(v)][dnaToNum(w)];
}

void ScoringMatrix::printMatrix() {
	for (int i = 0; i < _xrange; i++) {
		for (int j = 0; j < _yrange; j++) {
			std::cout << std::fixed;
			std::cout << std::setw(4) << std::left 
					  << std::setprecision(2) << _matrix[i][j] << "\t";
		}
		std::cout << std::endl;
	}
}

ScoringMatrix::~ScoringMatrix() {
	for (int i = 0; i < _xrange; i++) {
		delete[] _matrix[i];
	}
	delete[] _transTable;
}
