#include "ScoringMatrix.h"
#include <sstream>
#include <vector>
#include <cmath>


ScoringMatrix::ScoringMatrix(int xrange, int yrange) {	

	m_matrix = new double*[xrange];
	for (int i = 0; i < xrange; i++) {
		m_matrix[i] = new double[yrange];
	}
	m_xrange = xrange;
	m_yrange = yrange;

	m_matrix[m_xrange - 1][m_yrange - 1] = 0;
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
				int index = getIndex(items[1][0]); 
				double probability = std::stod(items[2]);
				m_matrix[index][index] = std::log(probability);
			}
			else if (items[0] == "mis") {
				int x = getIndex(items[1][0]);
				int y = getIndex(items[1][items[1].size() - 1]);
				double probability = std::stod(items[2]);
				m_matrix[x][y] = std::log(probability);
			}
			else if (items[0] == "del") {
				int x = getIndex(items[1][0]);
				int y = m_yrange - 1;
				double probability = std::stod(items[2]);
				m_matrix[x][y] = std::log(probability);
			}
			else if (items[0] == "ins") {
				int y = getIndex(items[1][0]);
				int x = m_xrange - 1;
				double probability = std::stod(items[2]);
				m_matrix[x][y] = std::log(probability);
			}	
		}
		file.close();
	}
}

double ScoringMatrix::getScore(char v, char w) {
	int x = getIndex(v);
	int y = getIndex(w);
	return m_matrix[x][y];
}

int ScoringMatrix::getIndex(char c) {
	c = toupper(c);
	int value = -1;
	if (c == 'A') {
		value = 0;
	}
	else if (c == 'C') {
		value = 1;
	}
	else if (c == 'T') {
		value = 2;
	}
	else if (c == 'G') {
		value = 3;
	}
	else if (c == '-') {
		value = 4;
	}

	if (value == -1) {
		std::cout << "Error in get Index!\n";
	}

	return value;
}

void ScoringMatrix::printMatrix() {
	for (int i = 0; i < m_xrange; i++) {
		for (int j = 0; j < m_yrange; j++) {
			std::cout << std::fixed;
			std::cout << std::setw(4) << std::left << std::setprecision(2) << m_matrix[i][j] << "\t";
		}
		std::cout << std::endl;
	}
}

ScoringMatrix::~ScoringMatrix() {
	for (int i = 0; i < m_xrange; i++) {
		delete[] m_matrix[i];
	}
}
