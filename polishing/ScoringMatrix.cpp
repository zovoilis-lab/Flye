#include "ScoringMatrix.h"
#include <sstream>
#include <vector>
#include <cmath>
#include <stdexcept>

namespace
{
	int dnaToNum(char c)
	{
		c = toupper(c);
		switch (c)
		{
		case 'A':
			return 0;
		case 'C':
			return 1;
		case 'T':
			return 2;
		case 'G':
			return 3;
		case '-':
			return 4;
		}
		throw std::runtime_error("Unknown DNA base");
		return -1;
	}
}

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
				int index = dnaToNum(items[1][0]); 
				double probability = std::stod(items[2]);
				m_matrix[index][index] = std::log(probability);
			}
			else if (items[0] == "mis") {
				int x = dnaToNum(items[1][0]);
				int y = dnaToNum(items[1][items[1].size() - 1]);
				double probability = std::stod(items[2]);
				m_matrix[x][y] = std::log(probability);
			}
			else if (items[0] == "del") {
				int x = dnaToNum(items[1][0]);
				int y = m_yrange - 1;
				double probability = std::stod(items[2]);
				m_matrix[x][y] = std::log(probability);
			}
			else if (items[0] == "ins") {
				int y = dnaToNum(items[1][0]);
				int x = m_xrange - 1;
				double probability = std::stod(items[2]);
				m_matrix[x][y] = std::log(probability);
			}	
		}
		file.close();
	}
}

double ScoringMatrix::getScore(char v, char w) 
{
	return m_matrix[dnaToNum(v)][dnaToNum(w)];
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
