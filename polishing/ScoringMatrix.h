#pragma once

#include <string>
#include <fstream>
#include <iostream>
#include <iomanip>

class ScoringMatrix 
{
public:
	ScoringMatrix(int xrange, int yrange);
	void loadMatrix(std::string path);
	double getScore(char v, char w);
	void printMatrix();
	~ScoringMatrix();

private:
	int m_xrange;
	int m_yrange;
	double** m_matrix;
};

