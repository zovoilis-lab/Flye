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
	int dnaToNum(char dna)
		{return _transTable[(size_t)dna];}

	int _xrange;
	int _yrange;
	int* _transTable;
	double** _matrix;
};

