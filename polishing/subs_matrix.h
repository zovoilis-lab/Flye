//(c) 2016 by Authors
//This file is a part of ABruijn program.
//Released under the BSD license (see LICENSE file)

#pragma once

#include <string>
#include <fstream>
#include <iostream>
#include <iomanip>
#include <vector>

class SubstitutionMatrix 
{
public:
	SubstitutionMatrix();
	void loadMatrix(std::string path);
	double getScore(char v, char w) const;
	void printMatrix() const;

private:
	static int dnaToNum(char dna)
		{return _transTable[(size_t)dna];}
	static std::vector<int> _transTable;	

	const int X_SIZE = 5;
	const int Y_SIZE = 5;
	std::vector<std::vector<double>> _matrix;
};
