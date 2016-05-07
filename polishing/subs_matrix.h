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
	SubstitutionMatrix(const std::string& path);
	double getScore(char v, char w) const;
	void printMatrix() const;
private:
	void loadMatrix(const std::string& path);

	const int X_SIZE = 5;
	const int Y_SIZE = 5;
	std::vector<std::vector<double>> _matrix;
};

class HopoMatrix
{
public:
	struct State
	{
		State():
			nucl(0), length(0), hash(0)
		{}
		State(char nucl, char length);
		State(const std::string& str, size_t start = 0,
			  size_t end = std::string::npos);

		char nucl;
		char length;
		uint16_t hash;
	};
	typedef uint16_t Observation;

	HopoMatrix(const std::string& fileName);
	double getScore(State state, Observation observ) const
		{return _matrix[state.hash][observ];}

	static Observation strToObs(const std::string& str, size_t start = 0,
					   		    size_t end = std::string::npos);
	static std::string obsToStr(Observation obs);
private:
	void loadMatrix(const std::string& filaName);
	std::vector<std::vector<double>> _matrix;
};
