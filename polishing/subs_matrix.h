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
			nucl(0), length(0), id(0)
		{}
		State(char nucl, uint32_t length);
		//State(uint32_t id);
		State(const std::string& str, size_t start = 0,
			  size_t end = std::string::npos);

		char nucl;
		uint32_t length;
		uint32_t id;
	};
	struct Observation
	{
		uint32_t id;
		bool extactMatch;
	};

	HopoMatrix(const std::string& fileName);
	double getObsProb(State state, Observation observ) const
		{return _observationProbs[state.id][observ.id];}
	double getGenomeProb(State state) const
		{return _genomeProbs[state.id];}

	static Observation strToObs(char mainNucl, const std::string& dnaStr, 
								size_t start = 0, 
								size_t end = std::string::npos);
	//static std::string obsToStr(Observation obs);
private:
	void loadMatrix(const std::string& filaName);

	std::vector<std::vector<double>> _observationProbs;
	std::vector<double> 			 _genomeProbs;
};
