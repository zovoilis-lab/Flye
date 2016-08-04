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
	float getScore(char v, char w) const;
private:
	void loadMatrix(const std::string& path);

	const int X_SIZE = 5;
	const int Y_SIZE = 5;
	std::vector<std::vector<float>> _matrix;
};

class HopoMatrix
{
public:
	//State represents a homopolymer that is observed in
	//reference sequence
	struct State
	{
		State():
			nucl(0), length(0), id(0)
		{}

		State(char nucl, uint32_t length);

		State(const std::string& str, size_t start = 0,
			  size_t end = std::string::npos);

		char nucl;
		uint32_t length;
		uint32_t id;
	};
	//Observation represents the read segment that corresponds
	//to a homopolymer in the reference (State). Might not be
	//a homopolymer, e.g. contain some other nucleotides, like
	//5A2X
	struct Observation
	{
		Observation(uint32_t id, bool extactMatch = false):
			id(id), extactMatch(extactMatch)
		{}
		uint32_t id;
		bool extactMatch;
	};
	typedef std::vector<Observation> ObsVector;

	HopoMatrix(const std::string& fileName);
	float getObsProb(State state, Observation observ) const
		{return _observationProbs[state.id][observ.id];}
	float getGenomeProb(State state) const
		{return _genomeProbs[state.id];}
	ObsVector knownObservations(State state) const;
	static Observation strToObs(char mainNucl, const std::string& dnaStr, 
								size_t start = 0, 
								size_t end = std::string::npos);

	//static std::string obsToStr(Observation obs);
private:
	void loadMatrix(const std::string& filaName);

	std::vector<std::vector<float>> _observationProbs;
	std::vector<float> 			 	_genomeProbs;
};
