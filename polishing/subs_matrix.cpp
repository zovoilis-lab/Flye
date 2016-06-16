//(c) 2016 by Authors
//This file is a part of ABruijn program.
//Released under the BSD license (see LICENSE file)

#include <sstream>
#include <vector>
#include <cmath>
#include <stdexcept>
#include <cassert>

#include "subs_matrix.h"
#include "utility.h"

namespace
{
	static std::vector<size_t> table;
	size_t dnaToId(char c)
	{
		return table[(size_t)c];
	}
	struct TableFiller
	{
		TableFiller()
		{
			static bool tableFilled = false;
			if (!tableFilled)
			{
				tableFilled = true;
				table.assign(256, -1);	//256 chars
				table[(size_t)'A'] = 0;
				table[(size_t)'a'] = 0;
				table[(size_t)'C'] = 1;
				table[(size_t)'c'] = 1;
				table[(size_t)'T'] = 2;
				table[(size_t)'t'] = 2;
				table[(size_t)'G'] = 3;
				table[(size_t)'g'] = 3;
				table[(size_t)'-'] = 4;
			}
		}
	};
	TableFiller filler;

	static const size_t NUM_STATES = 128;
	static const size_t NUM_OBS = 65536;
	static const size_t MAX_RUN = 15;
	static const double ZERO_PROB = 0.001f;
}

SubstitutionMatrix::SubstitutionMatrix(const std::string& path)
{	
	for (int i = 0; i < X_SIZE; i++) 
	{
		_matrix.push_back(std::vector<double>(Y_SIZE, 0));
	}
	this->loadMatrix(path);
}

void SubstitutionMatrix::loadMatrix(const std::string& path) 
{
	std::string line;
	std::ifstream file(path);

	if (!file.is_open()) 
	{
		throw std::runtime_error("Can't open substitution matrix");
	}
	while (std::getline(file, line))
	{
		if (line.empty()) continue;
		if (line[line.length() - 1] == '\r') line.pop_back();

		auto items = splitString(line, '\t');

		if (items[0] == "mat") 
		{
			int index = dnaToId(items[1][0]); 
			double probability = std::stod(items[2]);
			_matrix[index][index] = std::log(probability);
		}
		else if (items[0] == "mis") 
		{
			int x = dnaToId(items[1][0]);
			int y = dnaToId(items[1][items[1].size() - 1]);
			double probability = std::stod(items[2]);
			_matrix[x][y] = std::log(probability);
		}
		else if (items[0] == "del") 
		{
			int x = dnaToId(items[1][0]);
			int y = Y_SIZE - 1;
			double probability = std::stod(items[2]);
			_matrix[x][y] = std::log(probability);
		}
		else if (items[0] == "ins") 
		{
			int y = dnaToId(items[1][0]);
			int x = X_SIZE - 1;
			double probability = std::stod(items[2]);
			_matrix[x][y] = std::log(probability);
		}	
	}
}

double SubstitutionMatrix::getScore(char v, char w) const 
{
	return _matrix[dnaToId(v)][dnaToId(w)];
}

void SubstitutionMatrix::printMatrix() const 
{
	for (int i = 0; i < X_SIZE; i++) {
		for (int j = 0; j < Y_SIZE; j++) {
			std::cout << std::fixed;
			std::cout << std::setw(4) << std::left 
					  << std::setprecision(2) << _matrix[i][j] << "\t";
		}
		std::cout << std::endl;
	}
}

namespace
{
	std::string expandHopo(const std::string& hopo)
	{
		std::string buf;
		std::string result;
		for (char c : hopo)
		{
			if (!std::isalpha(c))
			{
				buf += c;
			}
			else
			{
				int runLen = std::stoi(buf);
				result += std::string(runLen, c);
				buf.clear();
			}
		}
		//std::cerr << hopo << " " << result << std::endl;
		return result;
	}
}


HopoMatrix::State::State(char nucl, char length):
	nucl(nucl), length(length)
{
	 hash = std::min(length, (char)MAX_RUN) + 
					(MAX_RUN + 1) * dnaToId(nucl);

}
		
HopoMatrix::State::State(const std::string& str, size_t start, size_t end)
{
	if (end == std::string::npos) end = str.length();
	assert(end - start > 0);

	//std::cerr << str.substr(start, end - start) << std::endl;
	size_t runLength = 0;
	char runNucl = -1;
	for (size_t i = start; i < end; ++i)
	{
		if (str[i] != '-')
		{
			if (runNucl == -1) runNucl = str[i];
			if (str[i] != runNucl) throw std::runtime_error("Wrong homopolymer");
			++runLength;
		}
	}
	if (runLength == 0) throw std::runtime_error("Wrong homopolymer");

	hash = std::min(runLength, MAX_RUN) + (MAX_RUN + 1) * dnaToId(runNucl);
	nucl = runNucl;
	length = runLength;
}


HopoMatrix::Observation HopoMatrix::strToObs(const std::string& str, 
											 size_t begin, size_t end)
{
	std::vector<size_t> counts = {0, 0, 0, 0, 0};
	uint16_t result = 0;
	if (end == std::string::npos) end = str.length();
	for (size_t pos = begin; pos < end; ++pos) ++counts[dnaToId(str[pos])];

	for (size_t i = 0; i < 4; ++i)
	{
		counts[i] = std::min(counts[i], MAX_RUN);
		result <<= 4;
		result += (uint16_t)counts[i];
	}
	//std::cerr << str << " " << result << std::endl;
	return result;
}

std::string HopoMatrix::obsToStr(HopoMatrix::Observation obs)
{
	//ACTG
	char NUCLS[] = "GTCA";
	std::string result;
	for (size_t i = 0; i < 4; ++i)
	{
		int num = obs & MAX_RUN;
		obs >>= 4;
		if (num) result += std::to_string(num) + NUCLS[i];
	}
	return result;
}

HopoMatrix::HopoMatrix(const std::string& fileName)
{
	for (size_t i = 0; i < NUM_STATES; ++i)
	{
		_matrix.push_back(std::vector<double>(NUM_OBS, 0.0f));
	}
	this->loadMatrix(fileName);
}


void HopoMatrix::loadMatrix(const std::string& fileName)
{
	std::ifstream fin(fileName);
	std::string buffer;
	if (!fin.is_open()) 
	{
		throw std::runtime_error("Can't open homopolymer matrix");
	}

	std::vector<std::vector<int>> frequencies;
	for (size_t i = 0; i < NUM_STATES; ++i)
	{
		frequencies.push_back(std::vector<int>(NUM_OBS, 0));
	}

	while (std::getline(fin, buffer))
	{
		if (buffer.empty()) continue;
		if (buffer[buffer.length() - 1] == '\r') buffer.pop_back();

		auto tokens = splitString(buffer, '\t');
		tokens[0].pop_back();
		State state = State(expandHopo(tokens[0]));
		for (size_t i = 1; i < tokens.size(); ++i)
		{
			auto obsTokens = splitString(tokens[i], '=');
			Observation obs = strToObs(expandHopo(obsTokens[0]));
			frequencies[state.hash][obs] = std::stoi(obsTokens[1]);
		}
	}

	for (size_t i = 0; i < NUM_STATES; ++i)
	{
		double sumFreq = 0;
		for (size_t j = 0; j < NUM_OBS; ++j)
		{
			sumFreq += frequencies[i][j];
		}
		for (size_t j = 0; j < NUM_OBS; ++j)
		{
			double prob = (sumFreq > 0.0f) ? 
						  (double)frequencies[i][j] / sumFreq : 0.0f;
			_matrix[i][j] = std::log(std::max(prob, ZERO_PROB));
		}
	}
}
