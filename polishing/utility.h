//(c) 2013-2016 by Authors
//This file is a part of Ragout program.
//Released under the BSD license (see LICENSE file)

#pragma once

#include <iterator>
#include <utility>
#include <iostream>
#include <ctime>
#include <sstream>
#include <vector>


inline std::vector<std::string> 
splitString(const std::string &s, char delim) 
{
	std::vector<std::string> elems;
	std::stringstream ss(s);
	std::string item;
	while (std::getline(ss, item, delim)) elems.push_back(item);
	return elems;
}

class ProgressPercent
{
public:
	ProgressPercent(size_t finalCount = 0):
		_finalCount(finalCount), _curCount(0), _prevPercent(-1),
		_stopped(false)
	{}

	void setFinalCount(size_t finalCount) {_finalCount = finalCount;}
	void setValue(size_t value)
	{
		this->advance(value - _curCount);
	}
	void setDone()
	{
		this->setValue(_finalCount);
	}
	void advance(size_t step = 1)
	{
		if (_stopped) return;

		_curCount += step;
		int percent = 10UL * _curCount / _finalCount;
		if (percent > _prevPercent)
		{
			std::cerr << percent * 10 << "% ";
			_prevPercent = percent;
		}

		if (_prevPercent >= 10)
		{
			std::cerr << std::endl;
			_stopped = true;
		}
	}

private:
	size_t _finalCount;
	size_t _curCount;
	int  _prevPercent;
	bool _stopped;
};
