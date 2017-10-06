//(c) 2016 by Authors
//This file is a part of ABruijn program.
//Released under the BSD license (see LICENSE file)

#pragma once

#include <vector>
#include <algorithm>

template<class T>
void vecRemove(std::vector<T>& v, T val)
{
	v.erase(std::remove(v.begin(), v.end(), val), v.end()); 
}

struct pairhash 
{
public:
	template <typename T, typename U>
	std::size_t operator()(const std::pair<T, U> &x) const
	{
		return std::hash<T>()(x.first) ^ std::hash<U>()(x.second);
	}
};

template<typename T>
T median(std::vector<T>& vec)
{
	std::sort(vec.begin(), vec.end());
	//NOTE: there's a bug in libstdc++ nth_element, 
	//that sometimes leads to a segfault
	//std::nth_element(vec.begin(), vec.begin() + vec.size() / 2, 
	//				 vec.end());
	if (vec.size() % 2 == 1)
	{
		return vec[vec.size() / 2];
	}
	else
	{
		return (vec[vec.size() / 2 - 1] + vec[vec.size() / 2]) / 2;
	}
}

template<typename T>
T q75(std::vector<T>& vec)
{
	std::sort(vec.begin(), vec.end());
	//NOTE: there's a bug in libstdc++ nth_element, 
	//that sometimes leads to a segfault
	//std::nth_element(vec.begin(), vec.begin() + vec.size() / 2, 
	//				 vec.end());
	return vec[vec.size() * 3 / 4];
}
