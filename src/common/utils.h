//(c) 2016 by Authors
//This file is a part of ABruijn program.
//Released under the BSD license (see LICENSE file)

#pragma once

#include <vector>
#include <algorithm>
#include <sstream>

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
T quantile(const std::vector<T>& vec, int percent)
{
	//NOTE: there's a bug in libstdc++ nth_element, 
	//that sometimes leads to a segfault. This is why
	//we have this inefficient impleemntation here
	//std::nth_element(vec.begin(), vec.begin() + vec.size() / 2, 
	//				 vec.end());
	auto sortedVec = vec;
	std::sort(sortedVec.begin(), sortedVec.end());
	size_t targetId = std::min(vec.size() * (size_t)percent / 100, 
							   vec.size() - 1);
	return sortedVec[targetId];
}

template<typename T>
T median(const std::vector<T>& vec)
{
	return quantile(vec, 50);
}

inline std::vector<std::string> 
splitString(const std::string &s, char delim) 
{
	std::vector<std::string> elems;
	std::stringstream ss(s);
	std::string item;
	while (std::getline(ss, item, delim)) elems.push_back(item);
	return elems;
}
