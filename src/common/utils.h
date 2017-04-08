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
