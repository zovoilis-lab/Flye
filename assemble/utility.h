//(c) 2013-2016 by Authors
//This file is a part of Ragout program.
//Released under the BSD license (see LICENSE file)

#pragma once

#include <iterator>
#include <utility>
#include <iostream>
#include <ctime>
#include <algorithm>

#ifdef _WIN32
#include <direct.h>
#else
#include <sys/stat.h>
#endif

#ifdef _DEBUG_LOG
#define DEBUG_PRINT(x) do {std::cerr << timestamp() << " " << x << std::endl;} \
					   while(0)
#else
#define DEBUG_PRINT(x)
#endif

#ifdef _LOG
#define LOG_PRINT(x) do {std::cerr << timestamp() << " " << x << std::endl;} \
					   while(0)
#else
#define LOG_PRINT(x)
#endif

#define WARNING_PRINT(x) do {std::cerr << timestamp() << " [WARNING] " << x << std::endl;} \
					   	 while(0)

inline std::string timestamp(const char* format = "[%H:%M:%S]")
{
	std::time_t t = std::time(0);
	char cstr[128];
	std::strftime(cstr, sizeof(cstr), format, std::localtime(&t));
	return cstr;
}


inline bool makeDirectory(const std::string& name)
{
#ifdef _WIN32
	int ret = _mkdir(name.c_str());
#else
	int ret = mkdir(name.c_str(), 0777);
#endif
	return !ret;
}
