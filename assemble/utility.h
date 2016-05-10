//(c) 2013-2016 by Authors
//This file is a part of Ragout program.
//Released under the BSD license (see LICENSE file)

#pragma once

#include <iterator>
#include <utility>
#include <iostream>
#include <ctime>

#ifdef _DEBUG_LOG
#define TIMESTAMP << timestamp()
#else
#define TIMESTAMP
#endif

#ifdef _DEBUG_LOG
#define DEBUG_PRINT(x) do {std::cerr TIMESTAMP << " " << x << std::endl;} \
					   while(0)
#else
#define DEBUG_PRINT(x)
#endif

#ifdef _LOG
#define LOG_PRINT(x) do {std::cerr TIMESTAMP << " " << x << std::endl;} \
					   while(0)
#else
#define LOG_PRINT(x)
#endif

#define WARNING_PRINT(x) do {std::cerr TIMESTAMP << " [WARNING] " << x << std::endl;} \
					   	 while(0)

#define ERROR_PRINT(x) do {std::cerr TIMESTAMP << " [ERROR] " << x << std::endl;} \
					   	 while(0)

inline std::string timestamp(const char* format = "[%H:%M:%S]")
{
	std::time_t t = std::time(0);
	char cstr[128];
	std::strftime(cstr, sizeof(cstr), format, std::localtime(&t));
	return cstr;
}
