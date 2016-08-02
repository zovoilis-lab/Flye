//(c) 2013-2016 by Authors
//This file is a part of Ragout program.
//Released under the BSD license (see LICENSE file)

#pragma once

#include <iostream>
#include <fstream>
#include <ctime>
#include <stdexcept>
#include <atomic>

class Logger
{
public:
	static Logger& get()
	{
		static Logger instance;
		return instance;
	}

	void setOutputFile(const std::string& filename)
	{
		_logFile.open(filename, std::ofstream::out | std::ofstream::app);
		_logFileSet = true;
		if (!_logFile.is_open()) 
		{
			throw std::runtime_error("Can't open log file");
		}
		_logFile << "-----------Begin assembly log------------\n";
	}

	void setDebugging(bool debug) {_debug = debug;}

	class StreamWriter
	{
	public:
		StreamWriter(const std::string& level, 
					 std::ostream* consoleStream = nullptr,
					 std::ostream* fileStream = nullptr):
			_fileStream(fileStream), _consoleStream(consoleStream)
		{
			if (_fileStream) *_fileStream << timestamp() << " " << level << " ";
			if (_consoleStream) *_consoleStream << timestamp() 
												<< " " << level << " ";
		}
		~StreamWriter()
		{
			if (_fileStream) *_fileStream << std::endl;
			if (_consoleStream) *_consoleStream << std::endl;
		}
		
		template <class T>
		Logger::StreamWriter& operator<< (const T& val)
		{
			if (_fileStream) *_fileStream << val;
			if (_consoleStream) *_consoleStream << val;
			return *this;
		}

	private:
		std::ostream* _fileStream;
		std::ostream* _consoleStream;
	};

	StreamWriter info()
	{
		std::ostream* logPtr = _logFileSet ? &_logFile : nullptr;
		return StreamWriter("INFO:", &std::cerr, logPtr);
	}

	StreamWriter warning()
	{
		std::ostream* logPtr = _logFileSet ? &_logFile : nullptr;
		return StreamWriter("WARNING:", &std::cerr, logPtr);
	}

	StreamWriter error()
	{
		std::ostream* logPtr = _logFileSet ? &_logFile : nullptr;
		return StreamWriter("ERROR:", &std::cerr, logPtr);
	}

	StreamWriter debug()
	{
		std::ostream* logPtr = _logFileSet ? &_logFile : nullptr;
		std::ostream* consolePtr = _debug ? &std::cerr : nullptr;
		return StreamWriter("DEBUG:", consolePtr, logPtr);
	}

private:
	static std::string timestamp(const char* format = "[%H:%M:%S]")
	{
		std::time_t t = std::time(0);
		char cstr[128];
		std::strftime(cstr, sizeof(cstr), format, std::localtime(&t));
		return cstr;
	}

	Logger():
		_debug(false), _logFileSet(false)
	{}
	~Logger()
	{
		if (_logFileSet)
		{
			_logFile << "-----------End assembly log------------\n";
		}
	}

	bool _debug;
	bool _logFileSet;
	std::ofstream _logFile;
};

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
			int expected = percent - 1;
			if (_prevPercent.compare_exchange_weak(expected, percent))
			{
				std::cerr << percent * 10 << "% ";
				if (percent >= 10)
				{
					std::cerr << std::endl;
					_stopped = true;
				}
			}
		}
	}

private:
	size_t 			    _finalCount;
	std::atomic<size_t> _curCount;
	std::atomic<int>  	_prevPercent;
	bool 				_stopped;
};
