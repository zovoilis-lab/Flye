//(c) 2016 by Authors
//This file is a part of ABruijn program.
//Released under the BSD license (see LICENSE file)

#include <chrono>
#include <thread>

#include "bubble_processor.h"

namespace
{
	size_t fileSize(const std::string& filename)
	{
		std::ifstream in(filename);
		in.ignore(std::numeric_limits<std::streamsize>::max());
	    return in.gcount();
	}
}

BubbleProcessor::BubbleProcessor(const std::string& subsMatPath,
								 const std::string& hopoMatrixPath):
	_subsMatrix(subsMatPath),
	_hopoMatrix(hopoMatrixPath),
	_generalPolisher(_subsMatrix),
	_homoPolisher(_subsMatrix, _hopoMatrix),
	_verbose(false)
{
}


void BubbleProcessor::polishAll(const std::string& inBubbles, 
								const std::string& outConsensus,
			   					int numThreads)
{
	_cachedBubbles.clear();
	_cachedBubbles.reserve(BUBBLES_CACHE);

	size_t fileLength = fileSize(inBubbles);
	_bubblesFile.open(inBubbles);
	if (!_bubblesFile.is_open())
	{
		throw std::runtime_error("Error opening bubbles file");
	}

	_progress.setFinalCount(fileLength);

	_consensusFile.open(outConsensus);
	if (!_consensusFile.is_open())
	{
		throw std::runtime_error("Error opening consensus file");
	}

	std::vector<std::thread> threads(numThreads);
	for (size_t i = 0; i < threads.size(); ++i)
	{
		threads[i] = std::thread(&BubbleProcessor::parallelWorker, this);
	}
	for (size_t i = 0; i < threads.size(); ++i)
	{
		threads[i].join();
	}
	_progress.setDone();
}


void BubbleProcessor::parallelWorker()
{
	_stateMutex.lock();
	while (true)
	{
		if (_cachedBubbles.empty())
		{
			if(_bubblesFile.eof())
			{
				_stateMutex.unlock();
				return;
			}
			else
			{
				this->cacheBubbles(BUBBLES_CACHE);
			}
		}

		Bubble bubble = _cachedBubbles.back();
		_cachedBubbles.pop_back();

		_stateMutex.unlock();
		_generalPolisher.polishBubble(bubble);
		_homoPolisher.polishBubble(bubble);
		_stateMutex.lock();
		
		this->writeBubbles({bubble});
		if (_verbose) this->writeLog({bubble});
	}
}


void BubbleProcessor::writeBubbles(const std::vector<Bubble>& bubbles)
{
	for (auto& bubble : bubbles)
	{
		_consensusFile << ">" << bubble.header << " " << bubble.position
			 		   << " " << bubble.branches.size() << std::endl
			 		   << bubble.candidate << std::endl;
	}
}

void BubbleProcessor::enableVerboseOutput(const std::string& filename)
{
	_verbose = true;
	_logFile.open(filename);
	if (!_logFile.is_open())
	{
		throw std::runtime_error("Error opening log file");
	}
}

void BubbleProcessor::writeLog(const std::vector<Bubble>& bubbles)
{
	std::vector<std::string> methods = {"None", "Insertion", "Substitution",
										"Deletion", "Homopolymer"};

	for (auto& bubble : bubbles)
	{
		for (auto& stepInfo : bubble.polishSteps)
		{
			 _logFile << std::fixed
				 << std::setw(22) << std::left << "Consensus: " 
				 << std::right << stepInfo.sequence << std::endl
				 << std::setw(22) << std::left << "Score: " << std::right 
				 << std::setprecision(2) << stepInfo.score << std::endl
				 << std::setw(22) << std::left << "Last method applied: " 
				 << std::right << methods[stepInfo.methodUsed] << std::endl;

			if (stepInfo.methodUsed == StepDel)
				_logFile << "Char at pos: " << stepInfo.changedIndex << " was deleted. \n";
			else if (stepInfo.methodUsed == StepSub)
				_logFile << "Char at pos " << stepInfo.changedIndex << " was substituted with " 
					<< "'" << stepInfo.changedLetter << "'.\n";
			else if (stepInfo.methodUsed == StepIns)
				_logFile << "'"<< stepInfo.changedLetter << "'" 
					 << " was inserted at pos " << stepInfo.changedIndex << ".\n";

			_logFile << std::endl;
		}
		_logFile << "-----------------\n";
	}
}


void BubbleProcessor::cacheBubbles(int maxRead)
{
	std::string buffer;
	std::string candidate;

	int readBubbles = 0;
	while (!_bubblesFile.eof() && readBubbles < maxRead)
	{
		std::getline(_bubblesFile, buffer);
		if (buffer.empty()) break;

		std::vector<std::string> elems = splitString(buffer, ' ');
		if (elems.size() < 3 || elems[0][0] != '>')
		{
			throw std::runtime_error("Error parsing bubbles file");
		}
		std::getline(_bubblesFile, candidate);
		std::transform(candidate.begin(), candidate.end(), 
				       candidate.begin(), ::toupper);
		
		Bubble bubble;
		bubble.candidate = candidate;
		bubble.header = elems[0].substr(1, std::string::npos);
		bubble.position = std::stoi(elems[1]);
		int numOfReads = std::stoi(elems[2]);

		int count = 0;
		while (count < numOfReads) 
		{
			if (buffer.empty()) break;

			std::getline(_bubblesFile, buffer);
			std::getline(_bubblesFile, buffer);
			std::transform(buffer.begin(), buffer.end(), 
				       	   buffer.begin(), ::toupper);
			bubble.branches.push_back(buffer);
			count++;
		}
		if (count != numOfReads)
		{
			throw std::runtime_error("Error parsing bubbles file");
		}

		_cachedBubbles.push_back(std::move(bubble));
		++readBubbles;
	}

	int filePos = _bubblesFile.tellg();
	if (filePos > 0)
	{
		_progress.setValue(filePos);
	}
}
