//(c) 2016 by Authors
//This file is a part of ABruijn program.
//Released under the BSD license (see LICENSE file)

#include <chrono>
#include <thread>

#include "bubble_processor.h"


BubbleProcessor::BubbleProcessor(const std::string& subsMatPath,
								 const std::string& hopoMatrixPath):
	_subsMatrix(subsMatPath),
	_hopoMatrix(hopoMatrixPath),
	_generalPolisher(_subsMatrix),
	_homoPolisher(_subsMatrix, _hopoMatrix),
	_nextJob(0)
{
}


void BubbleProcessor::polishAll(const std::string& dataPath, int numThreads) 
{
	this->readBubbles(dataPath);
	_progress.setFinalCount(_bubbles.size());
	_nextJob = 0;

	std::vector<std::thread> threads(numThreads);
	for (size_t i = 0; i < threads.size(); ++i)
	{
		threads[i] = std::thread(&BubbleProcessor::parallelWorker, this);
	}
	for (size_t i = 0; i < threads.size(); ++i)
	{
		threads[i].join();
	}
}


void BubbleProcessor::parallelWorker()
{
	_stateMutex.lock();
	while (true)
	{
		if (_nextJob == _bubbles.size())
		{
			_stateMutex.unlock();
			return;
		}

		_progress.advance();
		size_t curJob = _nextJob++;
		Bubble bubble = _bubbles[curJob];

		_stateMutex.unlock();
		_generalPolisher.polishBubble(bubble);
		_homoPolisher.polishBubble(bubble);
		_stateMutex.lock();
		
		_bubbles[curJob] = std::move(bubble);
	}
}


void BubbleProcessor::writeConsensuses(const std::string& fileName)
{
	std::ofstream fout(fileName);
	for (auto& bubble : _bubbles)
	{
		fout << ">" << bubble.header << " " << bubble.position
			 << " " << bubble.branches.size() << std::endl
			 << bubble.candidate << std::endl;
	}
}

void BubbleProcessor::writeLog(const std::string& fileName)
{
	std::ofstream fout(fileName);
	std::vector<std::string> methods = {"None", "Insertion", "Substitution",
										"Deletion", "Homopolymer"};

	for (auto& bubble : _bubbles)
	{
		for (auto& stepInfo : bubble.polishSteps)
		{
			fout << std::fixed
				 << std::setw(22) << std::left << "Consensus: " 
				 << std::right << stepInfo.sequence << std::endl
				 << std::setw(22) << std::left << "Score: " << std::right 
				 << std::setprecision(2) << stepInfo.score << std::endl
				 << std::setw(22) << std::left << "Last method applied: " 
				 << std::right << methods[stepInfo.methodUsed] << std::endl;

			if (stepInfo.methodUsed == StepDel)
				fout << "Char at pos: " << stepInfo.changedIndex << " was deleted. \n";
			else if (stepInfo.methodUsed == StepSub)
				fout << "Char at pos " << stepInfo.changedIndex << " was substituted with " 
					<< "'" << stepInfo.changedLetter << "'.\n";
			else if (stepInfo.methodUsed == StepIns)
				fout << "'"<< stepInfo.changedLetter << "'" 
					 << " was inserted at pos " << stepInfo.changedIndex << ".\n";

			fout << std::endl;
		}
		fout << "-----------------\n";
	}
}


void BubbleProcessor::readBubbles(const std::string& fileName)
{
	std::string buffer;
	std::ifstream file(fileName);
	std::string candidate;
	_bubbles.clear();

	if (!file.is_open())
	{
		throw std::runtime_error("Error opening bubble file");
	}

	while (!file.eof())
	{
		std::getline(file, buffer);
		if (buffer.empty()) break;

		std::vector<std::string> elems = splitString(buffer, ' ');
		if (elems.size() < 3 || elems[0][0] != '>')
		{
			throw std::runtime_error("Error parsing bubbles file");
		}
		std::getline(file, candidate);
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

			std::getline(file, buffer);
			std::getline(file, buffer);
			std::transform(buffer.begin(), buffer.end(), 
				       	   buffer.begin(), ::toupper);
			bubble.branches.push_back(buffer);
			count++;
		}
		if (count != numOfReads)
		{
			throw std::runtime_error("Error parsing bubbles file");
		}

		_bubbles.push_back(std::move(bubble));
	}
}
