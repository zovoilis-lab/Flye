// SeqAlignment.cpp : Defines the entry point for the console application.
//

#include <ctime>
#include <iostream>

#include "Worker.h"

int main(int argc, char* argv[]) 
{
	const int longParamCount = 3 * 2 + 1;

	std::string scoringMatPath = "";
	std::string reads = "";
	std::string type = "";
	std::string format = "";

	if (argc != longParamCount) {
		std::cerr << "The number of parameters is incorrect." << std::endl;
		return 1;
	}

	for (int i = 1; i < argc; ++i) 
	{
		//Scoring matrix
		if (std::string(argv[i]) == "--scoringMatrix" 
			|| std::string(argv[i]) == "-sm") 
		{
			if (i + 1 < argc)
				scoringMatPath = std::string(argv[++i]);
			else
				std::cerr << "Wrong number of parameters!" << std::endl;
			continue;
		}

		//Reads
		if (std::string(argv[i]) == "--reads" 
			|| std::string(argv[i]) == "-r") 
		{
			if (i + 1 < argc)
				reads = std::string(argv[++i]);
			else
				std::cerr << "Wrong number of parameters!" << std::endl;
			continue;
		}

		//Format
		if (std::string(argv[i]) == "--outputFormat" 
			|| std::string(argv[i]) == "-of") 
		{
			if (i + 1 < argc)
				format = std::string(argv[++i]);
			else
				std::cerr << "Wrong number of parameters!" << std::endl;

			if (format != "short" && format != "verbose") 
			{
				std::cerr << "The parameter passed to the 'outputFormat' is incorrect! "
							 "It should be 'short' or 'verbose'." << std::endl;
				return 1;
			}
			continue;
		}
	}

	std::clock_t startTimer;
	startTimer = std::clock();

	if (scoringMatPath.empty()) {
		std::cerr << "The scoring matrix was not given." << std::endl;
		return 1;
	}
	Worker worker(scoringMatPath);

	if (reads.empty()) {
		std::cerr << "The file with reads was not given." << std::endl;
		return 1;
	}

	worker.run(reads, format); 

	std::cout << "Run time: " 
			  << (std::clock() - startTimer) / (double)(CLOCKS_PER_SEC / 1000) 
			  << " ms" << std::endl;
	return 0;
}


