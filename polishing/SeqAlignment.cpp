// SeqAlignment.cpp : Defines the entry point for the console application.
//

#include "Worker.h"
#include <ctime>
#include <iostream>

#include <stdio.h>
#include <time.h>
#include <omp.h>


int main(int argc, char* argv[]) {
	const int shortParamCount = 4 * 2 + 1;
	const int longParamCount = 6 * 2 + 1;

	std::string scoringMatPath = "";
	std::string reads = "";
	std::string type = "";
	std::string format = "";
	size_t start = std::numeric_limits<size_t>::max();
	size_t end = std::numeric_limits<size_t>::max();
	

	if (argc != shortParamCount && argc != longParamCount) {
		std::cerr << "The number of parameters is incorrect." << std::endl;
		return 0;
	}

	for (int i = 1; i < argc; ++i) {
		//Scoring matrix
		if (std::string(argv[i]) == "--scoringMatrix" || std::string(argv[i]) == "-sm") {
			if (i + 1 < argc) { // Make sure we aren't at the end of argv!
				scoringMatPath = std::string(argv[++i]);
			}
			else { 
				std::cerr << "Wrong number of parameters!" << std::endl;
			}
			continue;
		}

		//Reads
		if (std::string(argv[i]) == "--reads" || std::string(argv[i]) == "-r") {
			if (i + 1 < argc) { // Make sure we aren't at the end of argv!
				reads = std::string(argv[++i]);
			}
			else {
				std::cerr << "Wrong number of parameters!" << std::endl;
			}
			continue;
		}

		//Type
		if (std::string(argv[i]) == "--type" || std::string(argv[i]) == "-t") {
			if (i + 1 < argc) { // Make sure we aren't at the end of argv!
				type = std::string(argv[++i]);
			}
			else {
				std::cerr << "Wrong number of parameters!" << std::endl;
			}

			if (type != "single" && type != "many") {
				std::cerr << "The parameter passed to the 'type' is incorrect! "
							 "It should be 'single' or 'many'." << std::endl;
				return 0;
			}
			continue;
		}

		//Format
		if (std::string(argv[i]) == "--outputFormat" || std::string(argv[i]) == "-of") {
			if (i + 1 < argc) { // Make sure we aren't at the end of argv!
				format = std::string(argv[++i]);
			}
			else {
				std::cerr << "Wrong number of parameters!" << std::endl;
			}

			if (format != "short" && format != "verbose") 
			{
				std::cerr << "The parameter passed to the 'outputFormat' is incorrect! "
							 "It should be 'short' or 'verbose'." << std::endl;
				return 0;
			}

			continue;
		}

		//Optional start/ends params
		//Start
		if (std::string(argv[i]) == "--start" || std::string(argv[i]) == "-s") {
			int temp;
			if (i + 1 < argc) { // Make sure we aren't at the end of argv!
				temp = stoi(std::string(argv[++i]));
				if (temp < 0) {
					std::cerr << "The parameter passed to the 'start' is less than zero!" << std::endl;
					return 0;
				}
				start = temp;
			}
			else {
				std::cerr << "Wrong number of parameters!" << std::endl;
			}
			continue;
		}

		//End
		if (std::string(argv[i]) == "--end" || std::string(argv[i]) == "-e") {
			int temp;
			if (i + 1 < argc) { // Make sure we aren't at the end of argv!
				temp = stoi(std::string(argv[++i]));
				if (temp < 1) {
					std::cerr << "The parameter passed to the 'end' is less than one!" << std::endl;
					return 0;
				}
				end = temp;
			}
			else {
				std::cerr << "Wrong number of parameters!" << std::endl;
			}
			continue;
		}
	}

	std::clock_t startTimer;
	startTimer = std::clock();

	if (scoringMatPath.compare("") == 0) {
		std::cerr << "The scoring matrix was not given." << std::endl;
		return 0;
	}

	Worker worker(scoringMatPath); //"../input/alignment_probs_newchem_20k_noins.txt"

	if (reads.compare("") == 0) {
		std::cerr << "The file with reads was not given." << std::endl;
		return 0;
	}

	if (type.compare("single") == 0) {
		std::cout << "Running " << argv[0] << "..... \n";
		worker.run(reads, format); //"../input/nhood.2.fasta"
	}

	else if (type.compare("many") == 0){
		std::cout << "Running " << argv[0] << "..... \n";
		if (start == std::numeric_limits<size_t>::max() || 
			end == std::numeric_limits<size_t>::max()) 
		{
			std::cerr << "Start or end or both were not set!" << std::endl;
			return 0;
		}
		worker.run(start, end, reads, format); //"../input/reads.20k.fasta.15.8.24.v0.fasta.all.bubble1.fasta"
	}
	else {
		std::cerr << "Incorrect type. Should be 'single'/'many'." << std::endl;
		return 0;
	}


	std::cout << "Run time: " 
			  << (std::clock() - startTimer) / (double)(CLOCKS_PER_SEC / 1000) 
			  << " ms" << std::endl;
	cin.get();
	return 0;
}


