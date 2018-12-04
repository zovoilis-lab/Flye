//(c) 2016 by Authors
//This file is a part of ABruijn program.
//Released under the BSD license (see LICENSE file)

#include <ctime>
#include <iostream>
#include <getopt.h>

#include "../polishing/bubble_processor.h"


bool parseArgs(int argc, char** argv, std::string& bubblesFile, 
			   std::string& scoringMatrix, std::string& hopoMatrix,
			   std::string& outConsensus, std::string& outVerbose,
			   int& numThreads, bool& quiet)
{
	auto printUsage = []()
	{
		std::cerr << "Usage: polish bubbles_file subs_matrix "
				  << "hopo_matrix out_file \n\t\t"
				  << "[-v verbose_log] [-t num_threads] [-q]\n\n"
				  << "positional arguments:\n"
				  << "\tbubbles_file\tpath to bubbles file\n"
				  << "\tsubs_matrix\tpath to substitution matrix\n"
				  << "\thopo_matrix\tpath to homopolymer matrix\n"
				  << "\tout_file\tpath to output file\n"
				  << "\noptional arguments:\n"
				  << "\t-v verbose_log\tpath to the file "
				  << "with verbose log [default = not set]\n"
				  << "\t-t num_threads\tnumber of parallel threads "
				  << "[default = 1]\n"
				  << "\t-q quiet\tdo not output progress to terminal "
				  << "[default = false]\n";
	};

	numThreads = 1;
	quiet = false;

	const char* optString = "v:t:hq";
	int opt = 0;

	while ((opt = getopt(argc, argv, optString)) != -1)
	{
		switch(opt)
		{
		case 'v':
			outVerbose = optarg;
			break;
		case 't':
			numThreads = atoi(optarg);
			break;
		case 'q':
			quiet = true;
			break;
		case 'h':
			printUsage();
			exit(0);
		}
	}
	if (argc - optind != 4)
	{
		printUsage();
		return false;
	}
	bubblesFile = *(argv + optind);
	scoringMatrix = *(argv + optind + 1);
	hopoMatrix = *(argv + optind + 2);
	outConsensus = *(argv + optind + 3);
	return true;
}

int main(int argc, char* argv[]) 
{
	std::string bubblesFile;
	std::string scoringMatrix;
	std::string hopoMatrix;
	std::string outConsensus;
	std::string outVerbose;
	int  numThreads = 0;
	bool quiet = false;
	if (!parseArgs(argc, argv, bubblesFile, scoringMatrix, 
				   hopoMatrix, outConsensus, outVerbose, numThreads,
				   quiet))
		return 1;

	BubbleProcessor bp(scoringMatrix, hopoMatrix, !quiet);
	if (!outVerbose.empty())
		bp.enableVerboseOutput(outVerbose);
	bp.polishAll(bubblesFile, outConsensus, numThreads); 

	return 0;
}
