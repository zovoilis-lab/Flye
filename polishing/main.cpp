#include <ctime>
#include <iostream>
#include <getopt.h>

#include "Worker.h"


bool parseArgs(int argc, char** argv, std::string& bubblesFile, 
			   std::string& scoringMatrix, std::string& outConsensus, 
			   std::string& outVerbose)
{
	auto printUsage = []()
	{
		std::cerr << "Usage: polish bubbles_file scoring_matrix out_file "
				  << "[-v verbose_log]\n\n"
				  << "positional arguments:\n"
				  << "\tbubbles_file\tpath to bubbles file\n"
				  << "\tscoring_matrix\tpath to scoring matrix\n"
				  << "\tout_file\tpath to output file\n"
				  << "\noptional arguments:\n"
				  << "\t-v verbose_log\tpath to the file "
				  << "with verbose log [default = not set]\n";
	};

	const char* optString = "v:h";
	int opt = 0;
	while ((opt = getopt(argc, argv, optString)) != -1)
	{
		switch(opt)
		{
		case 'v':
			outVerbose = optarg;
			break;
		case 'h':
			printUsage();
			return false;
		}
	}
	if (argc - optind != 3)
	{
		printUsage();
		return false;
	}
	bubblesFile = *(argv + optind);
	scoringMatrix = *(argv + optind + 1);
	outConsensus = *(argv + optind + 2);
	return true;
}

int main(int argc, char* argv[]) 
{
	std::string bubblesFile;
	std::string scoringMatrix;
	std::string outConsensus;
	std::string outVerbose;
	if (!parseArgs(argc, argv, bubblesFile, scoringMatrix, 
				   outConsensus, outVerbose))
		return 1;

	Worker worker(scoringMatrix);
	worker.run(bubblesFile); 
	worker.writeConsensuses(outConsensus);
	if (!outVerbose.empty())
		worker.writeLog(outVerbose);

	return 0;
}
