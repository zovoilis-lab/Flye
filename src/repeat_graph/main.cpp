//(c) 2016 by Authors
//This file is a part of ABruijn program.
//Released under the BSD license (see LICENSE file)

#include <iostream>
#include <getopt.h>

#include "../sequence/vertex_index.h"
#include "../sequence/sequence_container.h"
#include "../sequence/overlap.h"
#include "../common/config.h"
#include "../common/logger.h"
#include "repeat_graph.h"
#include "graph_processing.h"
#include "repeat_resolver.h"

bool parseArgs(int argc, char** argv, std::string& readsFasta, 
			   std::string& outFolder, std::string& logFile, 
			   std::string& inAssembly, int& kmerSize,
			   int& minOverlap, bool& debug, size_t& numThreads)
{
	auto printUsage = [argv]()
	{
		std::cerr << "Usage: " << argv[0]
				  << "\tin_assembly reads_file out_folder \n\t\t\t\t"
				  << "[-l log_file] [-t num_threads] [-v min_overlap]\n\t\t\t\t"
				  << "[-k kmer_size] [-d]\n\n"
				  << "positional arguments:\n"
				  << "\tin_assembly\tpath to input assembly\n"
				  << "\treads file\tpath to fasta with reads\n"
				  << "\tout_assembly\tpath to output assembly\n"
				  << "\noptional arguments:\n"
				  << "\t-k kmer_size\tk-mer size [default = 15] \n"
				  << "\t-v min_overlap\tminimum overlap between reads "
				  << "[default = 5000] \n"
				  << "\t-d \t\tenable debug output "
				  << "[default = false] \n"
				  << "\t-l log_file\toutput log to file "
				  << "[default = not set] \n"
				  << "\t-t num_threads\tnumber of parallel threads "
				  << "[default = 1] \n";
	};

	numThreads = 1;
	debug = false;
	logFile = "";
	minOverlap = 5000;
	kmerSize = 15;

	const char optString[] = "l:t:k:v:hd";
	int opt = 0;
	while ((opt = getopt(argc, argv, optString)) != -1)
	{
		switch(opt)
		{
		case 't':
			numThreads = atoi(optarg);
			break;
		case 'v':
			minOverlap = atoi(optarg);
			break;
		case 'k':
			kmerSize = atoi(optarg);
			break;
		case 'l':
			logFile = optarg;
			break;
		case 'd':
			debug = true;
			break;
		case 'h':
			printUsage();
			exit(0);
		}
	}
	if (argc - optind != 3)
	{
		printUsage();
		return false;
	}

	inAssembly = *(argv + optind);
	readsFasta = *(argv + optind + 1);
	outFolder = *(argv + optind + 2);

	return true;
}

bool fileExists(const std::string& path)
{
	std::ifstream fin(path);
	return fin.good();
}

int main(int argc, char** argv)
{
	bool debugging = false;
	size_t numThreads;
	int kmerSize;
	int minOverlap;
	std::string readsFasta;
	std::string inAssembly;
	std::string outFolder;
	std::string logFile;

	if (!parseArgs(argc, argv, readsFasta, outFolder, logFile, inAssembly,
				   kmerSize, minOverlap, debugging, numThreads)) 
	{
		return 1;
	}
	Parameters::get().minimumOverlap = minOverlap;
	Parameters::get().kmerSize = kmerSize;
	Parameters::get().numThreads = numThreads;

	try
	{
		Logger::get().setDebugging(debugging);
		if (!logFile.empty()) Logger::get().setOutputFile(logFile);

		Logger::get().debug() << "Build date: " << __DATE__ << " " << __TIME__;

		SequenceContainer seqAssembly; 
		seqAssembly.readFasta(inAssembly);
		SequenceContainer seqReads;
		seqReads.readFasta(readsFasta);

		RepeatGraph rg(seqAssembly);
		rg.build();
		rg.outputDot(outFolder + "/graph_before.dot");

		GraphProcessor proc(rg, seqAssembly, seqReads);
		proc.simplify();

		RepeatResolver resolver(rg, seqAssembly, seqReads);
		resolver.alignReads();
		resolver.correctEdgesMultiplicity();

		rg.outputDot(outFolder + "/graph_simplified.dot");

		resolver.resolveRepeats();
		rg.outputDot(outFolder + "/graph_after.dot");

		proc.generateContigs();
		proc.outputContigsGraph(outFolder + "/graph_condensed.dot");
		proc.outputContigsFasta(outFolder + "/graph_edges.fasta");

		return 0;

	}
	catch (std::runtime_error& e)
	{
		Logger::get().error() << e.what();
		return 1;
	}

	return 0;
}
