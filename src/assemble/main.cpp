//(c) 2016 by Authors
//This file is a part of ABruijn program.
//Released under the BSD license (see LICENSE file)

#include <iostream>
#include <signal.h>
#include <stdlib.h>
#include <unistd.h>
#include <cmath>
#include <execinfo.h>

#include "../sequence/vertex_index.h"
#include "../sequence/sequence_container.h"
#include "../sequence/overlap.h"
#include "../sequence/consensus_generator.h"
#include "../common/config.h"
#include "chimera.h"
#include "extender.h"
#include "parameters_estimator.h"
#include "../common/logger.h"
#include "../common/utils.h"

#include <getopt.h>

bool parseArgs(int argc, char** argv, std::string& readsFasta, 
			   std::string& outAssembly, std::string& logFile, size_t& genomeSize,
			   int& minKmer, int& maxKmer, int& kmerSize, bool& debug,
			   size_t& numThreads, int& minOverlap, std::string& configPath,
			   int& minReadLength)
{
	auto printUsage = [argv]()
	{
		std::cerr << "Usage: " << argv[0]
				  << "\treads_files out_assembly genome_size config_file\n\t"
				  << "[-m min_kmer_cov] [-r min_read_length]\n\t"
				  << "[-x max_kmer_cov] [-l log_file] [-t num_threads] [-d]\n\n"
				  << "positional arguments:\n"
				  << "\treads file\tcomma-separated list of read files\n"
				  << "\tout_assembly\tpath to output file\n"
				  << "\tgenome_size\tgenome size in bytes\n"
				  << "\tconfig_file\tpath to the config file\n"
				  << "\noptional arguments:\n"
				  << "\t-k kmer_size\tk-mer size [default = 15] \n"
				  << "\t-m min_kmer_cov\tminimum k-mer coverage "
				  << "[default = auto] \n"
				  << "\t-x max_kmer_cov\tmaximum k-mer coverage "
				  << "[default = auto] \n"
				  << "\t-v min_overlap\tminimum overlap between reads "
				  << "[default = 5000] \n"
				  << "\t-d \t\tenable debug output "
				  << "[default = false] \n"
				  << "\t-l log_file\toutput log to file "
				  << "[default = not set] \n"
				  << "\t-t num_threads\tnumber of parallel threads "
				  << "[default = 1] \n";
	};

	//kmerSize = -1;
	minKmer = -1;
	maxKmer = -1;
	kmerSize = 15;
	numThreads = 1;
	debug = false;
	minOverlap = -1;
	minReadLength = 0;

	const char optString[] = "m:x:l:t:v:k:r:hds";
	int opt = 0;
	while ((opt = getopt(argc, argv, optString)) != -1)
	{
		switch(opt)
		{
		case 'k':
			kmerSize = atoi(optarg);
			break;
		case 'r':
			minReadLength = atoi(optarg);
			break;
		case 'm':
			minKmer = atoi(optarg);
			break;
		case 'x':
			maxKmer = atoi(optarg);
			break;
		case 't':
			numThreads = atoi(optarg);
			break;
		case 'v':
			minOverlap = atoi(optarg);
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
	if (argc - optind != 4)
	{
		printUsage();
		return false;
	}
	readsFasta = *(argv + optind);
	outAssembly = *(argv + optind + 1);
	genomeSize = atoll(*(argv + optind + 2));
	configPath = *(argv + optind + 3);

	return true;
}

bool fileExists(const std::string& path)
{
	std::ifstream fin(path);
	return fin.good();
}

void segfaultHandler(int signal)
{
	void *stackArray[20];
	size_t size = backtrace(stackArray, 10);
	Logger::get().error() << "Segmentation fault! Backtrace:";
	char** backtrace = backtrace_symbols(stackArray, size);
	for (size_t i = 0; i < size; ++i)
	{
		Logger::get().error() << "\t" << backtrace[i];
	}
	exit(1);
}

void exceptionHandler()
{
	static bool triedThrow = false;
	try
	{
        if (!triedThrow)
		{
			triedThrow = true;
			throw;
		}
    }
    catch (const std::exception &e) 
	{
        Logger::get().error() << "Caught unhandled exception: " << e.what();
    }
	catch (...) {}

	void *stackArray[20];
	size_t size = backtrace(stackArray, 10);
	char** backtrace = backtrace_symbols(stackArray, size);
	for (size_t i = 0; i < size; ++i)
	{
		Logger::get().error() << "\t" << backtrace[i];
	}
	exit(1);
}

int main(int argc, char** argv)
{
	#ifndef _DEBUG
	signal(SIGSEGV, segfaultHandler);
	std::set_terminate(exceptionHandler);
	#endif

	int minKmerCov = 0;
	int maxKmerCov = 0;
	int kmerSize = 0;
	int minReadLength = 0;
	size_t genomeSize = 0;
	int minOverlap = 5000;
	bool debugging = false;
	size_t numThreads = 1;
	std::string readsFasta;
	std::string outAssembly;
	std::string logFile;
	std::string configPath;

	if (!parseArgs(argc, argv, readsFasta, outAssembly, logFile, genomeSize,
				   minKmerCov, maxKmerCov, kmerSize, debugging, numThreads,
				   minOverlap, configPath, minReadLength)) return 1;

	Logger::get().setDebugging(debugging);
	if (!logFile.empty()) Logger::get().setOutputFile(logFile);
	Logger::get().debug() << "Build date: " << __DATE__ << " " << __TIME__;
	std::ios::sync_with_stdio(false);

	Config::load(configPath);
	Parameters::get().numThreads = numThreads;
	Parameters::get().kmerSize = kmerSize;
	Logger::get().debug() << "Running with k-mer size: " << 
		Parameters::get().kmerSize; 

	SequenceContainer readsContainer;
	std::vector<std::string> readsList = splitString(readsFasta, ',');
	Logger::get().info() << "Reading sequences";
	try
	{
		for (auto& readsFile : readsList)
		{
			readsContainer.loadFromFile(readsFile, minReadLength);
		}
	}
	catch (SequenceContainer::ParseException& e)
	{
		Logger::get().error() << e.what();
		return 1;
	}
	readsContainer.buildPositionIndex();
	VertexIndex vertexIndex(readsContainer, 
							(int)Config::get("assemble_kmer_sample"));
	vertexIndex.outputProgress(true);

	Parameters::get().minimumOverlap = minOverlap;
	Logger::get().debug() << "Running with minimum overlap " << minOverlap;

	int64_t sumLength = 0;
	for (auto& seq : readsContainer.iterSeqs())
	{
		sumLength += seq.sequence.length();
	}
	int coverage = sumLength / 2 / genomeSize;
	Logger::get().debug() << "Expected read coverage: " << coverage;

	Logger::get().info() << "Generating solid k-mer index";
	//size_t hardThreshold = std::min(5, std::max(2, 
	//		coverage / (int)Config::get("hard_min_coverage_rate")));
	vertexIndex.countKmers(/*hard threshold*/ 2, genomeSize);

	ParametersEstimator estimator(readsContainer, vertexIndex, genomeSize);
	estimator.estimateMinKmerCount();
	if (minKmerCov == -1)
	{
		minKmerCov = estimator.minKmerCount();
	}

	vertexIndex.setRepeatCutoff(minKmerCov);
	//vertexIndex.buildIndex(minKmerCov);
	vertexIndex.buildIndexUnevenCoverage(/*min coverage*/ 2);

	OverlapDetector ovlp(readsContainer, vertexIndex,
						 (int)Config::get("maximum_jump"), 
						 Parameters::get().minimumOverlap,
						 (int)Config::get("maximum_overhang"),
						 /*no max ovlp*/ 0, 
						 /*store alignment*/ false,
						 /*only max*/ true,
						 (float)Config::get("assemble_ovlp_ident"), 
						 /* bad end adjustment*/ 0.0f);
	OverlapContainer readOverlaps(ovlp, readsContainer);

	Extender extender(readsContainer, readOverlaps);
	extender.assembleContigs();
	vertexIndex.clear();

	ConsensusGenerator consGen;
	auto contigsFasta = consGen.generateConsensuses(extender.getContigPaths());
	SequenceContainer::writeFasta(contigsFasta, outAssembly);

	return 0;
}
