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
#include "../assemble/extender.h"
#include "../assemble/parameters_estimator.h"
#include "../common/logger.h"
#include "../common/utils.h"
#include "../common/memory_info.h"

#include <getopt.h>

bool parseArgs(int argc, char** argv, std::string& readsFasta, 
			   std::string& outAssembly, std::string& logFile, size_t& genomeSize,
			   int& kmerSize, bool& debug, size_t& numThreads, int& minOverlap, 
			   std::string& configPath, int& minReadLength, bool& unevenCov)
{
	auto printUsage = [argv]()
	{
		std::cerr << "Usage: " << argv[0]
				  << " reads_files out_assembly genome_size config_file\n\t"
				  << "[-r min_read_length] [-l log_file] [-t num_threads]\n\t"
				  << "[-u] [-d] [-h]\n\n"
				  << "positional arguments:\n"
				  << "\treads file\tcomma-separated list of read files\n"
				  << "\tout_assembly\tpath to output file\n"
				  << "\tgenome_size\tgenome size in bytes\n"
				  << "\tconfig_file\tpath to the config file\n"
				  << "\noptional arguments:\n"
				  << "\t-k kmer_size\tk-mer size [default = 15] \n"
				  << "\t-v min_overlap\tminimum overlap between reads "
				  << "[default = 5000] \n"
				  << "\t-d \t\tenable debug output "
				  << "[default = false] \n"
				  << "\t-u \t\tenable uneven coverage (metagenome) mode "
				  << "[default = false] \n"
				  << "\t-l log_file\toutput log to file "
				  << "[default = not set] \n"
				  << "\t-t num_threads\tnumber of parallel threads "
				  << "[default = 1] \n";
	};

	const char optString[] = "l:t:v:k:r:hdu";
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
		case 'u':
			unevenCov = true;
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

int main(int argc, char** argv)
{
	#ifdef NDEBUG
	signal(SIGSEGV, segfaultHandler);
	std::set_terminate(exceptionHandler);
	#endif

	int kmerSize = 15;
	int minReadLength = 0;
	size_t genomeSize = 0;
	int minOverlap = 5000;
	bool debugging = false;
	bool unevenCov = false;
	size_t numThreads = 1;
	std::string readsFasta;
	std::string outAssembly;
	std::string logFile;
	std::string configPath;

	if (!parseArgs(argc, argv, readsFasta, outAssembly, logFile, genomeSize,
				   kmerSize, debugging, numThreads, minOverlap, configPath, 
				   minReadLength, unevenCov)) return 1;

	Logger::get().setDebugging(debugging);
	if (!logFile.empty()) Logger::get().setOutputFile(logFile);
	Logger::get().debug() << "Build date: " << __DATE__ << " " << __TIME__;
	std::ios::sync_with_stdio(false);

	Logger::get().debug() << "Total RAM: " 
		<< getMemorySize() / 1024 / 1024 / 1024 << " Gb";
	Logger::get().debug() << "Available RAM: " 
		<< getFreeMemorySize() / 1024 / 1024 / 1024 << " Gb";
	Logger::get().debug() << "Total CPUs: " << std::thread::hardware_concurrency();

	Config::load(configPath);
	Parameters::get().numThreads = numThreads;
	Parameters::get().kmerSize = kmerSize;
	Parameters::get().minimumOverlap = minOverlap;
	Parameters::get().unevenCoverage = unevenCov;
	Logger::get().debug() << "Running with k-mer size: " << 
		Parameters::get().kmerSize; 
	Logger::get().debug() << "Running with minimum overlap " << minOverlap;
	Logger::get().debug() << "Metagenome mode: " << "NY"[unevenCov];

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

	int64_t sumLength = 0;
	for (auto& seq : readsContainer.iterSeqs())
	{
		sumLength += seq.sequence.length();
	}
	int coverage = sumLength / 2 / genomeSize;
	Logger::get().debug() << "Expected read coverage: " << coverage;

	Logger::get().info() << "Generating solid k-mer index";
	if (!Parameters::get().unevenCoverage)
	{
		size_t hardThreshold = std::min(5, std::max(2, 
				coverage / (int)Config::get("hard_min_coverage_rate")));
		vertexIndex.countKmers(hardThreshold, genomeSize);
	}
	else
	{
		vertexIndex.countKmers(/*hard threshold*/ 2, genomeSize);
	}
	ParametersEstimator estimator(readsContainer, vertexIndex, genomeSize);
	estimator.estimateMinKmerCount();
	int minKmerCov = estimator.minKmerCount();
	vertexIndex.setRepeatCutoff(minKmerCov);
	if (!Parameters::get().unevenCoverage)
	{
		vertexIndex.buildIndex(minKmerCov);
	}
	else
	{
		static const float SELECT_RATE = 0.25;
		static const int TANDEM_FREQ = 10;
		vertexIndex.buildIndexUnevenCoverage(/*min coverage*/ 2, SELECT_RATE, 
											 TANDEM_FREQ);
	}

	Logger::get().debug() << "Peak RAM usage: " 
		<< getPeakRSS() / 1024 / 1024 / 1024 << " Gb";

	//int maxOverlapsNum = !Parameters::get().unevenCoverage ? 5 * coverage : 0;
	OverlapDetector ovlp(readsContainer, vertexIndex,
						 (int)Config::get("maximum_jump"), 
						 Parameters::get().minimumOverlap,
						 (int)Config::get("maximum_overhang"),
						 /*no max overlaps*/ 0, 
						 /*store alignment*/ false,
						 /*only max*/ true,
						 (float)Config::get("assemble_ovlp_divergence"),
						 /* bad end adjustment*/ 0.0f,
						 /* nucl alignent*/ false);
	OverlapContainer readOverlaps(ovlp, readsContainer);

	Extender extender(readsContainer, readOverlaps);
	extender.assembleDisjointigs();
	vertexIndex.clear();

	ConsensusGenerator consGen;
	auto disjointigsFasta = 
		consGen.generateConsensuses(extender.getDisjointigPaths());
	SequenceContainer::writeFasta(disjointigsFasta, outAssembly);

	Logger::get().debug() << "Peak RAM usage: " 
		<< getPeakRSS() / 1024 / 1024 / 1024 << " Gb";

	return 0;
}
