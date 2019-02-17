//(c) 2016 by Authors
//This file is a part of ABruijn program.
//Released under the BSD license (see LICENSE file)

#include <iostream>
#include <signal.h>
#include <stdlib.h>
#include <unistd.h>
#include <cmath>

#include "../sequence/sequence_container.h"
#include "../common/config.h"
#include "../common/logger.h"
#include "../common/utils.h"
#include "../common/memory_info.h"

#include "../repeat_graph/repeat_graph.h"
#include "../repeat_graph/multiplicity_inferer.h"
#include "../repeat_graph/graph_processing.h"
#include "../repeat_graph/repeat_resolver.h"
#include "../repeat_graph/output_generator.h"

#include <getopt.h>

bool parseArgs(int argc, char** argv, std::string& readsFasta, 
			   std::string& outFolder, std::string& logFile, 
			   std::string& inAssembly, int& kmerSize,
			   int& minOverlap, bool& debug, size_t& numThreads, 
			   std::string& configPath, bool& unevenCov)
{
	auto printUsage = [argv]()
	{
		std::cerr << "Usage: " << argv[0]
				  << "\tin_assembly reads_files out_folder config_path\n\t"
				  << "[-l log_file] [-t num_threads] [-v min_overlap]\n\t"
				  << "[-d] [-u]\n\n"
				  << "positional arguments:\n"
				  << "\tin_assembly\tpath to input assembly\n"
				  << "\treads_files\tcomma-separated list with reads\n"
				  << "\tout_assembly\tpath to output assembly\n"
				  << "\tconfig_path\tpath to config file\n"
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

	const char optString[] = "l:t:v:k:hdu";
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

	inAssembly = *(argv + optind);
	readsFasta = *(argv + optind + 1);
	outFolder = *(argv + optind + 2);
	configPath = *(argv + optind + 3);

	return true;
}

int main(int argc, char** argv)
{
	#ifdef NDEBUG
	signal(SIGSEGV, segfaultHandler);
	std::set_terminate(exceptionHandler);
	#endif

	bool debugging = false;
	size_t numThreads = 1;
	int kmerSize = 15;
	int minOverlap = 5000;
	bool unevenCov = false;
	std::string readsFasta;
	std::string inAssembly;
	std::string outFolder;
	std::string logFile;
	std::string configPath;
	if (!parseArgs(argc, argv, readsFasta, outFolder, logFile, inAssembly,
				   kmerSize, minOverlap, debugging, 
				   numThreads, configPath, unevenCov))  return 1;
	
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
	Logger::get().debug() << "Selected minimum overlap " << minOverlap;
	Logger::get().debug() << "Metagenome mode: " << "NY"[unevenCov];

	Logger::get().info() << "Reading sequences";
	SequenceContainer seqAssembly; 
	SequenceContainer seqReads;
	std::vector<std::string> readsList = splitString(readsFasta, ',');
	try
	{
		seqAssembly.loadFromFile(inAssembly);
		for (auto& readsFile : readsList)
		{
			seqReads.loadFromFile(readsFile);
		}
	}
	catch (SequenceContainer::ParseException& e)
	{
		Logger::get().error() << e.what();
		return 1;
	}
	seqReads.buildPositionIndex();
	seqAssembly.buildPositionIndex();

	SequenceContainer edgeSequences;
	RepeatGraph rg(seqAssembly, &edgeSequences);
	GraphProcessor proc(rg, seqAssembly);
	ReadAligner aligner(rg, seqReads);
	OutputGenerator outGen(rg, aligner, seqReads);

	Logger::get().info() << "Building repeat graph";
	rg.build();
	//outGen.outputDot(proc.getEdgesPaths(), outFolder + "/graph_raw.gv");

	Logger::get().info() << "Aligning reads to the graph";
	aligner.alignReads();

	MultiplicityInferer multInf(rg, aligner, seqAssembly, seqReads);
	multInf.estimateCoverage();

	//remove edges/connections with low coverage
	multInf.removeUnsupportedEdges();
	multInf.removeUnsupportedConnections();

	//collapse graph structures cause by heterogenity
	multInf.collapseHeterozygousLoops();
	multInf.collapseHeterozygousBulges();

	Logger::get().info() << "Resolving repeats";
	RepeatResolver resolver(rg, seqAssembly, seqReads, aligner, multInf);
	resolver.findRepeats();
	outGen.outputDot(proc.getEdgesPaths(), outFolder + "/graph_before_rr.gv");
	//outGen.outputGfa(proc.getEdgesPaths(), outFolder + "/graph_before_rr.gfa");
	outGen.outputFasta(proc.getEdgesPaths(), outFolder + "/graph_before_rr.fasta");
	resolver.resolveRepeats();

	//collapse again after repeat resolution
	multInf.collapseHeterozygousLoops();
	multInf.collapseHeterozygousBulges();
	//do tip trimming only after repeat resolution, since those
	//tips might actually help to resolve some repeats,
	//and help to identify repeat boundaries
	multInf.trimTips();

	resolver.findRepeats();
	resolver.finalizeGraph();

	outGen.outputDot(proc.getEdgesPaths(), outFolder + "/graph_after_rr.gv");
	rg.storeGraph(outFolder + "/repeat_graph_dump");
	aligner.storeAlignments(outFolder + "/read_alignment_dump");
	SequenceContainer::writeFasta(edgeSequences.iterSeqs(), 
								  outFolder + "/repeat_graph_edges.fasta",
								  /*only pos strand*/ true);

	Logger::get().debug() << "Peak RAM usage: " 
		<< getPeakRSS() / 1024 / 1024 / 1024 << " Gb";

	return 0;
}
