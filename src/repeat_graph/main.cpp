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
#include "../common/config.h"
#include "../common/logger.h"
#include "repeat_graph.h"
#include "multiplicity_inferer.h"
#include "graph_processing.h"
#include "repeat_resolver.h"
#include "output_generator.h"
#include "contig_extender.h"

#include <getopt.h>

bool parseArgs(int argc, char** argv, std::string& readsFasta, 
			   std::string& outFolder, std::string& logFile, 
			   std::string& inAssembly, int& kmerSize,
			   int& minOverlap, bool& debug, size_t& numThreads, 
			   std::string& configPath, size_t& genomeSize)
{
	auto printUsage = [argv]()
	{
		std::cerr << "Usage: " << argv[0]
				  << "\tin_assembly reads_files out_folder genome_size config_path\n\t"
				  << "[-l log_file] [-t num_threads] [-v min_overlap]\n\t"
				  << "[-d]\n\n"
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
				  << "\t-l log_file\toutput log to file "
				  << "[default = not set] \n"
				  << "\t-t num_threads\tnumber of parallel threads "
				  << "[default = 1] \n";
	};

	numThreads = 1;
	debug = false;
	minOverlap = -1;
	kmerSize = 15;

	const char optString[] = "l:t:v:k:hd";
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
	if (argc - optind != 5)
	{
		printUsage();
		return false;
	}

	inAssembly = *(argv + optind);
	readsFasta = *(argv + optind + 1);
	outFolder = *(argv + optind + 2);
	genomeSize = atoll(*(argv + optind + 3));
	configPath = *(argv + optind + 4);

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

	bool debugging = false;
	size_t numThreads = 1;
	int kmerSize = 15;
	int minOverlap = 5000;
	size_t genomeSize = 0;
	std::string readsFasta;
	std::string inAssembly;
	std::string outFolder;
	std::string logFile;
	std::string configPath;
	if (!parseArgs(argc, argv, readsFasta, outFolder, logFile, inAssembly,
				   kmerSize, minOverlap, debugging, 
				   numThreads, configPath, genomeSize))  return 1;
	
	Logger::get().setDebugging(debugging);
	if (!logFile.empty()) Logger::get().setOutputFile(logFile);
	Logger::get().debug() << "Build date: " << __DATE__ << " " << __TIME__;
	std::ios::sync_with_stdio(false);
	
	Config::load(configPath);
	Parameters::get().numThreads = numThreads;
	Parameters::get().kmerSize = kmerSize;
	Logger::get().debug() << "Running with k-mer size: " << 
		Parameters::get().kmerSize; 

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

	Parameters::get().minimumOverlap = minOverlap;
	Logger::get().debug() << "Selected minimum overlap " << minOverlap;

	RepeatGraph rg(seqAssembly);
	GraphProcessor proc(rg, seqAssembly, seqReads);
	ReadAligner aligner(rg, seqAssembly, seqReads);
	OutputGenerator outGen(rg, aligner, seqAssembly, seqReads);

	Logger::get().info() << "Building repeat graph";
	rg.build();
	//outGen.outputDot(proc.getEdgesPaths(), outFolder + "/graph_raw.gv");
	proc.simplify();

	Logger::get().info() << "Aligning reads to the graph";
	aligner.alignReads();

	MultiplicityInferer multInf(rg, aligner, seqAssembly, seqReads);
	multInf.estimateCoverage();
	multInf.removeUnsupportedEdges();
	multInf.removeUnsupportedConnections();

	//for diploid genomes, turned off by default
	multInf.collapseHeterozygousLoops();
	multInf.collapseHeterozygousBulges();

	Logger::get().info() << "Resolving repeats";
	RepeatResolver resolver(rg, seqAssembly, seqReads, aligner, multInf);
	resolver.findRepeats();

	outGen.outputDot(proc.getEdgesPaths(), outFolder + "/graph_before_rr.gv");
	//outGen.outputGfa(proc.getEdgesPaths(), outFolder + "/graph_before_rr.gfa");
	outGen.outputFasta(proc.getEdgesPaths(), outFolder + "/graph_before_rr.fasta");
	//outGen.detailedFasta(outFolder + "/before_rr_detailed.fasta");

	resolver.resolveRepeats();
	resolver.fixLongEdges();
	//outGen.outputDot(proc.getEdgesPaths(), outFolder + "/graph_after_rr.gv");

	Logger::get().info() << "Generating contigs";

	ContigExtender extender(rg, aligner, seqAssembly, seqReads, 
							multInf.getMeanCoverage());
	extender.generateUnbranchingPaths();
	extender.generateContigs();
	extender.outputContigs(outFolder + "/graph_paths.fasta");
	extender.outputStatsTable(outFolder + "/contigs_stats.txt");
	extender.outputScaffoldConnections(outFolder + "/scaffolds_links.txt");

	outGen.dumpRepeats(extender.getUnbranchingPaths(),
					   outFolder + "/repeats_dump.txt");
	outGen.outputDot(extender.getUnbranchingPaths(),
					 outFolder + "/graph_final.gv");
	outGen.outputFasta(extender.getUnbranchingPaths(),
					   outFolder + "/graph_final.fasta");
	outGen.outputGfa(extender.getUnbranchingPaths(),
					 outFolder + "/graph_final.gfa");
	return 0;
}
