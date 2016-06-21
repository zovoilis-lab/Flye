//(c) 2016 by Authors
//This file is a part of ABruijn program.
//Released under the BSD license (see LICENSE file)

#include <iostream>
#include <getopt.h>

#include "vertex_index.h"
#include "sequence_container.h"
#include "overlap.h"
#include "chimera.h"
#include "extender.h"
#include "contig.h"
#include "parameters_estimator.h"
#include "logger.h"

bool parseArgs(int argc, char** argv, std::string& readsFasta, 
			   std::string& outAssembly, std::string& logFile, int& coverage,
			   int& kmerSize, int& minKmer, int& maxKmer, bool& debug,
			   size_t& numThreads)
{
	auto printUsage = [argv]()
	{
		std::cerr << "Usage: " << argv[0]
				  << "\treads_file out_assembly coverage \n\t\t\t\t"
				  << "[-k kmer_size] [-m min_kmer_cov] \n\t\t\t\t"
				  << "[-x max_kmer_cov] [-l log_file] [-t num_threads] [-d]\n\n"
				  << "positional arguments:\n"
				  << "\treads file\tpath to fasta with reads\n"
				  << "\tout_assembly\tpath to output file\n"
				  << "\tcoverage\testimated assembly coverage\n"
				  << "\noptional arguments:\n"
				  << "\t-k kmer_size\tk-mer size [default = 15] \n"
				  << "\t-m min_kmer_cov\tminimum k-mer coverage "
				  << "[default = auto] \n"
				  << "\t-x max_kmer_cov\tmaximum k-mer coverage "
				  << "[default = auto] \n"
				  << "\t-d \t\tenable debug output "
				  << "[default = false] \n"
				  << "\t-l log_file\toutput log to file "
				  << "[default = not set] \n"
				  << "\t-t num_threads\tnumber of parallel threads "
				  << "[default = 1] \n";
	};

	kmerSize = 15;
	minKmer = -1;
	maxKmer = -1;
	numThreads = 1;
	debug = false;
	logFile = "";

	const char optString[] = "k:m:x:l:t:hd";
	int opt = 0;
	while ((opt = getopt(argc, argv, optString)) != -1)
	{
		switch(opt)
		{
		case 'k':
			kmerSize = atoi(optarg);
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
	readsFasta = *(argv + optind);
	outAssembly = *(argv + optind + 1);
	coverage = atoi(*(argv + optind + 2));

	return true;
}

int main(int argc, char** argv)
{
	static const int MAX_JUMP = 1500;
	static const int MIN_OVERLAP = 5000;
	static const int MAX_OVERHANG = 1500;
	static const int MAGIC_10 = 10;

	int kmerSize = 0;
	int minKmerCov = 0;
	int maxKmerCov = 0;
	int coverage = 0;
	bool debugging = false;
	size_t numThreads;
	std::string readsFasta;
	std::string outAssembly;
	std::string logFile;

	if (!parseArgs(argc, argv, readsFasta, outAssembly, logFile, coverage,
				   kmerSize, minKmerCov, maxKmerCov, debugging, numThreads)) return 1;

	try
	{
		Logger::get().setDebugging(debugging);
		if (!logFile.empty()) Logger::get().setOutputFile(logFile);

		SequenceContainer& seqContainer = SequenceContainer::get();
		Logger::get().debug() << "Reading FASTA";
		seqContainer.readFasta(readsFasta);
		VertexIndex& vertexIndex = VertexIndex::get();
		vertexIndex.setKmerSize(kmerSize);

		//rough estimate
		size_t hardThreshold = std::max(1, coverage / MAGIC_10);
		//vertexIndex.buildKmerIndex(seqContainer, hardThreshold);
		vertexIndex.countKmers(seqContainer, hardThreshold);

		if (maxKmerCov == -1)
		{
			maxKmerCov = MAGIC_10 * coverage;
		}
		if (minKmerCov == -1)
		{
			ParametersEstimator estimator;
			minKmerCov = estimator.estimateMinKmerCount(coverage, maxKmerCov);
		}

		vertexIndex.buildIndex(minKmerCov, maxKmerCov);

		OverlapDetector ovlp(MAX_JUMP, MIN_OVERLAP, MAX_OVERHANG);
		ovlp.findAllOverlaps(numThreads);

		ChimeraDetector chimDetect(MAX_OVERHANG, MAX_JUMP, MIN_OVERLAP,
								   coverage, ovlp);
		chimDetect.detectChimeras();

		Extender extender(ovlp, chimDetect, MAX_JUMP, coverage);
		extender.assembleContigs();

		ContigGenerator contGen(MAX_JUMP, extender, ovlp);
		contGen.generateContigs();
		contGen.outputContigs(outAssembly);
	}
	catch (std::runtime_error& e)
	{
		Logger::get().error() << e.what();
		return 1;
	}

	return 0;
}
