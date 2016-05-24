//(c) 2016 by Authors
//This file is a part of ABruijn program.
//Released under the BSD license (see LICENSE file)

#include <iostream>
#include <getopt.h>

#include "kmer_index.h"
#include "fasta.h"
#include "overlap.h"
#include "utility.h"
#include "chimera.h"
#include "extender.h"
#include "contig.h"
#include "parameters_estimator.h"

bool parseArgs(int argc, char** argv, std::string& readsFasta, 
			   std::string& outAssembly, int& coverage,
			   int& kmerSize, int& minKmer, int& maxKmer)
{
	auto printUsage = [argv]()
	{
		std::cerr << "Usage: " << argv[0]
				  << " reads_file out_assembly coverage "
				  << "[-k kmer_size] [-m min_kmer_cov] "
				  << "[-x max_kmer_cov] \n\n"
				  << "positional arguments:\n"
				  << "\treads file\tpath to fasta with reads\n"
				  << "\tout_assembly\tpath to output file\n"
				  << "\tcoverage\testimated assembly coverage\n"
				  << "\noptional arguments:\n"
				  << "\t-k kmer_size\tk-mer size [default = 15] \n"
				  << "\t-m min_kmer_cov\tminimum k-mer coverage "
				  << "[default = auto] \n"
				  << "\t-x max_kmer_cov\tmaximum k-mer coverage "
				  << "[default = not set] \n";
	};

	kmerSize = 15;
	minKmer = -1;
	maxKmer = -1;

	const char optString[] = "k:m:x:h";
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
	std::string readsFasta;
	std::string outAssembly;

	if (!parseArgs(argc, argv, readsFasta, outAssembly, coverage,
				   kmerSize, minKmerCov, maxKmerCov))
	{
		return 1;
	}

	try
	{
		SequenceContainer& seqContainer = SequenceContainer::getInstance();
		LOG_PRINT("Reading FASTA");
		seqContainer.readFasta(readsFasta);
		VertexIndex& vertexIndex = VertexIndex::getInstance();
		vertexIndex.setKmerSize(kmerSize);

		//rough estimate
		size_t hardThreshold = std::max(1, coverage / MAGIC_10);
		vertexIndex.buildKmerIndex(seqContainer, hardThreshold);

		if (maxKmerCov == -1)
		{
			maxKmerCov = MAGIC_10 * coverage;
		}
		if (minKmerCov == -1)
		{
			ParametersEstimator estimator(vertexIndex, seqContainer);
			minKmerCov = estimator.estimateMinKmerCount(coverage, maxKmerCov);
		}

		vertexIndex.applyKmerThresholds(minKmerCov, maxKmerCov);
		vertexIndex.buildReadIndex();

		OverlapDetector ovlp(MAX_JUMP, MIN_OVERLAP, MAX_OVERHANG,
							 vertexIndex, seqContainer);
		ovlp.findAllOverlaps();

		ChimeraDetector chimDetect(MAX_OVERHANG, MAX_JUMP, MIN_OVERLAP,
								   coverage, ovlp, seqContainer);
		chimDetect.detectChimeras();

		Extender extender(ovlp, chimDetect, seqContainer, MAX_JUMP);
		extender.assembleContigs();

		ContigGenerator contGen(MAX_JUMP, extender, ovlp, 
								vertexIndex, seqContainer);
		contGen.generateContigs();
		contGen.outputContigs(outAssembly);
	}
	catch (std::runtime_error& e)
	{
		ERROR_PRINT(e.what());
		return 1;
	}

	return 0;
}
