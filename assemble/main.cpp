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

bool parseArgs(int argc, char** argv, std::string& readsFasta, 
			   std::string& outAssembly, int& kmerSize, int& minKmer,
			   int& maxKmer)
{
	auto printUsage = [argv]()
	{
		std::cerr << "Usage: " << argv[0]
				  << " reads_file out_assembly "
				  << "[-k kmer_size] [-m min_kmer_cov] "
				  << "[-x max_kmer_cov] \n\n"
				  << "positional arguments:\n"
				  << "\treads file\tpath to fasta with reads\n"
				  << "\tout_assembly\tpath to output file\n"
				  << "\noptional arguments:\n"
				  << "\t-k kmer_size\tk-mer size [default = 15] \n"
				  << "\t-m min_kmer_cov\tminimum k-mer coverage "
				  << "[default = auto] \n"
				  << "\t-x max_kmer_cov\tmaximum k-mer coverage "
				  << "[default = 24] \n";
	};

	kmerSize = 15;
	minKmer = -1;
	maxKmer = 24;

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
	if (argc - optind != 2)
	{
		printUsage();
		return false;
	}
	readsFasta = *(argv + optind);
	outAssembly = *(argv + optind + 1);
	return true;
}

int main(int argc, char** argv)
{
	static const int MAX_JUMP = 1500;
	static const int MIN_OVERLAP = 7000;
	static const int MAX_OVERHANG = 1500;

	int kmerSize = 0;
	int minKmerCov = 0;
	int maxKmerCov = 0;
	std::string readsFasta;
	std::string outAssembly;

	if (!parseArgs(argc, argv, readsFasta, outAssembly, 
				   kmerSize, minKmerCov, maxKmerCov))
	{
		return 1;
	}

	SequenceContainer& seqContainer = SequenceContainer::getInstance();
	LOG_PRINT("Reading FASTA");
	seqContainer.readFasta(readsFasta);
	VertexIndex& vertexIndex = VertexIndex::getInstance();
	vertexIndex.setKmerSize(kmerSize);

	LOG_PRINT("Building kmer index");
	vertexIndex.buildKmerIndex(seqContainer);

	LOG_PRINT("Trimming index");
	if (minKmerCov == -1)
	{
		minKmerCov = vertexIndex.estimateCoverageCutoff();
	}
	vertexIndex.applyKmerThresholds(minKmerCov, maxKmerCov);
	LOG_PRINT("Building read index");
	vertexIndex.buildReadIndex();

	OverlapDetector ovlp(MAX_JUMP, MIN_OVERLAP, MAX_OVERHANG);
	ovlp.findAllOverlaps(vertexIndex, seqContainer);

	ChimeraDetector chimDetect(MAX_OVERHANG, MAX_JUMP);
	chimDetect.detectChimeras(ovlp, seqContainer);

	Extender extender(ovlp, chimDetect, seqContainer);
	extender.extendReads();

	ContigGenerator contGen(MAX_JUMP, extender, ovlp, 
							vertexIndex, seqContainer);
	contGen.generateContigs();
	contGen.outputContigs(outAssembly);

	return 0;
}
