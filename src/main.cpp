
#include <iostream>

#include "kmer_index.h"
#include "fasta.h"
#include "overlap.h"
#include "utility.h"
#include "chimera.h"
#include "extender.h"
#include "contig.h"

int main(int argc, char** argv)
{
	static const int MAX_JUMP = 1500;
	static const int MIN_OVERLAP = 7000;
	static const int MAX_OVERHANG = 1500;

	//static const int MIN_KMER_COUNT = 3;
	static const int MIN_KMER_COUNT = 8;
	static const int MAX_KMER_COUNT = 24;
	//static const int COVERAGE = 10;
	static const int COVERAGE = 20;

	SequenceContainer& seqContainer = SequenceContainer::getInstance();
	std::cerr << "Reading FASTA\n";
	seqContainer.readFasta(argv[1]);
	VertexIndex& vertexIndex = VertexIndex::getInstance();
	vertexIndex.setKmerSize(15);

	std::cerr << "Building kmer index\n";
	for (auto rec : seqContainer.getIndex())
	{
		vertexIndex.addFastaSequence(rec.second);
	}

	std::cerr << "Trimming index\n";
	vertexIndex.applyKmerThresholds(MIN_KMER_COUNT, MAX_KMER_COUNT);
	std::cerr << "Building read index\n";
	vertexIndex.buildReadIndex();

	OverlapDetector ovlp(MAX_JUMP, MIN_OVERLAP, MAX_OVERHANG);
	ovlp.findAllOverlaps(vertexIndex, seqContainer);

	ChimeraDetector chimDetect(MAX_OVERHANG, MAX_JUMP, COVERAGE);
	chimDetect.detectChimeras(ovlp, seqContainer);

	Extender extender(ovlp, chimDetect, seqContainer);
	extender.extendReads();

	ContigGenerator contGen(MAX_JUMP, extender, ovlp, 
							vertexIndex, seqContainer);
	contGen.generateContigs();
	contGen.outputContigs("genome.fasta");

	return 0;
}
