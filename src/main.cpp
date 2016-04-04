
#include <iostream>

#include "kmer_index.h"
#include "fasta.h"
#include "overlap.h"
#include "utility.h"
#include "chimera.h"
#include "extender.h"

int main(int argc, char** argv)
{
	static const int MAX_JUMP = 1500;
	static const int MIN_OVERLAP = 7000;
	static const int MAX_OVERHANG = 1500;

	//static const int MIN_KMER_COUNT = 3;
	static const int MIN_KMER_COUNT = 8;
	static const int MAX_KMER_COUNT = 24;
	static const int COVERAGE = 20;

	SequenceContainer& seqContainer = SequenceContainer::getInstance();
	DEBUG_PRINT("Reading FASTA");
	seqContainer.readFasta(argv[1]);
	VertexIndex& vertexIndex = VertexIndex::getInstance();
	vertexIndex.setKmerSize(15);

	DEBUG_PRINT("Building kmer index");
	for (auto rec : seqContainer.getIndex())
	{
		vertexIndex.addFastaSequence(rec.second);
	}

	DEBUG_PRINT("Trimming index");
	vertexIndex.applyKmerThresholds(MIN_KMER_COUNT, MAX_KMER_COUNT);
	DEBUG_PRINT("Building read index");
	vertexIndex.buildReadIndex();

	OverlapDetector ovlp(MAX_JUMP, MIN_OVERLAP, MAX_OVERHANG);
	ovlp.findAllOverlaps(vertexIndex, seqContainer);

	ChimeraDetector chimDetect(MAX_OVERHANG, MAX_JUMP, COVERAGE);
	chimDetect.detectChimeras(ovlp, seqContainer);

	Extender extender(ovlp, chimDetect, seqContainer);
	extender.extendReads();

	return 0;
}
