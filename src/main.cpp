
#include <iostream>

#include "kmer_index.h"
#include "fasta.h"
#include "overlap.h"
#include "utility.h"

int main(int argc, char** argv)
{
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
	vertexIndex.applyKmerThresholds(8, 24);
	DEBUG_PRINT("Building read index");
	vertexIndex.buildReadIndex();

	OverlapDetector ovlp(1500, 7000, 1500);
	ovlp.findAllOverlaps(vertexIndex, seqContainer);

	//vertexIndex.outputCounts();
	return 0;
}
