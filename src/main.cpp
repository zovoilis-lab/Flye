
#include <iostream>

#include "kmer_index.h"
#include "fasta.h"
#include "overlap.h"

int main(int argc, char** argv)
{
	SequenceContainer& seqContainer = SequenceContainer::getInstance();
	seqContainer.readFasta(argv[1]);
	VertexIndex& vertexIndex = VertexIndex::getInstance();
	vertexIndex.setKmerSize(15);

	{
		for (auto rec : seqContainer.getIndex())
		{
			vertexIndex.addFastaSequence(rec.second);
		}
	}

	vertexIndex.applyKmerThresholds(8, 24);
	vertexIndex.outputCounts();
	vertexIndex.buildReadIndex();

	OverlapDetector ovlp(1000, 1000, 1000);
	ovlp.findAllOverlaps(vertexIndex, seqContainer);

	return 0;
}
