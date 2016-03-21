
#include "kmer_index.h"
#include "fasta.h"
#include <iostream>

int main(int argc, char** argv)
{
	FastaReader fReader(argv[1]);
	VertexIndex& vertexIndex = VertexIndex::getInstance();
	vertexIndex.setKmerSize(16);

	{
		std::vector<FastaRecord> fastaRecords;
		fReader.GetSequences(fastaRecords);
		for (auto rec : fastaRecords)
		{
			vertexIndex.addFastaSequence(rec);
		}
	}

	vertexIndex.applyKmerThresholds(8, 20);
	vertexIndex.outputCounts();
	vertexIndex.buildReadIndex();

	std::cout << vertexIndex.getIndexByRead().at(100)[0].position << std::endl;
	return 0;
}
