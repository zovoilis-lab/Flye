
#include "kmer_index.h"
#include "fasta.h"

int main(int argc, char** argv)
{
	FastaReader fReader(argv[1]);
	std::vector<FastaRecord> fastaRecords;
	fReader.GetSequences(fastaRecords);

	KmerIndex& kmerIndex = KmerIndex::getIndex();
	kmerIndex.setKmerSize(16);
	for (auto rec : fastaRecords)
	{
		kmerIndex.addFastaSequence(rec);
	}
	kmerIndex.outputCounts();
	return 0;
}
