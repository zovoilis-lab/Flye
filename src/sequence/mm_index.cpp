#include "mm_index.h"

#include <cstdlib>

namespace
{
	int window_size = 5;
	int kmer_size = 19;
	bool is_hpc = true;
	int bucket_bits = 14;
}

MinimapIndex::MinimapIndex(const SequenceContainer &readsContainer): 
	_numOfSequences(readsContainer.getIndex().size() / 2),
	_pSequences(new const char*[_numOfSequences]),
	_pSequencesIds(new const char*[_numOfSequences]),
	_minimapIndex(nullptr)
{
	std::cout << "In MinimapIndex constructor" << std::endl;
	size_t i = 0;
	for (auto &hashPair : readsContainer.getIndex())
	{
		if (hashPair.first.strand())
		{
			_sequences.push_back(hashPair.second.sequence.str());
			_sequencesIds.push_back(hashPair.first.toString());
			_pSequences[i] = _sequences[i].data();
			_pSequencesIds[i] = _sequencesIds[i].data();
			i += 1;
		}
	}

	_minimapIndex = mm_idx_str(window_size, kmer_size, is_hpc, bucket_bits,
							   _numOfSequences, _pSequences, _pSequencesIds);
	std::cout << "MinimapIndex has been built!" << std::endl;
	//std::cout << "Sequences: " << std::endl;

	//for (size_t i = 0; i < _numOfSequences; ++i)
	//{
	//    std::cout << "id" << _sequencesIds[i] << ' ' << _sequences[i] << std::endl;
	//}
}

int32_t MinimapIndex::getSeqId(int32_t id) const
{
	return atoi(_minimapIndex->seq[id].name);
}

int32_t MinimapIndex::getSeqLength(int32_t id) const
{
	return _minimapIndex->seq[id].len;
}

mm_idx_t* MinimapIndex::get() const
{
	return _minimapIndex;
}

void MinimapIndex::clear()
{
	_sequences.clear();
	_sequencesIds.clear();
	delete [] _pSequences;
	delete [] _pSequencesIds;
	mm_idx_destroy(_minimapIndex);
}

MinimapIndex::~MinimapIndex()
{
	std::cout << "In MinimapIndex destructor" << std::endl;
	clear();
}