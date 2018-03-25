#include "mm_index.h"


MinimapIndex::MinimapIndex(const SequenceContainer &readsContainer): 
    _numOfSequences(readsContainer.getIndex().size()),
    _pSequences(new const char*[_numOfSequences]),
    _minimapIndex(nullptr)
{
    std::cout << "In MinimapIndex constructor" << std::endl;
    size_t i = 0;
    for (auto &hashPair : readsContainer.getIndex())
    {
        _sequences.push_back(hashPair.second.sequence.str());
        _pSequences[i] = _sequences[i].data();
        i += 1;
    }

    int window_size = 5;
    int kmer_size = 19;
    bool is_hpc = true;
    int bucket_bits = 14;

    _minimapIndex = mm_idx_str(window_size, kmer_size, is_hpc, bucket_bits, 
                               _numOfSequences, _pSequences, nullptr);
    std::cout << "MinimapIndex has been built!" << std::endl;
}

void MinimapIndex::clear()
{
    _sequences.clear();
    delete [] _pSequences;
    mm_idx_destroy(_minimapIndex);
}

MinimapIndex::~MinimapIndex()
{
    std::cout << "In MinimapIndex destructor" << std::endl;
    clear();
}