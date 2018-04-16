#include "mm_index.h"

#include <cstddef>

namespace
{
    int windows_size = 5;
    int kmer_size = 19;
    bool is_hpc = true;
    int bucket_bits = 14;
}

#include <iostream>

MinimapIndex::MinimapIndex(const SequenceContainer &readsContainer)
        : _numOfSequences(readsContainer.getIndex().size() / 2)
        , _pSequences(new const char*[_numOfSequences])
        , _pSequencesIds(new const char*[_numOfSequences])
        , _minimapIndex(nullptr)
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

    std::cout << "_numOfSeq: " << _numOfSequences << std::endl;
    std::cout << "i: " << i << std::endl;

    _minimapIndex = mm_idx_str(windows_size, kmer_size, is_hpc, bucket_bits,
                               _numOfSequences, _pSequences, _pSequencesIds);

    std::cout << "MinimapIndex has been built!" << std::endl;
}

mm_idx_t* MinimapIndex::get() const
{
    return _minimapIndex;
}

int32_t MinimapIndex::getSequenceId(size_t index) const
{
    return atoi(_minimapIndex->seq[index].name);
}

int32_t MinimapIndex::getSequenceLen(size_t index) const
{
    return _minimapIndex->seq[index].len;
}

MinimapIndex::~MinimapIndex()
{
    std::cout << "In MinimapIndex destructor" << std::endl;
    _sequences.clear();
    _sequencesIds.clear();
    delete [] _pSequences;
    delete [] _pSequencesIds;
    mm_idx_destroy(_minimapIndex);
}