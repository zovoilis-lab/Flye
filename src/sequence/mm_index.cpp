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
        , _minimapIndex(nullptr)
{
    std::cout << "In MinimapIndex constructor" << std::endl;

    size_t total_length = 0;
    for (auto &hashPair : readsContainer.getIndex())
    {
        if (hashPair.first.strand())
        {
            _sequences.push_back(hashPair.second.sequence.str());
            _sequencesIds.push_back(hashPair.first.toString());
            total_length += _sequences[_sequences.size() - 1].length();
        }
    }

    for (size_t i = 0; i < _numOfSequences; ++i)
    {
        _pSequences.push_back(_sequences[i].data());
        _pSequencesIds.push_back(_sequencesIds[i].data());
    }

    std::cout << "total length = " << total_length << std::endl;
    std::cout << "_numOfSeq: " << _numOfSequences << std::endl;

    _minimapIndex = mm_idx_str(windows_size, kmer_size, is_hpc, bucket_bits,
                               _numOfSequences, _pSequences.data(), _pSequencesIds.data());

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
    _pSequences.clear();
    _pSequencesIds.clear();
    mm_idx_destroy(_minimapIndex);
}
