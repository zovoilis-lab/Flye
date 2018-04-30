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

MinimapIndex::MinimapIndex(const SequenceContainer &readsContainer, const std::string &presetOptions)
        : _numOfSequences(readsContainer.getIndex().size() / 2)
        , _minimapIndex(nullptr)
        , _minimapOptions(new mm_mapopt_t())
{
    std::cout << "In MinimapIndex constructor" << std::endl;

    mm_mapopt_init(_minimapOptions);

    if (presetOptions == "ava-pb")
    {
        _minimapOptions->flag |= MM_F_ALL_CHAINS | MM_F_NO_DIAG | MM_F_NO_DUAL | MM_F_NO_LJOIN;
        _minimapOptions->min_chain_score = 100;
        _minimapOptions->pri_ratio = 0.0f;
        _minimapOptions->max_gap = 10000;
        _minimapOptions->max_chain_skip = 25;
    }

    if (presetOptions == "asm5")
    {
        // io->flag = 0; ? this option is not used by mm_idx_str
        // io->k = 19, io->w = 19; kmer size is the same as for ava-pb preset, window size is 19 (instead of 5)

        windows_size = 19;
        is_hpc = false; // this option is not used for asm5 preset, so it is false

        _minimapOptions->a = 1;
        _minimapOptions->b = 19;
        _minimapOptions->q = 39;
        _minimapOptions->q2 = 81;
        _minimapOptions->e = 3;
        _minimapOptions->e2 = 1;
        _minimapOptions->zdrop = _minimapOptions->zdrop_inv = 200;
        _minimapOptions->min_dp_max = 200;
        _minimapOptions->best_n = 50;
    }

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

    mm_mapopt_update(_minimapOptions, _minimapIndex);

    std::cout << "MinimapIndex has been built!" << std::endl;

    // clear sequences after we build an index
    _sequences.clear();
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

mm_mapopt_t* MinimapIndex::getOptions() const
{
    return _minimapOptions;
}

MinimapIndex::~MinimapIndex()
{
    std::cout << "In MinimapIndex destructor" << std::endl;
    _sequencesIds.clear();
    _pSequences.clear();
    _pSequencesIds.clear();
    mm_idx_destroy(_minimapIndex);
    delete _minimapOptions;
}
