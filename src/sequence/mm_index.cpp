#include "mm_index.h"

#include <cstddef>

#include <iostream>
#include <fstream>

///
#include <zlib.h>
#include "minimap.h"
///

MinimapIndex::MinimapIndex(const SequenceContainer &readsContainer, 
						   const std::string &presetOptions, bool align,
						   bool onlyMax, bool fromFile)
        : _numOfSequences(readsContainer.getIndex().size() / 2)
        , _indexForStringsInMemory(nullptr)
        , _mapOptions(new mm_mapopt_t())
        , _indexOptions(new mm_idxopt_t())
{
    std::cout << "In MinimapIndex constructor" << std::endl;
    mm_set_opt(0, _indexOptions, _mapOptions);

    if (align)
    {
        _mapOptions->flag |= MM_F_CIGAR;
    }
    else
    {
        _indexOptions->flag |= MM_I_NO_SEQ;
    }
    if (!onlyMax) _mapOptions->pri_ratio = 0;
    //_indexOptions->batch_size = 4000000000ULL; // 4.0Gb
    _indexOptions->batch_size = 120000000ULL;
    mm_set_opt(presetOptions.c_str(), _indexOptions, _mapOptions);
    mm_check_opt(_indexOptions, _mapOptions);

    if (fromFile)
    {
        std::cout << "building an index from file" << std::endl;
        builtFromFile = true;
        int n_threads = 4;

        //gzFile f = gzopen("__reads__.fasta", "r");

        //(void)kseq_read;
        //kseq_t *ks = kseq_init(f);
        mm_idx_reader_t *r = mm_idx_reader_open("input_reads.fasta", _indexOptions, 0);
        mm_idx_t *idx;

        int nIndexParts = 0;
        while ((idx = mm_idx_reader_read(r, n_threads)) != 0)
        {
            nIndexParts += 1;
            _indexes.push_back(idx);
            std::cout << "nIndexParts = " << nIndexParts << std::endl;
        }

        for (size_t i = 0; i < _indexes.size(); ++i)
        {
            std::cout << _indexes[i] << std::endl;
            mm_idx_stat(_indexes[i]);
        }

        mm_idx_reader_close(r);

        std::cout << "Done!" << std::endl;
    }
    else
    {
        std::cout << "building an index from sequences in memory" << std::endl;
        builtFromFile = false;
        std::string description;
        size_t descriptionLength;

        for (auto &hashPair : readsContainer.getIndex())
        {
            if (hashPair.first.strand())
            {
                descriptionLength = hashPair.second.description.length();
                description = hashPair.second.description.substr(1, descriptionLength - 1);
                _sequences.push_back(hashPair.second.sequence.str());
                _sequencesDescriptions.push_back(description);
            }
        }

        for (size_t i = 0; i < _numOfSequences; ++i)
        {
            _pSequences.push_back(_sequences[i].data());
            _pSequencesDescriptions.push_back(_sequencesDescriptions[i].data());
        }

        _indexForStringsInMemory = mm_idx_str(_indexOptions->w, _indexOptions->k, /*is hpc*/ false,
                                              /*bucket bits*/ 14, _numOfSequences, _pSequences.data(),
                                              _pSequencesDescriptions.data());

        mm_idx_stat(_indexForStringsInMemory);
        mm_mapopt_update(_mapOptions, _indexForStringsInMemory);

        _sequences.clear();

        std::cout << "Done!" << std::endl;
    }
}

/*
MinimapIndex::MinimapIndex(const SequenceContainer &readsContainer, const std::string &presetOptions)
        : _numOfSequences(readsContainer.getIndex().size() / 2)
        , _minimapIndex(nullptr)
        , _mapOptions(new mm_mapopt_t())
        , _idxOptions(new mm_idxopt_t())
{
    //std::cout << "In MinimapIndex constructor" << std::endl;

    //mm_mapopt_init(_minimapOptions);
	//mm_idxopt_t ipt;
	//mm_mapopt_t opt;

    std::cout << "preset options: " <<  presetOptions << std::endl;

    if (presetOptions == "ava-pb")
    {
		mm_set_opt("ava-pb", &ipt, _minimapOptions);
        _minimapOptions->flag |= MM_F_ALL_CHAINS | MM_F_NO_DIAG | MM_F_NO_DUAL | MM_F_NO_LJOIN;
        _minimapOptions->min_chain_score = 100;
        _minimapOptions->pri_ratio = 0.0f;
        _minimapOptions->max_gap = 10000;
        _minimapOptions->max_chain_skip = 25;


        std::cout << presetOptions << " is used" << std::endl;
    }

    if (presetOptions == "asm5")
    {
        // io->flag = 0; ? this option is not used by mm_idx_str
        // io->k = 19, io->w = 19; kmer size is the same as for ava-pb preset, window size is 19 (instead of 5)

        windows_size = 19;
        is_hpc = false; // this option is not used for asm5 preset, so it is false
		mm_set_opt("map-pb", &ipt, _minimapOptions);
		_minimapOptions->pri_ratio = 0;
		_minimapOptions->flag |= MM_F_CIGAR;


        _minimapOptions->pri_ratio = 0;
        _minimapOptions->a = 1;
        _minimapOptions->b = 19;
        _minimapOptions->q = 39;
        _minimapOptions->q2 = 81;
        _minimapOptions->e = 3;
        _minimapOptions->e2 = 1;
        _minimapOptions->zdrop = _minimapOptions->zdrop_inv = 200;
        _minimapOptions->min_dp_max = 200;
        _minimapOptions->best_n = 50;
        std::cout << presetOptions << " is used" << std::endl;
    }

    if (presetOptions == "asm10")
    {
        // io->flag = 0, io->k = 19, io->w = 19;
        windows_size = 19;
        kmer_size = 19;
		mm_set_opt("map-pb", &ipt, _minimapOptions);
		_minimapOptions->pri_ratio = 0;
		_minimapOptions->flag |= MM_F_CIGAR;

        is_hpc = false;
        _minimapOptions->flag = 0; // ?
        _minimapOptions->a = 1;
        _minimapOptions->b = 9;
        _minimapOptions->q = 16;
        _minimapOptions->q2 = 41;
        _minimapOptions->e = 2;
        _minimapOptions->e2 = 1;
        _minimapOptions->zdrop = _minimapOptions->zdrop_inv = 200;
        _minimapOptions->min_dp_max = 200;
        _minimapOptions->best_n = 50;

        std::cout << presetOptions << " is used" << std::endl;
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
    std::cout << "Clear sequences..." << std::endl;
    _sequences.clear();
}
*/

mm_idx_t* MinimapIndex::get(size_t index) const
{
    if (builtFromFile)
    {
        return _indexes[index];
    }
    else
    {
        return _indexForStringsInMemory;
    }
}

std::string MinimapIndex::getSequenceDescription(size_t indexNum, size_t index) const
{
    if (builtFromFile)
    {
        return std::string(_indexes[indexNum]->seq[index].name);
    }
    else
    {
        return std::string(_indexForStringsInMemory->seq[index].name);
    }
}

int32_t MinimapIndex::getSequenceLen(size_t indexNum, size_t index) const
{
    if (builtFromFile)
    {
        return _indexes[indexNum]->seq[index].len;
    }
    else
    {
        return _indexForStringsInMemory->seq[index].len;
    }
}

size_t MinimapIndex::getNumOfIndexes() const
{
    if (builtFromFile)
    {
        return _indexes.size();
    }
    else
    {
        return 1;
    }
}

mm_mapopt_t* MinimapIndex::getOptions() const
{
    return _mapOptions;
}

MinimapIndex::~MinimapIndex()
{
    std::cout << "In MinimapIndex destructor" << std::endl;
    _sequences.clear();
    _sequencesDescriptions.clear();
    _pSequences.clear();
    _pSequencesDescriptions.clear();

    for (size_t i = 0; i < _indexes.size(); ++i)
    {
        mm_idx_destroy(_indexes[i]);
    }
    mm_idx_destroy(_indexForStringsInMemory);

    delete _mapOptions;
    delete _indexOptions;
}
