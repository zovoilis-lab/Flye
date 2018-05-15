#pragma once

#include "../../lib/minimap2/minimap.h"
#include "sequence_container.h"

#include <vector>
#include <string>

class MinimapIndex
{
public:
    MinimapIndex(const SequenceContainer&,
                 const std::string &presetOptions,
                 bool align, bool onlyMax, bool fromFile);
    ~MinimapIndex();

    MinimapIndex(const MinimapIndex&) = delete;
    void operator=(const MinimapIndex&) = delete;

    mm_idx_t* get(size_t) const;

    std::string getSequenceDescription(size_t, size_t) const;
    int32_t getSequenceLen(size_t, size_t) const;
	size_t getNumOfIndexes() const;

    mm_mapopt_t* getOptions() const;

private:
    size_t _numOfSequences;
    std::vector<const char*> _pSequences;
    std::vector<const char*> _pSequencesDescriptions;
    std::vector<std::string> _sequences;
    std::vector<std::string> _sequencesDescriptions;

    mm_idx_t* _indexForStringsInMemory;
    mm_mapopt_t* _mapOptions;
    mm_idxopt_t* _indexOptions;
	std::vector<mm_idx_t*> _indexes;
    bool builtFromFile;
};
