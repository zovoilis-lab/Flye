#pragma once

#include "../../lib/minimap2/minimap.h"
#include "sequence_container.h"

#include <vector>
#include <string>

class MinimapIndex
{
public:
    explicit MinimapIndex(const SequenceContainer&, const std::string &presetOptions);
    ~MinimapIndex();

    MinimapIndex(const MinimapIndex&) = delete;
    void operator=(const MinimapIndex&) = delete;

    mm_idx_t* get() const;

    int32_t getSequenceId(size_t index) const;
    int32_t getSequenceLen(size_t index) const;

    const char* getSequence(size_t index) const;

    mm_mapopt_t* getOptions() const;

private:
    void saveSequencesToFile(const SequenceContainer&);

    size_t _numOfSequences;
    std::vector<const char*> _pSequences;
    std::vector<const char*> _pSequencesIds;
    std::vector<std::string> _sequences;
    std::vector<std::string> _sequencesIds;

    mm_idx_t* _minimapIndex;
    mm_mapopt_t* _mapOptions;
    mm_idxopt_t* _indexOptions;
};