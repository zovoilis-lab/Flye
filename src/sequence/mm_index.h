#pragma once

#include "../../lib/minimap2/minimap.h"
#include "sequence_container.h"

#include <vector>
#include <string>

class MinimapIndex
{
public:
    explicit MinimapIndex(const SequenceContainer&);
    ~MinimapIndex();

    MinimapIndex(const MinimapIndex&) = delete;
    void operator=(const MinimapIndex&) = delete;

    mm_idx_t* get() const;

    int32_t getSequenceId(size_t index) const;
    int32_t getSequenceLen(size_t index) const;

    const char* getSequence(size_t index) const;

private:
    size_t _numOfSequences;
    const char** _pSequences;
    const char** _pSequencesIds;

    std::vector<std::string> _sequences;
    std::vector<std::string> _sequencesIds;

    mm_idx_t* _minimapIndex;
};