#pragma once

#include "sequence_container.h"
#include "../../lib/minimap2/minimap.h"

#include <vector>
#include <string>


class MinimapIndex
{
public:
    explicit MinimapIndex(const SequenceContainer &);

    MinimapIndex(const MinimapIndex&) = delete;
    void operator=(const MinimapIndex&) = delete;

    void clear();
    mm_idx_t* get() const;

    int32_t getSeqId(int32_t id) const;
    int32_t getSeqLength(int32_t id) const;

    ~MinimapIndex();

private:
    const size_t _numOfSequences;
    const char **_pSequences;
    const char **_pSequencesIds;
    std::vector<std::string> _sequences;
    std::vector<std::string> _sequencesIds;
    mm_idx_t *_minimapIndex;
};