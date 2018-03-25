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

    ~MinimapIndex();

private:
    const size_t _numOfSequences;
    const char **_pSequences;
    std::vector<std::string> _sequences;
    mm_idx_t *_minimapIndex;
};