#pragma once

#include "../../lib/minimap2/minimap.h"
#include "mm_index.h"
#include "overlap.h"

#include <string>
#include <vector>

class MinimapAlignmentContainer
{
public:
    MinimapAlignmentContainer(mm_reg1_t*, int);
    ~MinimapAlignmentContainer();

    MinimapAlignmentContainer(const MinimapAlignmentContainer&) = delete;
    void operator=(const MinimapAlignmentContainer&) = delete;

    void printAllOverlaps(uint32_t, const std::string&,
                          const MinimapIndex&,
                          const SequenceContainer&) const;

    int     getNumOfAlignments()       const;

    int32_t getCurBegin (size_t index) const;
    int32_t getCurEnd   (size_t index) const;
    int32_t getExtBegin (size_t index) const;
    int32_t getExtEnd   (size_t index) const;

    int32_t getExtIndexId(size_t index) const;

    bool    curStrand(size_t index)  const;
    int32_t getScore (size_t index)   const;

    void clear();

private:
    mm_reg1_t* _pAlignments;
    int _numOfAlignments;
};