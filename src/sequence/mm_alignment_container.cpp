#include "mm_alignment_container.h"

#include <unordered_map>
#include <cstring>
#include <cstdlib>

MinimapAlignmentContainer
::MinimapAlignmentContainer(mm_reg1_t *pAlignments, int numOfAlignments)
        : _pAlignments(pAlignments)
        , _numOfAlignments(numOfAlignments)
{}

void MinimapAlignmentContainer::clear()
{
    if (_pAlignments)
    {
        for (int i = 0; i < _numOfAlignments; ++i)
        {
            free(_pAlignments[i].p);
        }
        free(_pAlignments);
    }

    _pAlignments = nullptr;
    _numOfAlignments = 0;
}

int MinimapAlignmentContainer::getNumOfAlignments() const
{
    return _numOfAlignments;
}

int32_t MinimapAlignmentContainer::getCurStart(int index) const
{
    return _pAlignments[index].qs;
}

int32_t MinimapAlignmentContainer::getCurEnd(int index) const
{
    return _pAlignments[index].qe;
}

bool MinimapAlignmentContainer::curStrand(int index) const
{
    return !_pAlignments[index].rev;
}

int32_t MinimapAlignmentContainer::getExtIndexId(int index) const
{
    return _pAlignments[index].rid;
}

int32_t MinimapAlignmentContainer::getExtStart(int index) const
{
    return _pAlignments[index].rs;
}

int32_t MinimapAlignmentContainer::getExtEnd(int index) const
{
    return _pAlignments[index].re;
}

int32_t MinimapAlignmentContainer::getScore(int index) const
{
    return _pAlignments[index].score;
}

MinimapAlignmentContainer::~MinimapAlignmentContainer()
{
    clear();
}

/*
std::vector<OverlapRange>
MinimapAlignmentContainer::
getBestOverlaps(int32_t curId, int32_t curLen, const MinimapIndex &minimapIndex)
{
    std::unordered_map<int32_t, OverlapRange> bestOverlapsHash;
    std::vector<OverlapRange> bestOverlapsVector;

    for (int i = 0; i < _numOfAlignments; ++i)
    {
        int32_t extId = minimapIndex.getSeqId(_pAlignments[i].rid);
        if (curId != extId)
        {
            int32_t curStart = _pAlignments[i].qs;
            int32_t curEnd = _pAlignments[i].qe;
            bool reverse = _pAlignments[i].rev;

            int32_t extStart = _pAlignments[i].rs;
            int32_t extEnd = _pAlignments[i].re;
            int32_t extLen = minimapIndex.getSeqLength(_pAlignments[i].rid);

            int32_t score = _pAlignments[i].mlen;
        }
    }

    return bestOverlapsVector;
}
*/