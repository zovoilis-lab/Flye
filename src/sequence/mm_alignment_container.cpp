#include "mm_alignment_container.h"
#include "mm_index.h"

#include <cstdlib>
#include <iostream>

MinimapAlignmentContainer::MinimapAlignmentContainer(mm_reg1_t* pAlignments,
                                                     int numOfAlignments)
        : _pAlignments(pAlignments)
        , _numOfAlignments(numOfAlignments)
{}

void MinimapAlignmentContainer::printAllOverlaps(uint32_t curId,
                                                 const std::string &cur,
                                                 const MinimapIndex &index,
                                                 const SequenceContainer &readsContainer) const
{
    for (size_t i = 0; i < _numOfAlignments; ++i)
    {
        int extId = index.getSequenceId(_pAlignments[i].rid);
        std::cout << "id" << curId << '\t';
        std::cout << cur.length() << '\t';
        std::cout << _pAlignments[i].qs << '\t';
        std::cout << _pAlignments[i].qe << '\t';
        std::cout << "+-"[_pAlignments[i].rev] << '\t';
        std::cout << "id" << extId << '\t';
        std::cout << index.getSequenceLen(_pAlignments[i].rid) << '\t';
        std::cout << _pAlignments[i].rs << '\t';
        std::cout << _pAlignments[i].re << '\t';
        std::cout << _pAlignments[i].mlen << std::endl;

        for (size_t j = _pAlignments[i].qs; j != _pAlignments[i].qs + 30; ++j)
        {
            std::cout << cur[j];
        }
        std::cout << "..........";

        for (size_t j = _pAlignments[i].qe - 30; j != _pAlignments[i].qe; ++j)
        {
            std::cout << cur[j];
        }
        std::cout << std::endl;

        auto ext = readsContainer.getSeq(FastaRecord::Id(extId));

        for (size_t j = _pAlignments[i].rs; j != _pAlignments[i].rs + 30; ++j)
        {
            std::cout << ext.at(j);
        }
        std::cout << "..........";

        for (size_t j = _pAlignments[i].re - 30; j != _pAlignments[i].re; ++j)
        {
            std::cout << ext.at(j);
        }
        std::cout << std::endl;
    }
}

int MinimapAlignmentContainer::getNumOfAlignments() const
{
    return _numOfAlignments;
}

int32_t MinimapAlignmentContainer::getCurBegin(size_t index) const
{
    return _pAlignments[index].qs;
}

int32_t MinimapAlignmentContainer::getCurEnd(size_t index) const
{
    return _pAlignments[index].qe;
}

int32_t MinimapAlignmentContainer::getExtBegin(size_t index) const
{
    return _pAlignments[index].rs;
}

int32_t MinimapAlignmentContainer::getExtEnd(size_t index) const
{
    return _pAlignments[index].re;
}

bool MinimapAlignmentContainer::curStrand(size_t index) const
{
    return !_pAlignments[index].rev;
}

int32_t MinimapAlignmentContainer::getScore(size_t index) const
{
    return _pAlignments[index].qe - _pAlignments[index].qs;
}

int32_t MinimapAlignmentContainer::getExtIndexId(size_t index) const
{
    return _pAlignments[index].rid;
}

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

MinimapAlignmentContainer::~MinimapAlignmentContainer()
{
    clear();
}