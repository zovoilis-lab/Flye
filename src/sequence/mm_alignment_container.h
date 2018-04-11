#pragma once

#include <vector>
#include "overlap.h"
#include "mm_index.h"

class MinimapAlignmentContainer
{
public:
	MinimapAlignmentContainer(mm_reg1_t*, int);
	~MinimapAlignmentContainer();

	MinimapAlignmentContainer(const MinimapAlignmentContainer&) = delete;
	void operator=(const MinimapAlignmentContainer&) = delete;

	void clear();

	int getNumOfAlignments() const;

	int32_t getCurStart(int index) const;
	int32_t getCurEnd(int index) const;
	bool curStrand(int index) const; // true -- positive, false -- negative

	int32_t getExtIndexId(int index) const;
	int32_t getExtStart(int index) const;
	int32_t getExtEnd(int index) const;

	int32_t getScore(int index) const;
private:
	mm_reg1_t* _pAlignments;
	int32_t    _numOfAlignments;
};