//(c) 2016 by Authors
//This file is a part of ABruijn program.
//Released under the BSD license (see LICENSE file)

#pragma once

#include <unordered_set>
#include <mutex>

#include <cuckoohash_map.hh>

#include "vertex_index.h"
#include "sequence_container.h"
#include "logger.h"


struct OverlapRange
{
	OverlapRange(FastaRecord::Id curId, FastaRecord::Id extId, 
				 int32_t curInit, int32_t extInit): 
		curId(curId), curBegin(curInit), curEnd(curInit), 
		extId(extId), extBegin(extInit), extEnd(extInit){}
	int32_t curRange() const {return curEnd - curBegin;}
	int32_t extRange() const {return extEnd - extBegin;}

	OverlapRange reverse() const
	{
		OverlapRange rev(*this);
		std::swap(rev.curId, rev.extId);
		std::swap(rev.curBegin, rev.extBegin);
		std::swap(rev.curEnd, rev.extEnd);
		rev.leftShift = -rev.leftShift;
		rev.rightShift = -rev.rightShift;
		return rev;
	}

	OverlapRange complement() const
	{
		int32_t curLen = SequenceContainer::get().seqLen(curId);
		int32_t extLen = SequenceContainer::get().seqLen(extId);

		OverlapRange comp(*this);
		std::swap(comp.leftShift, comp.rightShift);
		comp.leftShift = -comp.leftShift;
		comp.rightShift = -comp.rightShift;

		std::swap(comp.curBegin, comp.curEnd);
		comp.curBegin = curLen - comp.curBegin;
		comp.curEnd = curLen - comp.curEnd;

		std::swap(comp.extBegin, comp.extEnd);
		comp.extBegin = extLen - comp.extBegin;
		comp.extEnd = extLen - comp.extEnd;

		comp.curId = comp.curId.rc();
		comp.extId = comp.extId.rc();

		return comp;
	}

	bool contains(int32_t curPos, int32_t extPos) const
	{
		return curBegin <= curPos && curPos <= curEnd &&
			   extBegin <= extPos && extPos <= extEnd;
	}

	//current read
	FastaRecord::Id curId;
	int32_t curBegin;
	int32_t curEnd;
	int32_t leftShift;
	//extension read
	FastaRecord::Id extId;
	int32_t extBegin;
	int32_t extEnd;
	int32_t rightShift;
};

class OverlapDetector
{
public:
	OverlapDetector(int maximumJump, int minimumOverlap,
					int maximumOverhang, int coverage):
		_maximumJump(maximumJump), 
		_minimumOverlap(minimumOverlap),
		_maximumOverhang(maximumOverhang), 
		_coverage(coverage),
		_vertexIndex(VertexIndex::get()),
		_seqContainer(SequenceContainer::get()), 
		_progress(_seqContainer.getIndex().size()),
		_nextJob(0)
	{}

	typedef std::unordered_map<FastaRecord::Id, 
					   std::vector<OverlapRange>> OverlapIndex;
	
	void findAllOverlaps(size_t numThreads);
	const OverlapIndex& getOverlapIndex() const {return _overlapIndex;}
private:
	enum JumpRes {J_END, J_INCONS, J_CLOSE, J_FAR};

	std::vector<OverlapRange> getReadOverlaps(FastaRecord::Id readId) const;
	void 	addOverlapShifts(OverlapRange& ovlp) const;
	bool    goodStart(int32_t currentPos, int32_t extensionPos, 
				      FastaRecord::Id currentId, 
					  FastaRecord::Id extensionId) const;
	bool    overlapTest(const OverlapRange& ovlp, 
						int32_t curLen, int32_t extLen) const;
	JumpRes jumpTest(int32_t currentPrev, int32_t currentNext,
				     int32_t extensionPrev, int32_t extensionNext) const;
	void 	parallelWorker();

	typedef std::tuple<FastaRecord::Id, FastaRecord::Id> id_pair_t;
	struct key_hash : public std::unary_function<id_pair_t, std::size_t>
	{
		 size_t operator()(const id_pair_t& k) const
		 {
			size_t lhs = std::get<0>(k).hash();
			size_t rhs = std::get<1>(k).hash();
			lhs ^= rhs + 0x9ddfea08eb382d69ULL + (lhs << 6) + (lhs >> 2);
			return lhs;
		 }
	};

	const int _maximumJump;
	const int _minimumOverlap;
	const int _maximumOverhang;
	const int _coverage;

	OverlapIndex _overlapIndex;
	cuckoohash_map<id_pair_t, bool, key_hash> _overlapMatrix;

	const VertexIndex& _vertexIndex;
	const SequenceContainer& _seqContainer;

	std::mutex _fetchMutex;
	mutable std::mutex _logMutex;
	ProgressPercent _progress;
	std::vector<FastaRecord::Id> _jobQueue;
	size_t _nextJob;
};
