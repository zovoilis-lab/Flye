//(c) 2016 by Authors
//This file is a part of ABruijn program.
//Released under the BSD license (see LICENSE file)

#pragma once

#include <unordered_set>
#include <mutex>
#include <sstream>

#include <cuckoohash_map.hh>

#include "vertex_index.h"
#include "sequence_container.h"
#include "logger.h"


struct OverlapRange
{
	OverlapRange(FastaRecord::Id curId = FastaRecord::ID_NONE, 
				 FastaRecord::Id extId = FastaRecord::ID_NONE, 
				 int32_t curInit = 0, int32_t extInit = 0): 
		curId(curId), curBegin(curInit), curEnd(curInit), 
		extId(extId), extBegin(extInit), extEnd(extInit)
	{}
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

	std::string serialize() const
	{
		std::stringstream ss;
		ss << curId << " " << curBegin << " " << curEnd << " " 
		   << leftShift << " " << extId << " " << extBegin << " " 
		   << extEnd << " " << rightShift;
		return ss.str();
	}

	void unserialize(const std::string& str)
	{
		std::stringstream ss(str);
		ss >> curId >> curBegin >> curEnd >> leftShift 
		   >> extId >> extBegin >> extEnd >> rightShift;
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
	OverlapDetector(int coverage):
		_coverage(coverage),
		_vertexIndex(VertexIndex::get()),
		_seqContainer(SequenceContainer::get()), 
		_progress(_seqContainer.getIndex().size()),
		_nextJob(0)
	{}

	typedef std::unordered_map<FastaRecord::Id, 
					   std::vector<OverlapRange>> OverlapIndex;
	
	void findAllOverlaps(size_t numThreads);
	void saveOverlaps(const std::string filename);
	void loadOverlaps(const std::string filename);
	const OverlapIndex& getOverlapIndex() const {return _overlapIndex;}

private:
	enum JumpRes {J_END, J_INCONS, J_CLOSE, J_FAR};

	std::vector<OverlapRange> getReadOverlaps(FastaRecord::Id readId) const;
	void 	addOverlapShifts(OverlapRange& ovlp,
							 const std::vector<KmerPosition>& 
							 	solidKmersCache,
							 int32_t curLen, int32_t extLen) const;
	bool    goodStart(int32_t currentPos, int32_t extensionPos, 
				      int32_t curLen, int32_t extLen) const;
	bool    overlapTest(const OverlapRange& ovlp, 
						int32_t curLen, int32_t extLen) const;
	JumpRes jumpTest(int32_t currentPrev, int32_t currentNext,
				     int32_t extensionPrev, int32_t extensionNext) const;
	void 	parallelWorker();

	const int _coverage;

	OverlapIndex _overlapIndex;
	cuckoohash_map<FastaRecord::IdPair, bool> _overlapMatrix;

	const VertexIndex& _vertexIndex;
	const SequenceContainer& _seqContainer;

	std::mutex _fetchMutex;
	ProgressPercent _progress;
	std::vector<FastaRecord::Id> _jobQueue;
	size_t _nextJob;
	//mutable std::mutex _logMutex;
};
