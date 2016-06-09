//(c) 2016 by Authors
//This file is a part of ABruijn program.
//Released under the BSD license (see LICENSE file)

#pragma once

#include <unordered_set>
#include <mutex>

#include "kmer_index.h"
#include "fasta.h"
#include "logger.h"


struct OverlapRange
{
	OverlapRange(FastaRecord::Id curId, FastaRecord::Id extId, 
				 int32_t curInit, int32_t extInit): 
		curId(curId), curBegin(curInit), curEnd(curInit), 
		extId(extId), extBegin(extInit), extEnd(extInit){}
	int32_t curRange() const {return curEnd - curBegin;}
	int32_t extRange() const {return extEnd - extBegin;}

	void reverse()
	{
		std::swap(curId, extId);
		std::swap(curBegin, extBegin);
		std::swap(curEnd, extEnd);
		leftShift = -leftShift;
		rightShift = -rightShift;
	}

	void complement(int32_t curLen, int32_t extLen)
	{
		std::swap(leftShift, rightShift);
		leftShift = -leftShift;
		rightShift = -rightShift;

		std::swap(curBegin, curEnd);
		curBegin = curLen - curBegin;
		curEnd = curLen - curEnd;

		std::swap(extBegin, extEnd);
		extBegin = extLen - extBegin;
		extEnd = extLen - extEnd;

		curId = curId.rc();
		extId = extId.rc();
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
					int maximumOverhang, const VertexIndex& vi,
					const SequenceContainer& seqCont):
		_maximumJump(maximumJump), 
		_minimumOverlap(minimumOverlap),
		_maximumOverhang(maximumOverhang), 
		_vertexIndex(vi),
		_seqContainer(seqCont), 
		_progress(seqCont.getIndex().size()),
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
	//bool    goodStart(int32_t currentPos, int32_t extensionPos, 
	//			      int32_t currentLength, int32_t extensionLength) const;
	bool    overlapTest(const OverlapRange& ovlp, int32_t curLen, int32_t extLen) const;
	JumpRes jumpTest(int32_t currentPrev, int32_t currentNext,
				     int32_t extensionPrev, int32_t extensionNext) const;
	void 	parallelWorker();

	typedef std::tuple<FastaRecord::Id, FastaRecord::Id> id_pair_t;
	struct key_hash : public std::unary_function<id_pair_t, std::size_t>
	{
		 size_t operator()(const id_pair_t& k) const
		 {
		   	return ((size_t)std::get<0>(k).hash() << 32) ^ 
					std::get<1>(k).hash();
		 }
	};

	const int _maximumJump;
	const int _minimumOverlap;
	const int _maximumOverhang;

	OverlapIndex _overlapIndex;
	std::unordered_set<id_pair_t, key_hash> _overlapMatrix;

	const VertexIndex& _vertexIndex;
	const SequenceContainer& _seqContainer;

	std::mutex _fetchMutex;
	mutable std::mutex _logMutex;
	ProgressPercent _progress;
	std::vector<FastaRecord::Id> _jobQueue;
	size_t _nextJob;
};
