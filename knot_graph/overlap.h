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
#include "progress_bar.h"


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

	OverlapRange complement(int32_t curLen, int32_t extLen) const
	{
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

	int32_t curIntersect(const OverlapRange& other) const
	{
		return std::min(curEnd, other.curEnd) - 
			   std::max(curBegin, other.curBegin);
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
	OverlapDetector(const SequenceContainer& seqContainer,
					const VertexIndex& vertexIndex,
					int maxJump, int minOverlap, int maxOverhang):
		_maxJump(maxJump),
		_minOverlap(minOverlap),
		_maxOverhang(maxOverhang),
		_checkOverhang(maxOverhang > 0),
		_vertexIndex(vertexIndex),
		_seqContainer(seqContainer)
	{}

	std::vector<OverlapRange> 
	getSeqOverlaps(const std::string& sequence,
				   FastaRecord::Id idTrivial = FastaRecord::ID_NONE) const;

private:
	enum JumpRes {J_END, J_INCONS, J_CLOSE, J_FAR};

	
	bool    goodStart(int32_t currentPos, int32_t extensionPos, 
				      int32_t curLen, int32_t extLen) const;
	bool    overlapTest(const OverlapRange& ovlp, 
						int32_t curLen, int32_t extLen) const;
	JumpRes jumpTest(int32_t currentPrev, int32_t currentNext,
				     int32_t extensionPrev, int32_t extensionNext) const;

	const int _maxJump;
	const int _minOverlap;
	const int _maxOverhang;
	const bool _checkOverhang;

	const VertexIndex& _vertexIndex;
	const SequenceContainer& _seqContainer;
};


class OverlapContainer
{
public:
	OverlapContainer(const OverlapDetector& ovlpDetect,
					 const SequenceContainer& queryContainer):
		_ovlpDetect(ovlpDetect),
		_queryContainer(queryContainer)
	{}

	typedef std::unordered_map<FastaRecord::Id, 
					   std::vector<OverlapRange>> OverlapIndex;

	void saveOverlaps(const std::string& filename);
	void loadOverlaps(const std::string& filename);

	void findAllOverlaps();
	const OverlapIndex& getOverlapIndex() const {return _overlapIndex;}

private:
	//void 	addOverlapShifts(OverlapRange& ovlp,
	//						 const std::vector<KmerPosition>& 
	//						 	solidKmersCache,
	//						 int32_t curLen, int32_t extLen) const;

	const OverlapDetector& _ovlpDetect;
	const SequenceContainer& _queryContainer;

	OverlapIndex _overlapIndex;
	//cuckoohash_map<FastaRecord::IdPair, bool> _overlapMatrix;
};
