//(c) 2016 by Authors
//This file is a part of ABruijn program.
//Released under the BSD license (see LICENSE file)

#include <set>
#include <iostream>
#include <cassert>
#include <algorithm>
#include <thread>

#include "overlap.h"
#include "../common/config.h"
#include "../common/utils.h"
#include "../common/parallel.h"
#include "../common/disjoint_set.h"

#include "../../lib/minimap2/minimap.h"
#include "mm_index.h"
#include "mm_buffer.h"
#include "mm_alignment_container.h"

#include <cassert>

//reject overlaps early to speed everything up
bool OverlapDetector::goodStart(int32_t curPos, int32_t extPos, 
								int32_t curLen, int32_t extLen,
								FastaRecord::Id curId, 
								FastaRecord::Id extId) const
{	
	if (_checkOverhang)
	{
		//allow overlaps between a read and its complement to
		//bypass overhang checks to detect the typical PacBio chimera pattern
		if (std::min(curPos, extPos) > _maxOverhang &&
			curId != extId.rc()) return false;
	}

	if (extPos > extLen - _minOverlap ||
	    curPos > curLen - _minOverlap) return false;

	return true;
}


OverlapDetector::JumpRes 
OverlapDetector::jumpTest(int32_t curPrev, int32_t curNext,
						  int32_t extPrev, int32_t extNext) const
{
	//static const int CLOSE_JUMP = Config::get("close_jump_rate");
	static const int FAR_JUMP = Config::get("jump_divergence_rate");

	if (curNext - curPrev > _maxJump) return J_END;

	float jumpLength = ((curNext - curPrev) + (extNext - extPrev)) / 2.0f;
	float jumpDiv = abs((curNext - curPrev) - (extNext - extPrev));

	if (0 < curNext - curPrev && curNext - curPrev < _maxJump &&
		0 < extNext - extPrev && extNext - extPrev < _maxJump)
	{
		if (jumpLength < _gapSize &&
			jumpDiv < Parameters::get().kmerSize)
		{
			return J_CLOSE;
		}

		if (jumpDiv < _maxJump / FAR_JUMP)
		{
			return J_FAR;
		}
	}
	return J_INCONS;
}


//Check if it is a proper overlap
bool OverlapDetector::overlapTest(const OverlapRange& ovlp,
								  bool& outSuggestChimeric) const
{
	static const int OVLP_DIVERGENCE = Config::get("overlap_divergence_rate");
	if (ovlp.curRange() < _minOverlap || 
		ovlp.extRange() < _minOverlap) 
	{
		return false;
	}

	float lengthDiff = abs(ovlp.curRange() - ovlp.extRange());
	float meanLength = (ovlp.curRange() + ovlp.extRange()) / 2.0f;
	if (lengthDiff > meanLength / OVLP_DIVERGENCE)
	{
		return false;
	}

	if (ovlp.curId == ovlp.extId.rc()) outSuggestChimeric = true;
	if (_checkOverhang)
	{
		if (std::min(ovlp.curBegin, ovlp.extBegin) > 
			_maxOverhang) 
		{
			return false;
		}
		if (std::min(ovlp.curLen - ovlp.curEnd, ovlp.extLen - ovlp.extEnd) > 
			_maxOverhang)
		{
			return false;
		}
	}

	return true;
}

namespace
{
	struct DPRecord
	{
		DPRecord() {}
		DPRecord(const OverlapRange& ovlp):
			ovlp(ovlp) {}

		OverlapRange ovlp;
		std::vector<int32_t> shifts;
		//bool wasContinued;
	};
}

std::vector<OverlapRange> 
OverlapDetector::getSeqOverlaps(const FastaRecord& fastaRec, 
								bool uniqueExtensions,
								bool& outSuggestChimeric) const
{
    MinimapBuffer buffer;
    int numOfAlignments;

    std::string cur = fastaRec.sequence.str();
    int curLen = fastaRec.sequence.length();

    mm_reg1_t *pAlignments = mm_map(_minimapIndex.get(), curLen, cur.data(),
                                    &numOfAlignments, buffer.get(), _minimapIndex.getOptions(), 0);

    MinimapAlignmentContainer alignmentContainer(pAlignments, numOfAlignments);

    uint32_t curId = fastaRec.id.get();
    uint32_t curRevCompId = fastaRec.id.rc().get();

    std::vector<OverlapRange> overlaps;
    std::unordered_map<uint32_t, OverlapRange> bestOverlapsHash;

    for (int i = 0; i < alignmentContainer.getNumOfAlignments(); ++i)
    {
        uint32_t extId = _minimapIndex.getSequenceId(alignmentContainer.getExtIndexId(i));
		uint32_t extRevComp = FastaRecord::Id(extId).rc().get();

        if (curId == extId || curRevCompId == extId)
        {
            // outSuggestChimeric = true;
            continue;
        }

        int32_t extLen = _minimapIndex.getSequenceLen(alignmentContainer.getExtIndexId(i));
        int32_t curBegin = alignmentContainer.getCurBegin(i);
        int32_t curEnd = alignmentContainer.getCurEnd(i);
        int32_t extBegin = alignmentContainer.getExtBegin(i);
        int32_t extEnd = alignmentContainer.getExtEnd(i);
        int32_t score = alignmentContainer.getScore(i);

        bool curStrand = alignmentContainer.curStrand(i);

        if (!curStrand)
        {
            extId = extRevComp;
            int32_t newExtBegin = extLen - extEnd - 1;
            int32_t newExtEnd = extLen - extBegin - 1;
            extBegin = newExtBegin;
            extEnd = newExtEnd;
        }

        auto overlap = OverlapRange(curId, curBegin, curEnd, curLen,
                                    extId, extBegin, extEnd, extLen,
									curBegin - extBegin,
                                    (extLen - extEnd) - (curLen - curEnd),
                                    score);

        if (overlapTest(overlap, outSuggestChimeric))
        {
            if (uniqueExtensions)
            {
                auto overlapPair = bestOverlapsHash.find(extId);

                if (overlapPair != bestOverlapsHash.end())
                {
                    if (overlapPair->second.score < score)
                    {
                        overlapPair->second = overlap;
                    }
                }
                else
                {
                    bestOverlapsHash[extId] = overlap;
                }
            }
            else
            {
                overlaps.push_back(overlap);
            }
        }
    }

    if (uniqueExtensions)
    {
        for (auto &hashPair : bestOverlapsHash)
        {
            overlaps.push_back(hashPair.second);
        }
    }

    return overlaps;
}

std::vector<OverlapRange>
OverlapContainer::seqOverlaps(FastaRecord::Id seqId,
							  bool& outSuggestChimeric) const
{
	const FastaRecord& record = _queryContainer.getIndex().at(seqId);
	return _ovlpDetect.getSeqOverlaps(record, _onlyMax, outSuggestChimeric);
}


bool OverlapContainer::hasSelfOverlaps(FastaRecord::Id seqId)
{
	_indexMutex.lock();
	if (!_cached.count(seqId)) 
	{
		_indexMutex.unlock();
		this->lazySeqOverlaps(seqId);
		_indexMutex.lock();
	}
	bool selfOvlp = _suggestedChimeras.count(seqId);
	_indexMutex.unlock();
	return selfOvlp;
}

std::vector<OverlapRange>
	OverlapContainer::lazySeqOverlaps(FastaRecord::Id readId)
{
	_indexMutex.lock();
	if (!_cached.count(readId))
	{
		_indexMutex.unlock();
		bool suggestChimeric = false;
		auto overlaps = this->seqOverlaps(readId, suggestChimeric);
		_indexMutex.lock();
		this->storeOverlaps(overlaps, readId);
		if (suggestChimeric) 
		{
			_suggestedChimeras.insert(readId);
			_suggestedChimeras.insert(readId.rc());
		}
	}
	auto overlaps = _overlapIndex.at(readId);
	_indexMutex.unlock();
	return overlaps;
}

void OverlapContainer::storeOverlaps(const std::vector<OverlapRange>& overlaps,
									 FastaRecord::Id seqId)
{
	_cached.insert(seqId);
	_cached.insert(seqId.rc());

	auto& fwdOverlaps = _overlapIndex[seqId];
	auto& revOverlaps = _overlapIndex[seqId.rc()];

	std::unordered_set<FastaRecord::Id> extisting;
	if (_onlyMax)
	{
		for (auto& ovlp : fwdOverlaps) extisting.insert(ovlp.extId);
	}

	for (auto& ovlp : overlaps)
	{
		if (_onlyMax && extisting.count(ovlp.extId)) continue;

		auto revOvlp = ovlp.reverse();
		fwdOverlaps.push_back(ovlp);
		revOverlaps.push_back(ovlp.complement());
		_overlapIndex[revOvlp.curId].push_back(revOvlp);
		_overlapIndex[revOvlp.curId.rc()].push_back(revOvlp.complement());
	}
}

void OverlapContainer::findAllOverlaps()
{
	//Logger::get().info() << "Finding overlaps:";
	std::vector<FastaRecord::Id> allQueries;
	for (auto& hashPair : _queryContainer.getIndex())
	{
		allQueries.push_back(hashPair.first);
	}

	std::mutex indexMutex;
	std::function<void(const FastaRecord::Id&)> indexUpdate = 
	[this, &indexMutex] (const FastaRecord::Id& seqId)
	{
		auto& fastaRec = _queryContainer.getIndex().at(seqId);
		bool suggestChimeric = false;
		auto overlaps = _ovlpDetect.getSeqOverlaps(fastaRec, false, 
												   suggestChimeric);

		indexMutex.lock();
		this->storeOverlaps(overlaps, seqId);
		if (suggestChimeric) 
		{
			_suggestedChimeras.insert(seqId);
			_suggestedChimeras.insert(seqId.rc());
		}
		indexMutex.unlock();
	};

	processInParallel(allQueries, indexUpdate, 
					  Parameters::get().numThreads, true);

	int numOverlaps = 0;
	for (auto& seqOvlps : _overlapIndex) numOverlaps += seqOvlps.second.size();
	Logger::get().debug() << "Found " << numOverlaps << " overlaps";

	this->filterOverlaps();


	numOverlaps = 0;
	for (auto& seqOvlps : _overlapIndex) numOverlaps += seqOvlps.second.size();
	Logger::get().debug() << "Left " << numOverlaps 
		<< " overlaps after filtering";
}


void OverlapContainer::filterOverlaps()
{
	std::vector<FastaRecord::Id> seqIds;
	for (auto seqIndex : _queryContainer.getIndex())
	{
		seqIds.push_back(seqIndex.first);
	}

	std::function<void(const FastaRecord::Id& seqId)> filterParallel =
	[this] (const FastaRecord::Id& seqId)
	{
		auto& overlaps = _overlapIndex[seqId];
		
		std::vector<SetNode<OverlapRange*>*> overlapSets;
		for (auto& ovlp : overlaps) 
		{
			overlapSets.push_back(new SetNode<OverlapRange*>(&ovlp));
		}
		for (size_t i = 0; i < overlapSets.size(); ++i)
		{
			for (size_t j = 0; j < overlapSets.size(); ++j)
			{
				OverlapRange& ovlpOne = *overlapSets[i]->data;
				OverlapRange& ovlpTwo = *overlapSets[j]->data;

				if (ovlpOne.extId != ovlpTwo.extId) continue;
				int curDiff = ovlpOne.curRange() - ovlpOne.curIntersect(ovlpTwo);
				int extDiff = ovlpOne.extRange() - ovlpOne.extIntersect(ovlpTwo);

				if (curDiff < 100 && extDiff < 100) 
				{
					unionSet(overlapSets[i], overlapSets[j]);
				}
			}
		}
		std::unordered_map<SetNode<OverlapRange*>*, 
						   std::vector<OverlapRange>> clusters;
		for (auto& ovlp: overlapSets) 
		{
			clusters[findSet(ovlp)].push_back(*ovlp->data);
		}
		overlaps.clear();
		for (auto& cluster : clusters)
		{
			OverlapRange* maxOvlp = nullptr;
			for (auto& ovlp : cluster.second)
			{
				if (!maxOvlp || ovlp.score > maxOvlp->score)
				{
					maxOvlp = &ovlp;
				}
			}
			overlaps.push_back(*maxOvlp);
		}
		for (auto& ovlpNode : overlapSets) delete ovlpNode;

		std::sort(overlaps.begin(), overlaps.end(), 
				  [](const OverlapRange& o1, const OverlapRange& o2)
				  {return o1.curBegin < o2.curBegin;});

	};
	processInParallel(seqIds, filterParallel, 
					  Parameters::get().numThreads, false);
}


/*void OverlapContainer::saveOverlaps(const std::string& filename)
{
	std::ofstream fout(filename);
	if (!fout.is_open())
	{
		throw std::runtime_error("Can't open overlaps file");
	}

	for (auto& hashPair : _overlapIndex)
	{
		for (auto& ovlp : hashPair.second)
		{
			fout << ovlp.serialize() << std::endl;
		}
	}
}

void OverlapContainer::loadOverlaps(const std::string& filename)
{
	std::ifstream fin(filename);
	if (!fin.is_open())
	{
		throw std::runtime_error("Can't open overlaps file");
	}

	std::string buffer;
	while(!fin.eof())
	{
		std::getline(fin, buffer);
		if (buffer.empty()) break;
		OverlapRange ovlp;
		ovlp.unserialize(buffer);
		_overlapIndex[ovlp.curId].push_back(ovlp);
	}
}*/

void OverlapContainer::buildIntervalTree()
{
	Logger::get().debug() << "Building interval tree";
	for (auto& seqOvlps : _overlapIndex)
	{
		std::vector<Interval<OverlapRange*>> intervals;
		for (auto& ovlp : seqOvlps.second)
		{
			intervals.emplace_back(ovlp.curBegin, ovlp.curEnd, &ovlp);
		}
		_ovlpTree[seqOvlps.first] = IntervalTree<OverlapRange*>(intervals);
	}
}

std::vector<Interval<OverlapRange*>> 
	OverlapContainer::getOverlaps(FastaRecord::Id seqId, 
								  int32_t start, int32_t end) const
{
	return _ovlpTree.at(seqId).findOverlapping(start, end);
}
