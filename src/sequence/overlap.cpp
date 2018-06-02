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
	outSuggestChimeric = false;

	std::unordered_map<FastaRecord::Id, 
					   std::vector<DPRecord>> activePaths;
	std::unordered_map<FastaRecord::Id, 
					   std::vector<DPRecord>> completedPaths;
	std::set<size_t> eraseMarks;
	int32_t curLen = fastaRec.sequence.length();

	//for all kmers in this read
	for (auto curKmerPos : IterKmers(fastaRec.sequence))
	{
		if (!_vertexIndex.isSolid(curKmerPos.kmer)) continue;

		int32_t curPos = curKmerPos.position;

		//for all other occurences of this kmer (extension candidates)
		for (const auto& extReadPos : _vertexIndex.iterKmerPos(curKmerPos.kmer))
		{
			//no trivial matches
			if (extReadPos.readId == fastaRec.id &&
				extReadPos.position == curPos) continue;

			int32_t extLen = _seqContainer.seqLen(extReadPos.readId);
			if (extLen < (int32_t)_minOverlap) continue;

			int32_t extPos = extReadPos.position;
			auto& extPaths = activePaths[extReadPos.readId];

			size_t maxCloseId = 0;
			size_t maxFarId = 0;
			int32_t maxCloseScore = 0;
			int32_t maxFarScore = 0;
			bool extendsClose = false;
			bool extendsFar = false;
			eraseMarks.clear();

			//searching for longest possible extension
			for (size_t pathId = 0; pathId < extPaths.size(); ++pathId)
			{
				JumpRes jumpResult = 
					this->jumpTest(extPaths[pathId].ovlp.curEnd, curPos,
								   extPaths[pathId].ovlp.extEnd, extPos);

				static const int PENALTY_WND = Config::get("penalty_window");
				int32_t jumpLength = curPos - extPaths[pathId].ovlp.curEnd;
				int32_t gapScore = -(jumpLength - _gapSize) / PENALTY_WND;
				if (jumpLength < _gapSize) gapScore = 1;

				int32_t jumpScore = extPaths[pathId].ovlp.score + gapScore;

				switch (jumpResult)
				{
					case J_END:
						eraseMarks.insert(pathId);
						//if (this->overlapTest(extPaths[pathId].ovlp, 
						//					  outSuggestChimeric) &&
						//	!extPaths[pathId].wasContinued)
						if (this->overlapTest(extPaths[pathId].ovlp, 
											  outSuggestChimeric))
						{
							completedPaths[extReadPos.readId]
									.push_back(extPaths[pathId]);
						}
						break;
					case J_INCONS:
						break;
					case J_CLOSE:
						eraseMarks.insert(pathId);
						if (jumpScore > maxCloseScore)
						{
							extendsClose = true;
							maxCloseId = pathId;	
							maxCloseScore = jumpScore;
						}
						break;
					case J_FAR:
						if (jumpScore > maxFarScore)
						{
							extendsFar = true;
							maxFarId = pathId;
							maxFarScore = jumpScore;
						}
						break;
				}
			}
			//update the best close extension
			if (extendsClose)
			{
				eraseMarks.erase(maxCloseId);
				//extPaths[maxCloseId].wasContinued = false;
				extPaths[maxCloseId].ovlp.curEnd = curPos;
				extPaths[maxCloseId].ovlp.extEnd = extPos;
				extPaths[maxCloseId].ovlp.score = maxCloseScore;
				extPaths[maxCloseId].shifts.push_back(curPos - extPos);

				if (_keepAlignment)
				{
					auto& kmerMatches = extPaths[maxCloseId].ovlp.kmerMatches;
					if (kmerMatches.empty() || curPos - kmerMatches.back().first > 
											   (int32_t)Parameters::get().kmerSize)
					{
						kmerMatches.emplace_back(curPos, extPos);
					}
				}
			}
			//update the best far extension, keep the old path as a copy
			if (extendsFar)
			{
				//extPaths[maxFarId].wasContinued = true;
				extPaths.push_back(extPaths[maxFarId]);
				//extPaths.back().wasContinued = false;
				extPaths.back().ovlp.curEnd = curPos;
				extPaths.back().ovlp.extEnd = extPos;
				extPaths.back().ovlp.score = maxFarScore;
				extPaths.back().shifts.push_back(curPos - extPos);

				if (_keepAlignment)
				{
					auto& kmerMatches = extPaths.back().ovlp.kmerMatches;
					if (kmerMatches.empty() || curPos - kmerMatches.back().first > 
											   (int32_t)Parameters::get().kmerSize)
					{
						kmerMatches.emplace_back(curPos, extPos);
					}
				}
			}
			//if no extensions possible (or there are no active paths), start a new path
			if (!extendsClose && !extendsFar &&
				this->goodStart(curPos, extPos, curLen, extLen, 
								fastaRec.id, extReadPos.readId))
			{
				OverlapRange ovlp(fastaRec.id, extReadPos.readId,
								  curPos, extPos, curLen, extLen);
				extPaths.emplace_back(ovlp);
			}
			//cleaning up
			for (auto itEraseId = eraseMarks.rbegin(); 
				 itEraseId != eraseMarks.rend(); ++itEraseId)
			{
				extPaths[*itEraseId] = extPaths.back();
				extPaths.pop_back();
			}
		} //end loop over kmer occurences in other reads
	} //end loop over kmers in the current read

	//copy to coplete paths
	for (auto& ap : activePaths)
	{
		for (auto& dpRec : ap.second)
		{
			if (this->overlapTest(dpRec.ovlp, outSuggestChimeric))
			{
				completedPaths[ap.first].push_back(dpRec);
			}
		}
	}

	//leave only one overlap for each starting position
	for (auto& ap : completedPaths)
	{
		std::unordered_map<std::pair<int32_t, int32_t>, 
						   DPRecord, pairhash> maxByStart;
		for (auto& dpRec : ap.second)
		{
			auto& curMax = maxByStart[std::make_pair(dpRec.ovlp.curBegin, 
													 dpRec.ovlp.extBegin)];
			if (dpRec.ovlp.score > curMax.ovlp.score)
			{
				curMax = dpRec;
			}
		}
		ap.second.clear();
		for (auto& maxDp : maxByStart) ap.second.push_back(maxDp.second);
	}

	std::vector<OverlapRange> detectedOverlaps;
	for (auto& ap : completedPaths)
	{
		int32_t extLen = _seqContainer.seqLen(ap.first);
		DPRecord* maxRecord = nullptr;
		bool passedTest = false;
		for (auto& dpRec : ap.second)
		{
			if (!uniqueExtensions)
			{
				dpRec.ovlp.leftShift = median(dpRec.shifts);
				dpRec.ovlp.rightShift = extLen - curLen + 
										dpRec.ovlp.leftShift;
				detectedOverlaps.push_back(dpRec.ovlp);
			}
			else
			{
				passedTest = true;
				if (!maxRecord || dpRec.ovlp.score > maxRecord->ovlp.score)
				{
					maxRecord = &dpRec;
				}
			}
		}
		if (uniqueExtensions && passedTest)
		{
			maxRecord->ovlp.leftShift = median(maxRecord->shifts);
			maxRecord->ovlp.rightShift = extLen - curLen + 
									maxRecord->ovlp.leftShift;
			detectedOverlaps.push_back(maxRecord->ovlp);
		}
	}

	return detectedOverlaps;
}

bool OverlapContainer::hasSelfOverlaps(FastaRecord::Id readId)
{
	this->lazySeqOverlaps(readId);
	if (!readId.strand()) readId = readId.rc();
	return ((IndexVecWrapper)_overlapIndex[readId]).suggestChimeric;
}

const std::vector<OverlapRange>&
	OverlapContainer::lazySeqOverlaps(FastaRecord::Id readId)
{
	bool flipped = !readId.strand();
	if (flipped) readId = readId.rc();
	IndexVecWrapper wrapper;

	//upsert creates default value if it does not exist
	_overlapIndex.upsert(readId, 	
		[&wrapper](IndexVecWrapper& val)
			{wrapper = val;});
	if (wrapper.cached)
	{
		return !flipped ? *wrapper.fwdOverlaps : *wrapper.revOverlaps;
	}

	//otherwise, need to compute overlaps.
	//do it for forward strand to be distinct
	bool suggestChimeric;
	const FastaRecord& record = _queryContainer.getIndex().at(readId);
	auto overlaps = _ovlpDetect.getSeqOverlaps(record, _onlyMax, suggestChimeric);

	std::vector<OverlapRange> revOverlaps;
	revOverlaps.reserve(overlaps.size());
	for (auto& ovlp : overlaps) revOverlaps.push_back(ovlp.complement());

	_overlapIndex.update_fn(readId,
		[&wrapper, &overlaps, &revOverlaps, &suggestChimeric, &flipped]
		(IndexVecWrapper& val)
		{
			if (!val.cached)
			{
				*val.fwdOverlaps = std::move(overlaps);
				*val.revOverlaps = std::move(revOverlaps);
				val.suggestChimeric = suggestChimeric;
				val.cached = true;
			}
			wrapper = val;
		});

	return !flipped ? *wrapper.fwdOverlaps : *wrapper.revOverlaps;
}

void OverlapContainer::ensureTransitivity()
{
	Logger::get().debug() << "Computing transitive closure for overlaps";
	
	std::vector<FastaRecord::Id> allSeqs;
	for (auto& seqIt : _overlapIndex.lock_table()) 
	{
		allSeqs.push_back(seqIt.first);
		allSeqs.push_back(seqIt.first.rc());
	}

	int totalOverlaps = 0;
	for (auto& seq : allSeqs)
	{
		auto& curOvlps = this->unsafeSeqOverlaps(seq);
		totalOverlaps += curOvlps.size();
		for (auto& curOvlp : curOvlps)
		{
			auto& extOvlps = this->unsafeSeqOverlaps(curOvlp.extId);

			if (_onlyMax)
			{
				bool found = false;
				for (auto& extOvlp : extOvlps)
				{
					if (extOvlp.extId == curOvlp.curId)
					{
						if (curOvlp.score > extOvlp.score)
						{
							extOvlp = curOvlp.reverse();
						}
						found = true;
						break;
					}
				}
				if (!found)
				{
					extOvlps.push_back(curOvlp.reverse());
				}
			}
			else
			{
				extOvlps.push_back(curOvlp.reverse());
			}
		}
	}
	//Logger::get().debug() << "Total overlaps: " << totalOverlaps;
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
		this->lazySeqOverlaps(seqId);
	};
	processInParallel(allQueries, indexUpdate, 
					  Parameters::get().numThreads, true);
	this->ensureTransitivity();

	int numOverlaps = 0;
	for (auto& seqOvlps : _overlapIndex.lock_table()) 
	{
		numOverlaps += seqOvlps.second.fwdOverlaps->size() * 2;
	}
	Logger::get().debug() << "Found " << numOverlaps << " overlaps";

	this->filterOverlaps();

	numOverlaps = 0;
	for (auto& seqOvlps : _overlapIndex.lock_table()) 
	{
		numOverlaps += seqOvlps.second.fwdOverlaps->size() * 2;
	}
	Logger::get().debug() << "Left " << numOverlaps 
		<< " overlaps after filtering";
}

std::vector<OverlapRange>&
	OverlapContainer::unsafeSeqOverlaps(FastaRecord::Id seqId)
{
		FastaRecord::Id normId = seqId.strand() ? seqId : seqId.rc();
		_overlapIndex.insert(normId);	//ensure it's in the table
		IndexVecWrapper wrapper = _overlapIndex[normId];
		return seqId.strand() ? *wrapper.fwdOverlaps : 
								*wrapper.revOverlaps;
}

//TODO: potentially might become non-symmetric after filtering
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
		auto& overlaps = this->unsafeSeqOverlaps(seqId);
		
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


void OverlapContainer::buildIntervalTree()
{
	//Logger::get().debug() << "Building interval tree";
	std::vector<FastaRecord::Id> allSeqs;
	for (auto& seqIt : _overlapIndex.lock_table()) 
	{
		allSeqs.push_back(seqIt.first);
		allSeqs.push_back(seqIt.first.rc());
	}

	for (auto& seq : allSeqs)
	{
		std::vector<Interval<OverlapRange*>> intervals;
		auto& overlaps = this->unsafeSeqOverlaps(seq);
		for (auto& ovlp : overlaps)
		{
			intervals.emplace_back(ovlp.curBegin, ovlp.curEnd, &ovlp);
		}
		_ovlpTree[seq] = IntervalTree<OverlapRange*>(intervals);
	}
}

std::vector<Interval<OverlapRange*>> 
	OverlapContainer::getCoveringOverlaps(FastaRecord::Id seqId, 
								  int32_t start, int32_t end) const
{
	return _ovlpTree.at(seqId).findOverlapping(start, end);
}
