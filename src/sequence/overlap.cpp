//(c) 2016 by Authors
//This file is a part of ABruijn program.
//Released under the BSD license (see LICENSE file)

#include <set>
#include <iostream>
#include <cassert>
#include <algorithm>
#include <thread>
#include <ctime>
#include <queue>

#include "overlap.h"
#include "../common/config.h"
#include "../common/utils.h"
#include "../common/parallel.h"
#include "../common/disjoint_set.h"


//reject overlaps early to speed everything up
/*bool OverlapDetector::goodStart(int32_t curPos, int32_t extPos, 
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
}*/


/*OverlapDetector::JumpRes 
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
}*/


//Check if it is a proper overlap
bool OverlapDetector::overlapTest(const OverlapRange& ovlp,
								  bool& outSuggestChimeric) const
{
	static const float OVLP_DIVERGENCE = Config::get("overlap_divergence_rate");
	if (ovlp.curRange() < _minOverlap || 
		ovlp.extRange() < _minOverlap) 
	{
		return false;
	}

	float lengthDiff = abs(ovlp.curRange() - ovlp.extRange());
	float meanLength = (ovlp.curRange() + ovlp.extRange()) / 2.0f;
	if (lengthDiff > meanLength * OVLP_DIVERGENCE)
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

/*namespace
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
}*/

namespace
{
	struct KmerMatch
	{
		KmerMatch(int32_t cur = 0, int32_t ext = 0,
				  FastaRecord::Id extId = FastaRecord::ID_NONE): 
			curPos(cur), extPos(ext), extId(extId) {}
		int32_t curPos;
		int32_t extPos;
		FastaRecord::Id extId;
	};
}


//This implementation was inspired by Hen Li's minimap2 paper
std::vector<OverlapRange> 
OverlapDetector::getSeqOverlaps(const FastaRecord& fastaRec, 
								bool uniqueExtensions,
								bool& outSuggestChimeric) const
{
	const int MAX_LOOK_BACK = 50;
	const int MIN_SURV_KMERS = _minOverlap / 100;
	const int MAX_EXT_SEQS = 1000;

	static float totalDpTime = 0;
	static float totalKmerTime = 0;
	static float totalHashTime = 0;
	clock_t begin = clock();

	outSuggestChimeric = false;
	int32_t curLen = fastaRec.sequence.length();

	std::vector<unsigned char> seqHitCount(_seqContainer.getMaxSeqId(), 0);

	std::vector<KmerMatch> vecMatches;
	vecMatches.reserve(10000000);

	for (auto curKmerPos : IterKmers(fastaRec.sequence))
	{
		if (!_vertexIndex.isSolid(curKmerPos.kmer)) continue;
		//bool isRepetitive = _vertexIndex.kmerFreq(curKmerPos.kmer) > 500;

		FastaRecord::Id prevSeqId = FastaRecord::ID_NONE;
		for (const auto& extReadPos : _vertexIndex.iterKmerPos(curKmerPos.kmer))
		//auto kmerRange = _vertexIndex.iterKmerPos(curKmerPos.kmer);
		//for (auto ptr = kmerRange.begin; ptr < kmerRange.end; ++ptr)
		{
			//no trivial matches
			if ((extReadPos.readId == fastaRec.id &&
				extReadPos.position == curKmerPos.position)) continue;

			/*auto extReadPos = *ptr;
			if (kmerRange.rc)
			{
				int32_t extLen = _seqContainer.seqLen(extReadPos.readId);
				extReadPos.readId = extReadPos.readId.rc();
				extReadPos.position = extLen - extReadPos.position - 
							   Parameters::get().kmerSize;
			}*/

			//count one seq match for one unique k-mer
			//since k-mers in vector are stored relative to fwd strand,
			//check both read orientations
			if (prevSeqId != extReadPos.readId &&
				prevSeqId != extReadPos.readId.rc())
			{
				if (seqHitCount[extReadPos.readId.rawId()] <
					std::numeric_limits<unsigned char>::max())
				{
					++seqHitCount[extReadPos.readId.rawId()];
				}
			}
			prevSeqId = extReadPos.readId;

			vecMatches.emplace_back(curKmerPos.position, 
									extReadPos.position,
									extReadPos.readId);
		}
	}
	auto hashTime = clock();
	totalHashTime += double(hashTime - begin) / CLOCKS_PER_SEC;

	//std::vector<std::vector<KmerMatch>*> 
	//	seqMatches(_seqContainer.getMaxSeqId(), 0);
	cuckoohash_map<FastaRecord::Id, std::vector<KmerMatch>*> seqMatches;
	seqMatches.reserve(500);
	for (auto& match : vecMatches)
	{
		if (seqHitCount[match.extId.rawId()] < MIN_SURV_KMERS) continue;

		seqMatches.upsert(match.extId, 
						  [&seqHitCount, &match](std::vector<KmerMatch>*& vec)
						  {
							  if (!vec)
							  {
								  vec = new std::vector<KmerMatch>;
								  vec->reserve(seqHitCount[match.extId.rawId()]);
							  }
							  vec->push_back(match);
						  }, nullptr);

		/*if (!seqMatches[match.extId.rawId()])
		{
			seqMatches[match.extId.rawId()] = new std::vector<KmerMatch>();
			seqMatches[match.extId.rawId()]
				->reserve(seqHitCount[match.extId.rawId()]);
		}
		seqMatches[match.extId.rawId()]->push_back(match);*/
	}
	
  	clock_t end = clock();
  	double elapsed_secs = double(end - hashTime) / CLOCKS_PER_SEC;
	totalKmerTime += elapsed_secs;

	std::vector<OverlapRange> detectedOverlaps;
	int uniqueCandidates = 0;
	double sumMeanRepeatRatio = 0;
	int64_t numMeanRepeatRatio = 0;
	//for (auto& sm : seqMatches)
	for (auto& seqVec : seqMatches.lock_table())
	{
		//if (!sm) continue;
		const std::vector<KmerMatch>& matchesList = *seqVec.second;
		int32_t extLen = _seqContainer.seqLen(seqVec.first);

		//pre-filtering
		int32_t minCur = matchesList.front().curPos;
		int32_t maxCur = matchesList.back().curPos;
		int32_t minExt = std::numeric_limits<int32_t>::max();
		int32_t maxExt = std::numeric_limits<int32_t>::min();
		int32_t uniquePos = 0;
		int32_t prevPos = -1;
		for (auto& match : matchesList)
		{
			minExt = std::min(minExt, match.extPos);
			maxExt = std::max(maxExt, match.extPos);
			if (match.curPos != prevPos)
			{
				prevPos = match.curPos;
				++uniquePos;
			}
		}
		float repeatRate = (float)matchesList.size() / uniquePos;
		/*if (repeatRate > 1.1)
		{
			Logger::get().debug() << "\t" << repeatRate;
		}*/
		sumMeanRepeatRatio += repeatRate;
		++numMeanRepeatRatio;
		if (maxCur - minCur < _minOverlap || 
			maxExt - minExt < _minOverlap) continue;
		if (_checkOverhang)
		{
			if (std::min(minCur, minExt) > _maxOverhang) continue;
			if (std::min(curLen - maxCur, 
						 extLen - maxExt) > _maxOverhang) continue;
		}
		if (++uniqueCandidates > MAX_EXT_SEQS) break;
		
		//chain matiching positions with DP
		std::vector<int32_t> scoreTable(matchesList.size(), 0);
		std::vector<int32_t> backtrackTable(matchesList.size(), -1);
		int32_t skipCurPos = 0;
		int32_t skipCurId = 0;
		for (int32_t i = 1; i < (int32_t)scoreTable.size(); ++i)
		{
			int32_t maxScore = 0;
			int32_t maxId = 0;
			int32_t curNext = matchesList[i].curPos;
			int32_t extNext = matchesList[i].extPos;
			int32_t noImprovement = 0;

			if (curNext != skipCurPos)
			{
				skipCurPos = curNext;
				skipCurId = i - 1;
			}

			for (int32_t j = skipCurId; j >= 0; --j)
			{
				int32_t nextScore = 0;
				int32_t curPrev = matchesList[j].curPos;
				int32_t extPrev = matchesList[j].extPos;
				if (0 < curNext - curPrev && curNext - curPrev < _maxJump &&
					0 < extNext - extPrev && extNext - extPrev < _maxJump)
				{
					int32_t matchScore = std::min(curNext - curPrev, 
												  extNext - extPrev);
					int32_t jumpDiv = abs((curNext - curPrev) - 
										  (extNext - extPrev));
					int32_t gapCost = jumpDiv ? 0.1f * jumpDiv + log2(jumpDiv) : 0;

					nextScore = scoreTable[j] + matchScore - gapCost;

					if (nextScore > maxScore)
					{
						maxScore = nextScore;
						maxId = j;
						noImprovement = 0;
					}
					else
					{
						if (++noImprovement > MAX_LOOK_BACK) break;
					}
				}
				if (curNext - curPrev > _maxJump) break;
			}

			if (maxScore > 0)
			{
				scoreTable[i] = maxScore;
				backtrackTable[i] = maxId;
			}
		}

		//backtracking
		std::vector<OverlapRange> extOverlaps;
		for (int32_t chainStart = backtrackTable.size() - 1; 
			 chainStart > 0; --chainStart)
		{
			int32_t pos = chainStart;
			std::vector<int32_t> chainMatches;
			chainMatches.reserve(matchesList.size());
			std::vector<int32_t> shifts;
			shifts.reserve(matchesList.size());
			while (pos != -1)
			{
				chainMatches.push_back(pos);
				shifts.push_back(matchesList[pos].curPos - 
								 matchesList[pos].extPos);
				int32_t newPos = backtrackTable[pos];
				backtrackTable[pos] = -1;
				pos = newPos;
			}

			OverlapRange ovlp(fastaRec.id, matchesList.front().extId,
							  matchesList[chainMatches.back()].curPos, 
							  matchesList[chainMatches.back()].extPos,
							  curLen, extLen);
			ovlp.curEnd = matchesList[chainMatches.front()].curPos;
			ovlp.extEnd = matchesList[chainMatches.front()].extPos;
			ovlp.leftShift = median(shifts);
			ovlp.rightShift = extLen - curLen + ovlp.leftShift;
			ovlp.score = scoreTable[chainStart];
			if (this->overlapTest(ovlp, outSuggestChimeric))
			{
				extOverlaps.push_back(ovlp);
			}
		}
		
		OverlapRange* maxOvlp = nullptr;
		bool passedTest = false;
		for (auto& ovlp : extOverlaps)
		{
			if (!uniqueExtensions)
			{
				detectedOverlaps.push_back(ovlp);
			}
			else
			{
				passedTest = true;
				if (!maxOvlp || ovlp.score > maxOvlp->score)
				{
					maxOvlp = &ovlp;
				}
			}
		}
		if (uniqueExtensions && passedTest)
		{
			detectedOverlaps.push_back(*maxOvlp);
		}
	}

	clock_t ff = clock();
	double es = double(ff - end) / CLOCKS_PER_SEC;
	totalDpTime += es;
	/*Logger::get().debug() << "---------";
	Logger::get().debug() << " " << vecMatches.size() << " " 
		<< uniqueCandidates << " " << detectedOverlaps.size();
	Logger::get().debug() << "Mean repeat ratio: " 
		<< sumMeanRepeatRatio / (numMeanRepeatRatio + 1);
	Logger::get().debug() << "Matches per sequence: " 
		<< vecMatches.size() / (seqMatches.size() + 1);
	Logger::get().debug() << "hash: " << totalHashTime << " k-mer: " 
		<< totalKmerTime << " dp: " << totalDpTime;*/

	//for (auto& kmerCounts : seqMatches) delete kmerCounts;
	for (auto& kmerCounts : seqMatches.lock_table()) delete kmerCounts.second;
	return detectedOverlaps;
}

/*std::vector<OverlapRange> 
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
		//if (_vertexIndex.isRepetitive(curKmerPos.kmer)) continue;

		int32_t curPos = curKmerPos.position;

		//for all other occurences of this kmer (extension candidates)
		for (const auto& extReadPos : _vertexIndex.iterKmerPos(curKmerPos.kmer))
		{
			//no trivial matches
			if (extReadPos.readId == fastaRec.id &&
				extReadPos.position == curPos) continue;

			//don't process highly repetitive k-mers unless
			//we are continuing paths
			//if (_vertexIndex.isRepetitive(curKmerPos.kmer) &&
			//	!activePaths.count(extReadPos.readId)) continue;

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
				!_vertexIndex.isRepetitive(curKmerPos.kmer) &&
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
}*/

std::vector<OverlapRange>
OverlapContainer::seqOverlaps(FastaRecord::Id seqId,
							  bool& outSuggestChimeric) const
{
	const FastaRecord& record = _queryContainer.getRecord(seqId);
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
	for (auto& seq : _queryContainer.iterSeqs())
	{
		allQueries.push_back(seq.id);
	}

	std::mutex indexMutex;
	std::function<void(const FastaRecord::Id&)> indexUpdate = 
	[this, &indexMutex] (const FastaRecord::Id& seqId)
	{
		auto& fastaRec = _queryContainer.getRecord(seqId);
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
	const int MAX_ENDS_DIFF = 100;

	std::vector<FastaRecord::Id> seqIds;
	for (auto& seq : _queryContainer.iterSeqs())
	{
		seqIds.push_back(seq.id);
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

				if (curDiff < MAX_ENDS_DIFF && extDiff < MAX_ENDS_DIFF) 
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
