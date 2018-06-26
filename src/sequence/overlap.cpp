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
#include <cmath>

#include "overlap.h"
#include "../common/config.h"
#include "../common/utils.h"
#include "../common/parallel.h"
#include "../common/disjoint_set.h"


//Check if it is a proper overlap
bool OverlapDetector::overlapTest(const OverlapRange& ovlp,
								  bool& outSuggestChimeric) const
{
	if (ovlp.curRange() < _minOverlap || 
		ovlp.extRange() < _minOverlap) 
	{
		return false;
	}

	static const float OVLP_DIV = 0.5;
	float lengthDiff = abs(ovlp.curRange() - ovlp.extRange());
	//float meanLength = (ovlp.curRange() + ovlp.extRange()) / 2.0f;
	if (lengthDiff > OVLP_DIV * std::min(ovlp.curRange(), ovlp.extRange()))
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
	struct KmerMatch
	{
		KmerMatch(int32_t cur = 0, int32_t ext = 0,
				  FastaRecord::Id extId = FastaRecord::ID_NONE): 
			curPos(cur), extPos(ext), extId(extId) {}
		int32_t curPos;
		int32_t extPos;
		FastaRecord::Id extId;
	};

	struct MatchVecWrapper
	{
		MatchVecWrapper(){}
		MatchVecWrapper(const KmerMatch& match, size_t capacity):
			v(new std::vector<KmerMatch>)
		{
			v->reserve(capacity);
			v->push_back(match);
		}
		std::shared_ptr<std::vector<KmerMatch>> v;
		std::vector<KmerMatch>* operator->() {return v.get();}
		std::vector<KmerMatch>& operator*() {return *v;}
	};
}


//This implementation was inspired by Hen Li's minimap2 paper
//might be used in parallel
std::vector<OverlapRange> 
OverlapDetector::getSeqOverlaps(const FastaRecord& fastaRec, 
								bool& outSuggestChimeric) const
{
	const float MIN_KMER_SURV_RATE = 0.01;	//TODO: put into config
	const int MAX_LOOK_BACK = 50;
	const int kmerSize = Parameters::get().kmerSize;

	//static std::ofstream fout("../kmers.txt");

	//static float totalDpTime = 0;
	//static float totalKmerTime = 0;
	//static float totalHashTime = 0;
	//static float totalDpLoop = 0;
	//static float totalBackLoop = 0;
	//clock_t begin = clock();

	outSuggestChimeric = false;
	int32_t curLen = fastaRec.sequence.length();

	static thread_local cuckoohash_map<FastaRecord::Id, uint16_t> seqHitCount;
	auto lockedHitCount = seqHitCount.lock_table();
	thread_local static std::vector<KmerMatch> vecMatches;
	thread_local static std::vector<Kmer> curKmers;

	static thread_local std::vector<Kmer> curOvlpKmers;
	static thread_local std::vector<Kmer> extOvlpKmers;
	static thread_local std::vector<Kmer> intersect;

	static thread_local 
		std::vector<std::pair<FastaRecord::Id, uint16_t>> topSeqs;

	static thread_local std::vector<KmerMatch> matchesList;
	static thread_local std::vector<int32_t> scoreTable;
	static thread_local std::vector<int32_t> backtrackTable;

	static clock_t prevTime = clock();	//intentionally shared
	clock_t timeNow = clock();
	if (double(timeNow - prevTime) / CLOCKS_PER_SEC > 100)
	{
		prevTime = timeNow;
		Logger::get().debug() << "----------";
		Logger::get().debug() << "seqHitCount: " << lockedHitCount.size() << " " << lockedHitCount.capacity();
		Logger::get().debug() << "vecMatches: " << vecMatches.size() << " " << vecMatches.capacity();
		Logger::get().debug() << "curKmers: " << curKmers.size() << " " << curKmers.capacity();
		Logger::get().debug() << "curOvlpKmers: " << curOvlpKmers.size() << " " << curOvlpKmers.capacity();
		Logger::get().debug() << "extOvlpKmers: " << extOvlpKmers.size() << " " << extOvlpKmers.capacity();
		Logger::get().debug() << "intersect: " << intersect.size() << " " << intersect.capacity();
		Logger::get().debug() << "topSeqs: " << topSeqs.size() << " " << topSeqs.capacity();
		Logger::get().debug() << "matchesList: " << matchesList.size() << " " << matchesList.capacity();
		Logger::get().debug() << "scoreTable: " << scoreTable.size() << " " << scoreTable.capacity();
		Logger::get().debug() << "backtrackTable: " << backtrackTable.size() << " " << backtrackTable.capacity();
		Logger::get().debug() << "";
	}

	lockedHitCount.clear();
	vecMatches.clear();
	curKmers.clear();

	//count kmer hits
	for (auto curKmerPos : IterKmers(fastaRec.sequence))
	{
		curKmers.push_back(curKmerPos.kmer);
		if (!_vertexIndex.isSolid(curKmerPos.kmer)) continue;

		FastaRecord::Id prevSeqId = FastaRecord::ID_NONE;
		for (const auto& extReadPos : _vertexIndex.iterKmerPos(curKmerPos.kmer))
		{
			//no trivial matches
			if ((extReadPos.readId == fastaRec.id &&
				extReadPos.position == curKmerPos.position)) continue;

			//count one seq match for one unique k-mer
			//since k-mers in vector are stored relative to fwd strand,
			//check both read orientations
			if (prevSeqId != extReadPos.readId &&
				prevSeqId != extReadPos.readId.rc())
			{
				if (lockedHitCount[extReadPos.readId] <
					std::numeric_limits<unsigned char>::max())
				{
					++lockedHitCount[extReadPos.readId];
				}
			}
			prevSeqId = extReadPos.readId;
		}
	}
	//auto hashTime = clock();
	//totalHashTime += double(hashTime - begin) / CLOCKS_PER_SEC;

	//if there is a limit on the number of sequences to consider,
	//sort by the decreasing number of k-mer hits and filter
	//some out if needed
	if (_maxCurOverlaps > 0)
	{
		topSeqs.clear();
		for (auto seqCount : lockedHitCount)
		{
			if (seqCount.second >= MIN_KMER_SURV_RATE * _minOverlap)
			{
				topSeqs.emplace_back(seqCount.first, seqCount.second);
			}
		}
		std::sort(topSeqs.begin(), topSeqs.end(),
				  [](const std::pair<FastaRecord::Id, uint16_t> p1, 
					 const std::pair<FastaRecord::Id, uint16_t>& p2) 
					 {return p1.second > p2.second;});
		for (size_t i = (size_t)_maxCurOverlaps; i < topSeqs.size(); ++i)
		{
			lockedHitCount[topSeqs[i].first] = 0;
		}
	}

	//finally full the vector matches
	for (auto curKmerPos : IterKmers(fastaRec.sequence))
	{
		if (!_vertexIndex.isSolid(curKmerPos.kmer)) continue;
		for (const auto& extReadPos : _vertexIndex.iterKmerPos(curKmerPos.kmer))
		{
			//no trivial matches
			if ((extReadPos.readId == fastaRec.id &&
				extReadPos.position == curKmerPos.position)) continue;

			if (lockedHitCount[extReadPos.readId] >= 
				MIN_KMER_SURV_RATE * _minOverlap)
			{
				vecMatches.emplace_back(curKmerPos.position, 
										extReadPos.position,
										extReadPos.readId);
			}
		}
	}
	//group by extId
	std::stable_sort(vecMatches.begin(), vecMatches.end(),
			  		 [](const KmerMatch& k1, const KmerMatch& k2)
			  		 {return k1.extId < k2.extId;});
	
  	//clock_t end = clock();
  	//double elapsed_secs = double(end - hashTime) / CLOCKS_PER_SEC;
	//totalKmerTime += elapsed_secs;

	std::vector<OverlapRange> detectedOverlaps;
	//int uniqueCandidates = 0;
	size_t extRangeBegin = 0;
	size_t extRangeEnd = 0;
	while(true)
	{
		extRangeBegin = extRangeEnd;
		while (extRangeEnd < vecMatches.size() &&
			   vecMatches[extRangeBegin].extId == 
			   vecMatches[extRangeEnd].extId) ++extRangeEnd;
		if (extRangeEnd == vecMatches.size()) break;

		matchesList.assign(vecMatches.begin() + extRangeBegin,
						   vecMatches.begin() + extRangeEnd);
		FastaRecord::Id extId = vecMatches[extRangeBegin].extId;
		int32_t extLen = _seqContainer.seqLen(extId);

		//pre-filtering
		int32_t minCur = matchesList.front().curPos;
		int32_t maxCur = matchesList.back().curPos;
		int32_t minExt = std::numeric_limits<int32_t>::max();
		int32_t maxExt = std::numeric_limits<int32_t>::min();
		//int32_t uniquePos = 0;
		//int32_t prevPos = -1;
		for (auto& match : matchesList)
		{
			minExt = std::min(minExt, match.extPos);
			maxExt = std::max(maxExt, match.extPos);
			/*if (match.curPos != prevPos)
			{
				prevPos = match.curPos;
				++uniquePos;
			}*/
		}
		if (maxCur - minCur < _minOverlap || 
			maxExt - minExt < _minOverlap) continue;
		if (_checkOverhang)
		{
			if (std::min(minCur, minExt) > _maxOverhang) continue;
			if (std::min(curLen - maxCur, 
						 extLen - maxExt) > _maxOverhang) continue;
		}
		//++uniqueCandidates;

		//chain matiching positions with DP
		scoreTable.assign(matchesList.size(), 0);
		backtrackTable.assign(matchesList.size(), -1);

		bool extSorted = extLen > curLen;
		if (extSorted)
		{
			std::sort(matchesList.begin(), matchesList.end(),
					  [](const KmerMatch& k1, const KmerMatch& k2)
					  {return k1.extPos < k2.extPos;});
		}

  		//clock_t dpBegin = clock();
		for (int32_t i = 1; i < (int32_t)scoreTable.size(); ++i)
		{
			int32_t maxScore = 0;
			int32_t maxId = 0;
			int32_t curNext = matchesList[i].curPos;
			int32_t extNext = matchesList[i].extPos;
			int32_t noImprovement = 0;

			for (int32_t j = i - 1; j >= 0; --j)
			{
				int32_t curPrev = matchesList[j].curPos;
				int32_t extPrev = matchesList[j].extPos;
				if (0 < curNext - curPrev && curNext - curPrev < _maxJump &&
					0 < extNext - extPrev && extNext - extPrev < _maxJump)
				{
					int32_t matchScore = 
						std::min(std::min(curNext - curPrev, extNext - extPrev),
										  kmerSize);
					int32_t jumpDiv = abs((curNext - curPrev) - 
										  (extNext - extPrev));
					int32_t gapCost = jumpDiv ? 
							0.01f * kmerSize * jumpDiv + std::log2(jumpDiv) : 0;
					int32_t nextScore = scoreTable[j] + matchScore - gapCost;
					if (nextScore > maxScore)
					{
						maxScore = nextScore;
						maxId = j;
						noImprovement = 0;

						if (jumpDiv == 0 && curNext - curPrev < kmerSize) break;
					}
					else
					{
						if (++noImprovement > MAX_LOOK_BACK) break;
					}
				}
				if (extSorted && extNext - extPrev > _maxJump) break;
				if (!extSorted && curNext - curPrev > _maxJump) break;
			}

			scoreTable[i] = std::max(maxScore, kmerSize);
			if (maxScore > kmerSize)
			{
				backtrackTable[i] = maxId;
			}
		}
		//clock_t dpEnd = clock();
  		//totalDpLoop += double(dpEnd - dpBegin) / CLOCKS_PER_SEC;

		//backtracking
		std::vector<OverlapRange> extOverlaps;
		std::vector<int32_t> shifts;
		std::vector<std::pair<int32_t, int32_t>> kmerMatches;

		for (int32_t chainStart = backtrackTable.size() - 1; 
			 chainStart > 0; --chainStart)
		{
			if (backtrackTable[chainStart] == -1) continue;

			int32_t pos = chainStart;
			KmerMatch lastMatch = matchesList[pos];
			KmerMatch firstMatch = lastMatch;

			//int chainLength = 0;
			shifts.clear();
			kmerMatches.clear();
			//int totalMatch = kmerSize;
			while (pos != -1)
			{
				firstMatch = matchesList[pos];
				shifts.push_back(matchesList[pos].curPos - 
								 matchesList[pos].extPos);
				//++chainLength;

				/*int32_t prevPos = backtrackTable[pos];
				if (prevPos != -1)
				{
					int32_t curNext = matchesList[pos].curPos;
					int32_t extNext = matchesList[pos].extPos;
					int32_t curPrev = matchesList[prevPos].curPos;
					int32_t extPrev = matchesList[prevPos].extPos;
					int32_t matchScore = 
							std::min(std::min(curNext - curPrev, extNext - extPrev),
											  kmerSize);
					totalMatch += matchScore;
				}*/
				if (_keepAlignment)
				{
					if (kmerMatches.empty() || 
						kmerMatches.back().first - matchesList[pos].curPos >
						kmerSize)
					{
						kmerMatches.emplace_back(matchesList[pos].curPos,
								 				 matchesList[pos].extPos);
					}
				}

				int32_t newPos = backtrackTable[pos];
				backtrackTable[pos] = -1;
				pos = newPos;
			}

			OverlapRange ovlp(fastaRec.id, matchesList.front().extId,
							  firstMatch.curPos, firstMatch.extPos,
							  curLen, extLen);
			ovlp.curEnd = lastMatch.curPos + kmerSize - 1;
			ovlp.extEnd = lastMatch.extPos + kmerSize - 1;
			ovlp.score = scoreTable[chainStart];

			if (this->overlapTest(ovlp, outSuggestChimeric))
			{
				std::reverse(kmerMatches.begin(), kmerMatches.end());
				ovlp.kmerMatches = kmerMatches;
				ovlp.leftShift = median(shifts);
				ovlp.rightShift = extLen - curLen + ovlp.leftShift;
				extOverlaps.push_back(ovlp);
			}
		}
		
		//selecting the best
		std::vector<OverlapRange> ovlpCandidates;
		if (_onlyMaxExt)
		{
			OverlapRange* maxOvlp = nullptr;
			for (auto& ovlp : extOverlaps)
			{
				if (!maxOvlp || ovlp.score > maxOvlp->score)
				{
					maxOvlp = &ovlp;
				}
			}
			if (maxOvlp) ovlpCandidates.push_back(*maxOvlp);
		}
		else
		{
			//sort by decreasing score
			std::sort(extOverlaps.begin(), extOverlaps.end(),
					  [](const OverlapRange& r1, const OverlapRange& r2)
					  {return r1.score > r2.score;});
			
			for (auto& ovlp : extOverlaps)
			{
				bool isContained = false;
				for (auto& prim : ovlpCandidates)
				{
					if (ovlp.containedBy(prim))
					{
						isContained = true;
						break;
					}
				}
				if (!isContained)
				{
					ovlpCandidates.push_back(ovlp);
				}
			}
		}

		//computing divergence
		for (auto& ovlp : ovlpCandidates)
		{
			int rate = std::max(1, ovlp.curRange() / 5000);

			curOvlpKmers.clear();
			extOvlpKmers.clear();
			intersect.clear();
			
			for (auto extKmerPos : IterKmers(_seqContainer.getSeq(extId),
											 ovlp.extBegin, 
											 ovlp.extRange() - kmerSize))
			{
				if (extKmerPos.position % rate == 0)
				{
					extOvlpKmers.push_back(extKmerPos.kmer);
				}
			}
			for (int i = ovlp.curBegin; i < ovlp.curEnd - kmerSize; i += rate)
			{
				curOvlpKmers.push_back(curKmers[i]);
			}
			std::sort(curOvlpKmers.begin(), curOvlpKmers.end());
			std::sort(extOvlpKmers.begin(), extOvlpKmers.end());
			std::set_intersection(curOvlpKmers.begin(), curOvlpKmers.end(),
								  extOvlpKmers.begin(), extOvlpKmers.end(),
								  std::back_inserter(intersect));

			size_t uniqueCur = 1;
			for (size_t i = 1; i < curOvlpKmers.size(); ++i)
			{
				if (curOvlpKmers[i - 1] != curOvlpKmers[i]) ++uniqueCur;
			}
			size_t uniqueInt = 1;
			for (size_t i = 1; i < intersect.size(); ++i)
			{
				if (intersect[i - 1] != intersect[i]) ++uniqueInt;
			}
			float kmerDiv = (float)uniqueInt * rate / uniqueCur;
			float seqDiv = std::log(1 / kmerDiv) / kmerSize;

			//fout << ovlp.curRange() << " " << seqDiv << std::endl;

			ovlp.seqDivergence = seqDiv;
			if (seqDiv < _maxDivergence) detectedOverlaps.push_back(ovlp);
		}

		//totalBackLoop += double(clock() - dpEnd) / CLOCKS_PER_SEC;
	}

	/*Logger::get().debug() << "---------";
	Logger::get().debug() << " " << vecMatches.size() << " "
		<< seqMatches.size() << " " << uniqueCandidates
		<< " " << detectedOverlaps.size();
	Logger::get().debug() << "Repeat intensity: "
		<< vecMatches.size() / (solidPos.size() + 1);
	Logger::get().debug() << "hash: " << totalHashTime << " k-mer: "
		<< totalKmerTime << " dp: " << totalDpTime
		<< " dpLoop: " << totalDpLoop << " backLoop: " << totalBackLoop;*/


	//clock_t ff = clock();
	//double es = double(ff - end) / CLOCKS_PER_SEC;
	//totalDpTime += es;

	return detectedOverlaps;
}

bool OverlapContainer::hasSelfOverlaps(FastaRecord::Id readId)
{
	this->lazySeqOverlaps(readId);
	if (!readId.strand()) readId = readId.rc();
	return _overlapIndex.find(readId).suggestChimeric;
}


std::vector<OverlapRange> 
	OverlapContainer::quickSeqOverlaps(FastaRecord::Id readId) const
{
	bool suggestChimeric;
	const FastaRecord& record = _queryContainer.getRecord(readId);
	return _ovlpDetect.getSeqOverlaps(record, suggestChimeric);
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
	const FastaRecord& record = _queryContainer.getRecord(readId);
	auto overlaps = _ovlpDetect.getSeqOverlaps(record, suggestChimeric);

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

void OverlapContainer::ensureTransitivity(bool onlyMaxExt)
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
		std::vector<OverlapRange> ovlpsToAdd;
		for (auto& curOvlp : curOvlps)
		{
			auto& extOvlps = this->unsafeSeqOverlaps(curOvlp.extId);

			if (onlyMaxExt)
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
					ovlpsToAdd.push_back(curOvlp.reverse());
				}
			}
			else
			{
				ovlpsToAdd.push_back(curOvlp.reverse());
			}
		}
		for (auto& ovlp : ovlpsToAdd)
		{
			this->unsafeSeqOverlaps(ovlp.curId).push_back(ovlp);
		}
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
		this->lazySeqOverlaps(seqId);	//automatically stores overlaps
	};
	processInParallel(allQueries, indexUpdate, 
					  Parameters::get().numThreads, true);
	this->ensureTransitivity(false);

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
		IndexVecWrapper wrapper = _overlapIndex.find(normId);
		return seqId.strand() ? *wrapper.fwdOverlaps : 
								*wrapper.revOverlaps;
}

//TODO: potentially might become non-symmetric after filtering
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
		auto& overlaps = this->unsafeSeqOverlaps(seqId);
		
		SetVec<OverlapRange*> overlapSets;
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
		auto clusters = groupBySet(overlapSets);
		std::vector<OverlapRange> newOvlps;
		for (auto& cluster : clusters)
		{
			OverlapRange* maxOvlp = nullptr;
			for (auto& ovlp : cluster.second)
			{
				if (!maxOvlp || ovlp->score > maxOvlp->score)
				{
					maxOvlp = ovlp;
				}
			}
			newOvlps.push_back(*maxOvlp);
		}
		overlaps = std::move(newOvlps);

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
