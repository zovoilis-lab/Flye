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
	//static const float OVLP_DIVERGENCE = Config::get("overlap_divergence_rate");
	if (ovlp.curRange() < _minOverlap || 
		ovlp.extRange() < _minOverlap) 
	{
		return false;
	}

	/*float lengthDiff = abs(ovlp.curRange() - ovlp.extRange());
	float meanLength = (ovlp.curRange() + ovlp.extRange()) / 2.0f;
	if (lengthDiff > meanLength * OVLP_DIVERGENCE)
	{
		return false;
	}*/

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

	static float totalDpTime = 0;
	static float totalKmerTime = 0;
	static float totalHashTime = 0;
	static float totalDpLoop = 0;
	static float totalBackLoop = 0;
	clock_t begin = clock();

	outSuggestChimeric = false;
	int32_t curLen = fastaRec.sequence.length();

	thread_local static std::vector<uint16_t> seqHitCount;
	if (seqHitCount.size() != _seqContainer.getMaxSeqId())
	{
		seqHitCount = std::vector<uint16_t>(_seqContainer.getMaxSeqId(), 0);
	}
	thread_local static std::vector<KmerMatch> vecMatches;

	std::vector<int32_t> solidPos;
	for (auto curKmerPos : IterKmers(fastaRec.sequence))
	{
		if (!_vertexIndex.isSolid(curKmerPos.kmer)) continue;

		solidPos.push_back(curKmerPos.position);
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

	size_t maxExtSeq = seqHitCount.size() - 1;
	size_t selectedExtSeqs = 0;
	if (_maxCurOverlaps > 0)
	{
		for (size_t i = 0; i < seqHitCount.size(); ++i)
		{
			if (seqHitCount[i] >= MIN_KMER_SURV_RATE * _minOverlap)
			{
				++selectedExtSeqs;
				if ((int)selectedExtSeqs > 2 * _maxCurOverlaps)
				{
					maxExtSeq = i;
					break;
				}
			}
		}
	}

	thread_local static cuckoohash_map<FastaRecord::Id, 
									   MatchVecWrapper> seqMatches;
	for (auto& match : vecMatches)
	{
		if (match.extId.rawId() > maxExtSeq) continue;
		if (seqHitCount[match.extId.rawId()] < 
			MIN_KMER_SURV_RATE * _minOverlap) continue;

		seqMatches.upsert(match.extId, 
			  [&match](MatchVecWrapper& v)
			  {
				  v->push_back(match);
			  }, 
			  match, seqHitCount[match.extId.rawId()]);	//if key was not found
	}
	
  	clock_t end = clock();
  	double elapsed_secs = double(end - hashTime) / CLOCKS_PER_SEC;
	totalKmerTime += elapsed_secs;

	std::vector<OverlapRange> detectedOverlaps;
	int uniqueCandidates = 0;
	for (auto& seqVec : seqMatches.lock_table())
	{
		std::vector<KmerMatch>& matchesList = *seqVec.second;
		int32_t extLen = _seqContainer.seqLen(seqVec.first);

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
		++uniqueCandidates;
		

		//chain matiching positions with DP
		std::vector<int32_t> scoreTable(matchesList.size(), 0);
		std::vector<int32_t> backtrackTable(matchesList.size(), -1);

		bool extSorted = extLen > curLen;
		if (extSorted)
		{
			std::sort(matchesList.begin(), matchesList.end(),
					  [](const KmerMatch& k1, const KmerMatch& k2)
					  {return k1.extPos < k2.extPos;});
		}

  		clock_t dpBegin = clock();
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
		clock_t dpEnd = clock();
  		totalDpLoop += double(dpEnd - dpBegin) / CLOCKS_PER_SEC;

		//backtracking
		std::vector<OverlapRange> extOverlaps;
		std::vector<int32_t> shifts;
		shifts.reserve(1024);
		std::vector<std::pair<int32_t, int32_t>> kmerMatches;
		kmerMatches.reserve(1024);
		for (int32_t chainStart = backtrackTable.size() - 1; 
			 chainStart > 0; --chainStart)
		{
			if (backtrackTable[chainStart] == -1) continue;

			int32_t pos = chainStart;
			KmerMatch lastMatch = matchesList[pos];
			KmerMatch firstMatch = lastMatch;

			int chainLength = 0;
			shifts.clear();
			kmerMatches.clear();
			//int totalMatch = kmerSize;
			while (pos != -1)
			{
				firstMatch = matchesList[pos];
				shifts.push_back(matchesList[pos].curPos - 
								 matchesList[pos].extPos);
				++chainLength;

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
				int numCurKmers = 0;
				for (auto pos : solidPos)
				{
					if (ovlp.curBegin <= pos && pos <= ovlp.curEnd) ++numCurKmers;
				}
				float kmerMatch = std::min(1.0f, (float)chainLength * 
									_vertexIndex.getSampleRate() / numCurKmers);
				float divSet = std::log(1 / kmerMatch) / kmerSize;
				if (divSet < _maxDivergence)
				{
					std::reverse(kmerMatches.begin(), kmerMatches.end());
					ovlp.kmerMatches = kmerMatches;
					ovlp.leftShift = median(shifts);
					ovlp.rightShift = extLen - curLen + ovlp.leftShift;

					extOverlaps.push_back(ovlp);
				}
				//fout << divSet << std::endl;
				//Logger::get().debug() << ovlp.curRange() << " " <<
				//	" " << kmerMatch << " " << divSet;
			}
		}
		
		//selecting the best
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
			if (maxOvlp) detectedOverlaps.push_back(*maxOvlp);
		}
		else
		{
			//sort by decreasing score
			//std::sort(extOverlaps.begin(), extOverlaps.end(),
			//		  [](const OverlapRange& r1, const OverlapRange& r2)
			//		  {return r1.score > r2.score;});
			
			for (auto& ovlp : extOverlaps)
			{
				detectedOverlaps.push_back(ovlp);
				/*bool isContained = false;
				for (auto& prim : detectedOverlaps)
				{
					if (ovlp.containedBy(prim))
					{
						isContained = true;
						break;
					}
				}
				if (!isContained)
				{
					detectedOverlaps.push_back(ovlp);
				}*/
			}
		}
		totalBackLoop += double(clock() - dpEnd) / CLOCKS_PER_SEC;

		if (_maxCurOverlaps > 0 &&
			detectedOverlaps.size() > (size_t)_maxCurOverlaps) break;
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

	seqHitCount.assign(seqHitCount.size(), 0);
	vecMatches.clear();
	seqMatches.clear();

	clock_t ff = clock();
	double es = double(ff - end) / CLOCKS_PER_SEC;
	totalDpTime += es;

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
