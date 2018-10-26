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
#include <chrono>
#include <ctime>
#include <cstring>
#include <iomanip>

#include "ksw2.h"

#include "overlap.h"
#include "../common/config.h"
#include "../common/utils.h"
#include "../common/parallel.h"
#include "../common/disjoint_set.h"

namespace
{
	float kswAlign(const DnaSequence& seqT, const DnaSequence seqQ,
				   int matchScore, int misScore, int gapOpen, int gapExtend,
				   bool showAlignment)
	{
		static const int32_t MAX_JUMP = Config::get("maximum_jump");

		std::vector<uint8_t> tseq(seqT.length());
		for (size_t i = 0; i < (size_t)seqT.length(); ++i)
		{
			tseq[i] = seqT.atRaw(i);
		}
		std::vector<uint8_t> qseq(seqQ.length());
		for (size_t i = 0; i < seqQ.length(); ++i)
		{
			qseq[i] = seqQ.atRaw(i);
		}

		int seqDiff = abs((int)tseq.size() - (int)qseq.size());
		int bandWidth = seqDiff + MAX_JUMP * 2;
		//substitution matrix
		int8_t a = matchScore;
		int8_t b = misScore < 0 ? misScore : -misScore; // a > 0 and b < 0
		int8_t subsMat[] = {a, b, b, b, 0, 
						  	b, a, b, b, 0, 
						  	b, b, a, b, 0, 
						  	b, b, b, a, 0, 
						  	0, 0, 0, 0, 0};

		ksw_extz_t ez;
		memset(&ez, 0, sizeof(ksw_extz_t));
		const int NUM_NUCL = 5;
		const int Z_DROP = -1;
		const int FLAG = KSW_EZ_APPROX_MAX | KSW_EZ_APPROX_DROP;
		const int END_BONUS = 0;
		//ksw_extf2_sse(0, qseq.size(), &qseq[0], tseq.size(), &tseq[0], matchScore,
		//		 	  misScore, gapOpen, bandWidth, Z_DROP, &ez);
		ksw_extz2_sse(0, qseq.size(), &qseq[0], tseq.size(), &tseq[0], NUM_NUCL,
				 	  subsMat, gapOpen, gapExtend, bandWidth, Z_DROP, 
					  END_BONUS, FLAG, &ez);
		
		//decode CIGAR
		
		std::string strQ;
		std::string strT;
		std::string alnQry;
		std::string alnTrg;
		if (showAlignment)
		{
			for (auto x : qseq) strQ += "ACGT"[x];
			for (auto x : tseq) strT += "ACGT"[x];
		}

		size_t posQry = 0;
		size_t posTrg = 0;
        int matches = 0;
		int alnLength = 0;
		for (size_t i = 0; i < (size_t)ez.n_cigar; ++i)
		{
			int size = ez.cigar[i] >> 4;
			char op = "MID"[ez.cigar[i] & 0xf];
			alnLength += size;

        	if (op == 'M')
			{
				for (size_t i = 0; i < (size_t)size; ++i)
				{
					if (tseq[posTrg + i] == qseq[posQry + i]) ++matches;
				}
				if (showAlignment)
				{
					alnQry += strQ.substr(posQry, size);
					alnTrg += strT.substr(posTrg, size);
				}
				posQry += size;
				posTrg += size;
			}
            else if (op == 'I')
			{
				if (showAlignment)
				{
					alnQry += strQ.substr(posQry, size);
					alnTrg += std::string(size, '-');
				}
                posQry += size;
			}
            else //D
			{
				if (showAlignment)
				{
					alnQry += std::string(size, '-');
					alnTrg += strT.substr(posTrg, size);
				}
				posTrg += size;
			}
		}
        //float errRate = 1 - float(matches) / alnLength;
        float errRate = 1 - float(matches) / std::max(tseq.size(), qseq.size());
		free(ez.cigar);

		if (showAlignment)
		{
			for (size_t chunk = 0; chunk <= alnQry.size() / 100; ++chunk)
			{
				for (size_t i = chunk * 100; 
					 i < std::min((chunk + 1) * 100, alnQry.size()); ++i)
				{
					std::cout << alnQry[i];
				}
				std::cout << "\n";
				for (size_t i = chunk * 100; 
					 i < std::min((chunk + 1) * 100, alnQry.size()); ++i)
				{
					std::cout << alnTrg[i];
				}
				//std::cout << alnQry;
				std::cout << "\n\n";
				//std::cout << std::endl;
			}
			std::cout << "\n----------------\n" << std::endl;
		}

		return errRate;
	}
}

//Check if it is a proper overlap
bool OverlapDetector::overlapTest(const OverlapRange& ovlp,
								  bool& outSuggestChimeric) const
{
	if (ovlp.curRange() < _minOverlap || 
		ovlp.extRange() < _minOverlap) 
	{
		return false;
	}

	//filter overlaps that way to divergent in length.
	//theoretically, they should not pass sequence divergence filter,
	//but just in case
	static const float OVLP_DIV = 0.5;
	float lengthDiff = abs(ovlp.curRange() - ovlp.extRange());
	if (lengthDiff > OVLP_DIV * std::min(ovlp.curRange(), ovlp.extRange()))
	{
		return false;
	}

	if (ovlp.curId == ovlp.extId.rc()) 
	{
		int32_t projEnd = ovlp.extLen - ovlp.extEnd - 1;
		if (abs(ovlp.curEnd - projEnd) < _maxJump)
		{
			outSuggestChimeric = true;
		}
	}
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

	static const char LogTable256[256] = {
		#define LT(n) n, n, n, n, n, n, n, n, n, n, n, n, n, n, n, n
		-1, 0, 1, 1, 2, 2, 2, 2, 3, 3, 3, 3, 3, 3, 3, 3,
		LT(4), LT(5), LT(5), LT(6), LT(6), LT(6), LT(6),
		LT(7), LT(7), LT(7), LT(7), LT(7), LT(7), LT(7), LT(7)
	};

	static inline int ilog2_32(uint32_t v)
	{
		uint32_t t, tt;
		if ((tt = v>>16)) return (t = tt>>8) ? 24 + LogTable256[t] : 16 + LogTable256[tt];
		return (t = v>>8) ? 8 + LogTable256[t] : LogTable256[v];
	}
}

//This implementation was inspired by Heng Li's minimap2 paper
//might be used in parallel
std::vector<OverlapRange> 
OverlapDetector::getSeqOverlaps(const FastaRecord& fastaRec, 
								bool& outSuggestChimeric,
								OvlpDivStats& divStats) const
{
	//const int MAX_LOOK_BACK = 50;
	const int kmerSize = Parameters::get().kmerSize;
	//const float minKmerSruvivalRate = std::exp(-_maxDivergence * kmerSize);
	const float minKmerSruvivalRate = 0.01;

	//static std::ofstream fout("../kmers.txt");

	//static float totalDpTime = 0;
	//static float totalKmerTime = 0;
	//static float totalHashTime = 0;
	//static float totalDpLoop = 0;
	//static float totalBackLoop = 0;
	//clock_t begin = clock();

	//cache memory-intensive containers as
	//many parallel memory allocations slow us down significantly
	typedef uint32_t SeqCountType;
	static thread_local 
		cuckoohash_map<FastaRecord::Id, SeqCountType> seqHitCount;
	auto lockedHitCount = seqHitCount.lock_table();
	thread_local static std::vector<KmerMatch> vecMatches;
	static thread_local std::vector<KmerMatch> matchesList;
	static thread_local std::vector<int32_t> scoreTable;
	static thread_local std::vector<int32_t> backtrackTable;

	//although once in a while shrink allocated memory size
	static thread_local auto prevClenup = std::chrono::system_clock::now();
	if ((std::chrono::system_clock::now() - prevClenup) > 
		std::chrono::seconds(60))
	{
		prevClenup = std::chrono::system_clock::now();
		vecMatches.shrink_to_fit();
		matchesList.shrink_to_fit();
		scoreTable.shrink_to_fit();
		backtrackTable.shrink_to_fit();
	}
	lockedHitCount.clear();
	vecMatches.clear();

	outSuggestChimeric = false;
	int32_t curLen = fastaRec.sequence.length();
	std::vector<int32_t> curFilteredPos;

	//count kmer hits
	for (auto curKmerPos : IterKmers(fastaRec.sequence))
	{
		if (_vertexIndex.isRepetitive(curKmerPos.kmer))
		{
			curFilteredPos.push_back(curKmerPos.position);
		}
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
					std::numeric_limits<SeqCountType>::max())
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
		std::vector<std::pair<FastaRecord::Id, SeqCountType>> topSeqs;
		topSeqs.reserve(lockedHitCount.size());
		for (auto seqCount : lockedHitCount)
		{
			if (seqCount.second >= minKmerSruvivalRate * _minOverlap)
			{
				topSeqs.emplace_back(seqCount.first, seqCount.second);
			}
		}
		std::sort(topSeqs.begin(), topSeqs.end(),
				  [](const std::pair<FastaRecord::Id, SeqCountType> p1, 
					 const std::pair<FastaRecord::Id, SeqCountType>& p2) 
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
				minKmerSruvivalRate * _minOverlap)
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
	
	const int STAT_WND = 10000;
	std::vector<OverlapRange> divStatWindows(curLen / STAT_WND + 1);

	std::vector<OverlapRange> detectedOverlaps;
	size_t extRangeBegin = 0;
	size_t extRangeEnd = 0;
	while(extRangeEnd < vecMatches.size())
	{
		extRangeBegin = extRangeEnd;
		while (extRangeEnd < vecMatches.size() &&
			   vecMatches[extRangeBegin].extId == 
			   vecMatches[extRangeEnd].extId) ++extRangeEnd;
		matchesList.assign(vecMatches.begin() + extRangeBegin,
						   vecMatches.begin() + extRangeEnd);
		FastaRecord::Id extId = vecMatches[extRangeBegin].extId;
		int32_t extLen = _seqContainer.seqLen(extId);

		//pre-filtering
		int32_t minCur = matchesList.front().curPos;
		int32_t maxCur = matchesList.back().curPos;
		int32_t minExt = std::numeric_limits<int32_t>::max();
		int32_t maxExt = std::numeric_limits<int32_t>::min();
		for (auto& match : matchesList)
		{
			minExt = std::min(minExt, match.extPos);
			maxExt = std::max(maxExt, match.extPos);
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
			//int32_t noImprovement = 0;

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
					//convex gap seems useless..
					//int32_t gapCost = jumpDiv ? 
					//		kmerSize * jumpDiv + ilog2_32(jumpDiv) : 0;
					int32_t gapCost = (jumpDiv > 50) ? 2 * jumpDiv : 0;
					int32_t nextScore = scoreTable[j] + matchScore - gapCost;
					if (nextScore > maxScore)
					{
						maxScore = nextScore;
						maxId = j;
						//noImprovement = 0;

						if (jumpDiv == 0 && curNext - curPrev < kmerSize) break;
					}
					/*else
					{
						if (++noImprovement > MAX_LOOK_BACK) break;
					}*/
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
			//always start from a local maximum
			while (chainStart > 0 && 
				   scoreTable[chainStart - 1] > scoreTable[chainStart]) 
			{
				--chainStart;
			}

			if (backtrackTable[chainStart] == -1) continue;

			int32_t pos = chainStart;
			KmerMatch lastMatch = matchesList[pos];
			KmerMatch firstMatch = lastMatch;

			int32_t chainLength = 0;
			int32_t chainEnd = 0;
			shifts.clear();
			kmerMatches.clear();
			//int32_t totalMatch = kmerSize;
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
				if (newPos == -1)
				{
					chainEnd = pos;
				}
				pos = newPos;
			}

			OverlapRange ovlp(fastaRec.id, matchesList.front().extId,
							  firstMatch.curPos, firstMatch.extPos,
							  curLen, extLen);
			ovlp.curEnd = lastMatch.curPos + kmerSize - 1;
			ovlp.extEnd = lastMatch.extPos + kmerSize - 1;
			ovlp.score = scoreTable[chainStart] - scoreTable[chainEnd] + 
						 kmerSize - 1;

			if (this->overlapTest(ovlp, outSuggestChimeric))
			{
				if (_keepAlignment)
				{
					kmerMatches.emplace_back(ovlp.curBegin, ovlp.extBegin);
					std::reverse(kmerMatches.begin(), kmerMatches.end());
					kmerMatches.emplace_back(ovlp.curEnd, ovlp.extEnd);
					ovlp.kmerMatches = kmerMatches;
				}
				ovlp.leftShift = median(shifts);
				ovlp.rightShift = extLen - curLen + ovlp.leftShift;

				int32_t filteredPositions = 0;
				for (auto pos : curFilteredPos)
				{
					if (pos < ovlp.curBegin) continue;
					if (pos > ovlp.curEnd) break;
					++filteredPositions;
				}

				float matchRate = (float)chainLength / std::max(ovlp.curRange(), 
																ovlp.extRange());
				//float repeatRate = (float)filteredPositions / ovlp.curRange();
				ovlp.seqDivergence = std::log(1 / matchRate) / kmerSize;

				if(_nuclAlignment)
				{
					ovlp.seqDivergence = 
						kswAlign(fastaRec.sequence.substr(ovlp.curBegin, 
														 ovlp.curRange()),
								 _seqContainer.getSeq(extId).substr(ovlp.extBegin, 
								 									ovlp.extRange()),
								 /*match*/ 1, /*mm*/ -2, /*gap open*/ 2, 
								 /*gap ext*/ 1, false);
				}

				float divThreshold = _maxDivergence;
				if (ovlp.curBegin < _maxJump || curLen - ovlp.curEnd < _maxJump ||
					ovlp.extBegin < _maxJump || extLen - ovlp.extEnd < _maxJump)
				{
					divThreshold += _badEndAdjustment;
				}
				if (ovlp.seqDivergence < divThreshold)
				{
					extOverlaps.push_back(ovlp);
				}

				//collecting overlap statistics
				size_t wnd = ovlp.curBegin / STAT_WND;
				if (ovlp.curRange() > divStatWindows[wnd].curRange())
				{
					divStatWindows[wnd] = ovlp;
				}

				//benchmarking divergence
				/*float alnDiff = kswAlign(fastaRec.sequence
											.substr(ovlp.curBegin, ovlp.curRange()),
										 _seqContainer.getSeq(extId)
											.substr(ovlp.extBegin, ovlp.extRange()),
										 1, -2, 2, 1, false);
				fout << alnDiff << " " << ovlp.seqDivergence << std::endl;*/
				/*if (0.15 > alnDiff && ovlp.seqDivergence > 0.20)
				{
					kswAlign(fastaRec.sequence
								.substr(ovlp.curBegin, ovlp.curRange()),
							 _seqContainer.getSeq(extId)
								.substr(ovlp.extBegin, ovlp.extRange()),
							 1, -2, 2, 1, true);
					std::cout << alnDiff << " " << ovlp.seqDivergence << 
						" " << (float)filteredPositions / ovlp.curRange() << std::endl;
				}*/
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
			std::vector<OverlapRange> primOverlaps;
			//sort by decreasing score
			std::sort(extOverlaps.begin(), extOverlaps.end(),
					  [](const OverlapRange& r1, const OverlapRange& r2)
					  {return r1.score > r2.score;});
			
			for (auto& ovlp : extOverlaps)
			{
				bool isContained = false;
				for (auto& prim : primOverlaps)
				{
					if (ovlp.containedBy(prim) && prim.score > ovlp.score)
					{
						isContained = true;
						break;
					}
				}
				if (!isContained)
				{
					primOverlaps.push_back(ovlp);
				}
			}
			for (auto& ovlp : primOverlaps)
			{
				detectedOverlaps.push_back(ovlp);
			}
		}
		//totalBackLoop += double(clock() - dpEnd) / CLOCKS_PER_SEC;
	}

	for (auto ovlp : divStatWindows)
	{
		if (ovlp.curRange() > 0)
		{
			divStats.add(ovlp.seqDivergence);
		}
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
	OverlapContainer::quickSeqOverlaps(FastaRecord::Id readId)
{
	bool suggestChimeric;
	const FastaRecord& record = _queryContainer.getRecord(readId);
	return _ovlpDetect.getSeqOverlaps(record, suggestChimeric, 
									  _divergenceStats);
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
	auto overlaps = _ovlpDetect.getSeqOverlaps(record, suggestChimeric, 
											   _divergenceStats);
	overlaps.shrink_to_fit();

	std::vector<OverlapRange> revOverlaps;
	revOverlaps.reserve(overlaps.size());
	for (auto& ovlp : overlaps) revOverlaps.push_back(ovlp.complement());

	_overlapIndex.update_fn(readId,
		[&wrapper, &overlaps, &revOverlaps, &suggestChimeric, &flipped, this]
		(IndexVecWrapper& val)
		{
			if (!val.cached)
			{
				_indexSize += overlaps.size();
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
	static const int MAX_ENDS_DIFF = Parameters::get().kmerSize;

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


void OverlapContainer::overlapDivergenceStats()
{
	std::vector<float> ovlpDivergence(_divergenceStats.divVec.begin(),
									  _divergenceStats.divVec.begin() + 
									  		_divergenceStats.vecSize);
	const int HIST_LENGTH = 100;
	const int HIST_HEIGHT = 20;
	const float HIST_MIN = 0;
	const float HIST_MAX = 0.5;
	const float mult = HIST_LENGTH / (HIST_MAX * 100);
	std::vector<int> histogram(HIST_LENGTH, 0);
	for (float d : ovlpDivergence)
	{
		if (HIST_MIN <= d && d < HIST_MAX) 
		{
			++histogram[int(d * mult * 100)];
		}
	}
	int histMax = 1;
	for (int freq : histogram) histMax = std::max(histMax, freq);

	std::string histString = "\n";
	for (int height = HIST_HEIGHT - 1; height >= 0; --height)
	{
		histString += "    |";
		for (int i = 0; i < HIST_LENGTH; ++i)
		{
			if ((float)histogram[i] / histMax > (float)height / HIST_HEIGHT)
			{
				histString += '*';
			}
			else
			{
				histString += ' ';
			}
		}
		histString += '\n';
	}
	histString += "    " + std::string(HIST_LENGTH,  '-') + "\n";
	std::string footer(HIST_LENGTH, ' ');
	for (int i = 0; i < 10; ++i)
	{
		size_t startPos = i * HIST_LENGTH / 10;
		auto s = std::to_string(i * 5) + "%";
		for (size_t j = 0; j < s.size(); ++j) footer[j + startPos] = s[j];
	}
	histString += "    " + footer + "\n";

	Logger::get().info() << "Median overlap divergence: " 
		<< quantile(ovlpDivergence, 50); 
	Logger::get().debug() << "Sequence divergence distribution: \n" << histString
		<< "\n    Q25 = " << std::setprecision(2)
		<< quantile(ovlpDivergence, 25) << ", Q50 = " 
		<< quantile(ovlpDivergence, 50)
		<< ", Q75 = " << quantile(ovlpDivergence, 75) << "\n"
		<< std::setprecision(6);
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
