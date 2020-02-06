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

#define HAVE_KALLOC
#include "kalloc.h"
#include "ksw2.h"
#undef HAVE_KALLOC

#include "overlap.h"
#include "../common/config.h"
#include "../common/utils.h"
#include "../common/parallel.h"
#include "../common/disjoint_set.h"
#include "../common/bfcontainer.h"


namespace
{
	using namespace std::chrono;
	struct ThreadMemPool
	{
		ThreadMemPool():
			prevCleanup(system_clock::now() + seconds(rand() % 60))
		{
		   memPool = km_init();
		}
		~ThreadMemPool()
		{
		   km_destroy(memPool);
		}
		void cleanIter()
		{
			if ((system_clock::now() - prevCleanup) > seconds(60))
			{
				km_destroy(memPool);
               	memPool = km_init();
				prevCleanup = system_clock::now();
			}
		}

		time_point<system_clock> prevCleanup;
		void* memPool;
    };

	struct CigOp
	{
		char op;
		int len;
	};
	float kswAlign(const DnaSequence& trgSeq, size_t trgBegin, size_t trgLen,
				   const DnaSequence& qrySeq, size_t qryBegin, size_t qryLen,
				   int matchScore, int misScore, int gapOpen, int gapExtend,
				   std::vector<CigOp>& cigarOut)
	{
		static const int32_t MAX_JUMP = Config::get("maximum_jump");
		const int KMER_SIZE = Parameters::get().kmerSize;

		thread_local ThreadMemPool buf;
		thread_local std::vector<uint8_t> trgByte;
		thread_local std::vector<uint8_t> qryByte;
		buf.cleanIter();
		trgByte.assign(trgLen, 0);
		qryByte.assign(qryLen, 0);

		for (size_t i = 0; i < trgLen; ++i)
		{
			trgByte[i] = trgSeq.atRaw(i + trgBegin);
		}
		for (size_t i = 0; i < qryLen; ++i)
		{
			qryByte[i] = qrySeq.atRaw(i + qryBegin);
		}

		int seqDiff = abs((int)trgByte.size() - (int)qryByte.size());
		int bandWidth = seqDiff + MAX_JUMP;

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
		ksw_extz2_sse(buf.memPool, qryByte.size(), &qryByte[0], 
					  trgByte.size(), &trgByte[0], NUM_NUCL,
				 	  subsMat, gapOpen, gapExtend, bandWidth, Z_DROP, 
					  END_BONUS, FLAG, &ez);
		
        int numMatches = 0;
		int numMiss = 0;
		int numIndels = 0;

		cigarOut.clear();
		cigarOut.reserve((size_t)ez.n_cigar);

		//decode cigar
		size_t posQry = 0;
		size_t posTrg = 0;
		for (size_t i = 0; i < (size_t)ez.n_cigar; ++i)
		{
			int size = ez.cigar[i] >> 4;
			char op = "MID"[ez.cigar[i] & 0xf];
			//alnLength += size;

        	if (op == 'M')
			{
				for (size_t i = 0; i < (size_t)size; ++i)
				{
					char match = "X="[size_t(trgByte[posTrg + i] == 
											 qryByte[posQry + i])];
					if (i == 0 || (match != cigarOut.back().op))
					{
						cigarOut.push_back({match, 1});
					}
					else
					{
						++cigarOut.back().len;
					}
					numMatches += int(match == '=');
					numMiss += int(match == 'X');
				}
				posQry += size;
				posTrg += size;
			}
            else if (op == 'I')
			{
				cigarOut.push_back({'I', size});
                posQry += size;
				numIndels += std::min(size, KMER_SIZE);
			}
            else //D
			{
				cigarOut.push_back({'D', size});
				posTrg += size;
				numIndels += std::min(size, KMER_SIZE);
			}
		}
        float errRate = 1 - float(numMatches) / (numMatches + numMiss + numIndels);

		kfree(buf.memPool, ez.cigar);
		return errRate;
	}

	float getAlignmentIdy(const OverlapRange& ovlp,
						  const DnaSequence& trgSeq,
						  const DnaSequence& qrySeq,
						  bool showAlignment)
	{
		std::vector<CigOp> decodedCigar;
		float errRate = kswAlign(trgSeq, ovlp.curBegin, ovlp.curRange(),
								 qrySeq, ovlp.extBegin, ovlp.extRange(),
								 /*match*/ 1, /*mm*/ -2, /*gap open*/ 2, 
								 /*gap ext*/ 1, decodedCigar);

		//visualize alignents if needed
		if (showAlignment)
		{
			std::vector<uint8_t> trgByte;
			std::vector<uint8_t> qryByte;
			trgByte.assign(ovlp.curRange(), 0);
			qryByte.assign(ovlp.extRange(), 0);
			for (size_t i = 0; i < (size_t)ovlp.curRange(); ++i)
			{
				trgByte[i] = trgSeq.atRaw(i + ovlp.curBegin);
			}
			for (size_t i = 0; i < (size_t)ovlp.extRange(); ++i)
			{
				qryByte[i] = qrySeq.atRaw(i + ovlp.extBegin);
			}

			std::string strQ;
			std::string strT;
			std::string alnQry;
			std::string alnTrg;

			for (auto x : qryByte) strQ += "ACGT"[x];
			for (auto x : trgByte) strT += "ACGT"[x];

			size_t posQry = 0;
			size_t posTrg = 0;
			for (auto& op : decodedCigar)
			{
				if (op.op == '=' || op.op == 'X')
				{
					alnQry += strQ.substr(posQry, op.len);
					alnTrg += strT.substr(posTrg, op.len);
					posQry += op.len;
					posTrg += op.len;
				}
				else if (op.op == 'I')
				{
					alnQry += strQ.substr(posQry, op.len);
					alnTrg += std::string(op.len, '-');
                	posQry += op.len;
				}
				else
				{
					alnQry += std::string(op.len, '-');
					alnTrg += strT.substr(posTrg, op.len);
					posTrg += op.len;
				}
			}

			const int WIDTH = 100;
			for (size_t chunk = 0; chunk <= alnQry.size() / WIDTH; ++chunk)
			{
				for (size_t i = chunk * WIDTH; 
					 i < std::min((chunk + 1) * WIDTH, alnQry.size()); ++i)
				{
					std::cout << alnQry[i];
				}
				std::cout << "\n";
				for (size_t i = chunk * WIDTH; 
					 i < std::min((chunk + 1) * WIDTH, alnQry.size()); ++i)
				{
					std::cout << alnTrg[i];
				}
				std::cout << "\n\n";
			}
		}

		return errRate;
	}

	std::vector<OverlapRange> 
		checkIdyAndTrim(OverlapRange& ovlp, const DnaSequence& trgSeq,
						const DnaSequence& qrySeq, float maxDivergence,
						int32_t minOverlap, bool showAlignment)
	{
		std::vector<CigOp> decodedCigar;
		float errRate = kswAlign(trgSeq, ovlp.curBegin, ovlp.curRange(),
								 qrySeq, ovlp.extBegin, ovlp.extRange(),
								 /*match*/ 1, /*mm*/ -2, /*gap open*/ 2, 
								 /*gap ext*/ 1, decodedCigar);
		//the original alignment is already passing the threshold
		ovlp.seqDivergence = errRate;
		if (errRate < maxDivergence) 
		{
			return {ovlp};
		}

		//if not, then...

		//precomputing some stuff
		std::vector<int> sumMatches;
		sumMatches.reserve(decodedCigar.size() + 1);
		sumMatches.push_back(0);
		std::vector<int> sumLength;
		sumLength.reserve(decodedCigar.size() + 1);
		sumLength.push_back(0);
		for (auto op : decodedCigar)
		{
			sumLength.push_back(sumLength.back() + op.len);
			int match = (op.op == '=') ? op.len : 0;
			sumMatches.push_back(sumMatches.back() + match);
		}

		std::vector<std::pair<int, int>> goodIntervals;
		const float EPS = 0.005;

		for (int intLen = (int)decodedCigar.size(); intLen > 0; --intLen)
		{
			for (int intStart = 0; 
				intStart < (int)decodedCigar.size() - intLen + 1; ++intStart)
			{
				int i = intStart;
				int j = intStart + intLen - 1;
				int rangeLen = sumLength[j + 1] - sumLength[i];
				int rangeMatch = sumMatches[j + 1] - sumMatches[i];
				if (1.0f - float(rangeMatch) / rangeLen < maxDivergence - EPS)
				{
					if (j - i >= 1) goodIntervals.emplace_back(i, j);
				}
			}
		}

		//select non-intersecting set
		std::vector<std::pair<int, int>> nonIntersecting;
		for (auto& interval : goodIntervals)
		{
			bool intersects = false;
			for (auto& otherInt : nonIntersecting)
			{
				int ovl = std::min(interval.second, otherInt.second) - 
						std::max(interval.first, otherInt.first);
				if (ovl > 0) 
				{
					intersects = true;
					break;
				}
			}
			if (!intersects) nonIntersecting.push_back(interval);
		}

		//now, for each interesting interval ajust boundaries and check the
		//actual sequence length
		std::vector<OverlapRange> trimmedAlignments;
		for (auto intCand : nonIntersecting)
		{
			while (intCand.first < (int)decodedCigar.size() && 
				   decodedCigar[intCand.first].op != '=') ++intCand.first;
			while (intCand.second > 0 && 
				   decodedCigar[intCand.second].op != '=') ++intCand.second;
			if (intCand.second - intCand.first < 1) continue;

			int rangeLen = sumLength[intCand.second + 1] - sumLength[intCand.first];
			int rangeMatch = sumMatches[intCand.second + 1] - sumMatches[intCand.first];
			float newDivergence = 1.0f - float(rangeMatch) / rangeLen;

			OverlapRange newOvlp = ovlp;
			newOvlp.seqDivergence = newDivergence;
			size_t posQry = 0;
			size_t posTrg = 0;
			for (int i = 0; i < (int)decodedCigar.size(); ++i)
			{
				if (decodedCigar[i].op == '=' || decodedCigar[i].op == 'X')
				{
					posQry += decodedCigar[i].len;
					posTrg += decodedCigar[i].len;
				}
				else if (decodedCigar[i].op == 'I')
				{
					posQry += decodedCigar[i].len;
				}
				else
				{
					posTrg += decodedCigar[i].len;
				}
				if (i == intCand.first)
				{
					newOvlp.curBegin += posTrg;
					newOvlp.extBegin += posQry;
				}
				if (i == intCand.second)
				{
					newOvlp.curEnd = ovlp.curBegin + posTrg;
					newOvlp.extEnd = ovlp.extBegin + posQry;
				}
			}
			//TODO: updating score and k-mer matches?
			
			if (newOvlp.curRange() > minOverlap &&
				newOvlp.extRange() > minOverlap)
			{
				trimmedAlignments.push_back(newOvlp);
			}	
		}

		
		if (showAlignment)
		{
			Logger::get().debug() << "Adj from " << ovlp.curBegin 
				<< " " << ovlp.curRange() << 
				" " << ovlp.extBegin << " " << ovlp.extRange() 
				<< " " << errRate;

			for (auto& newOvlp : trimmedAlignments)
			{
				Logger::get().debug() << "      to " << newOvlp.curBegin 
					<< " " << newOvlp.curRange() << 
					" " << newOvlp.extBegin << " " << newOvlp.extRange() 
					<< " " << newOvlp.seqDivergence;
			}
		}

		return trimmedAlignments;
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

	//check if it's "almost trivial" match with intersecting sequence
	if (ovlp.curId == ovlp.extId)
	{
		int32_t intersect = std::min(ovlp.curEnd, ovlp.extEnd) - 
			   				std::max(ovlp.curBegin, ovlp.extBegin);
		if (intersect > ovlp.curRange() / 2) return false;
	}

	//check "strand skipping" PacBio pattern
	if (ovlp.curId == ovlp.extId.rc())
	{
		int32_t intersect = std::min(ovlp.curEnd, ovlp.extLen - ovlp.extBegin) - 
			   				std::max(ovlp.curBegin, ovlp.extLen - ovlp.extEnd);

		if (intersect > -_maxJump) outSuggestChimeric = true;
		if (intersect > ovlp.curRange() / 2) return false;
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

	template <class T>
	void shrinkAndClear(std::vector<T>& vec, float rate)
	{
		if (vec.empty()) return;
		//if ((float)vec.capacity() / vec.size() < rate) return;
		//Logger::get().debug() << "Shrink: " << vec.capacity() << " " << vec.size();

		size_t newCapacity = std::roundf((float)vec.capacity() / rate);
		vec = std::vector<T>();
		vec.reserve(newCapacity);
	}
}

//This implementation was inspired by Heng Li's minimap2 paper
//might be used in parallel
std::vector<OverlapRange> 
OverlapDetector::getSeqOverlaps(const FastaRecord& fastaRec, 
								bool& outSuggestChimeric,
								OvlpDivStats& divStats,
								int maxOverlaps) const
{
	//static std::ofstream fout("../kmers.txt");
	
	//const int MAX_LOOK_BACK = 50;
	const int kmerSize = Parameters::get().kmerSize;
	//const float minKmerSruvivalRate = std::exp(-_maxDivergence * kmerSize);
	const float minKmerSruvivalRate = 0.01;
	const float LG_GAP = 2;
	const float SM_GAP = 0.5;

	outSuggestChimeric = false;
	int32_t curLen = fastaRec.sequence.length();
	std::vector<int32_t> curFilteredPos;

	//cache memory-intensive containers as
	//many parallel memory allocations slow us down significantly
	//thread_local std::vector<KmerMatch> vecMatches;
	thread_local std::vector<KmerMatch> matchesList;
	thread_local std::vector<int32_t> scoreTable;
	thread_local std::vector<int32_t> backtrackTable;

	static ChunkPool<KmerMatch> sharedChunkPool;	//shared accoress threads
	BFContainer<KmerMatch> vecMatches(sharedChunkPool);

	//speed benchmarks
	thread_local float timeMemory = 0;
	thread_local float timeKmerIndexFirst = 0;
	thread_local float timeKmerIndexSecond = 0;
	thread_local float timeDp = 0;
	auto timeStart = std::chrono::system_clock::now();

	/*static std::mutex reportMutex;
	thread_local int threadId = rand() % 1000; 
	thread_local int numTicks = 0;
	++numTicks;
	thread_local auto prevReport = std::chrono::system_clock::now();
	if ((std::chrono::system_clock::now() - prevReport) > 
		std::chrono::seconds(60))
	{
		std::lock_guard<std::mutex> lock(reportMutex);
		prevReport = std::chrono::system_clock::now();
		Logger::get().debug() << ">Perf " << threadId
			<< " mem:" << timeMemory
			<< " kmFst:" << timeKmerIndexFirst 
			<< " kmSnd:" << timeKmerIndexSecond 
			<< " dp:" << timeDp << " ticks: " << numTicks; 
		Logger::get().debug() << ">Mem  " << threadId 
			<< " chunks:" << sharedChunkPool.numberChunks();

		timeMemory = 0;
		timeKmerIndexFirst = 0;
		timeKmerIndexSecond = 0;
		timeDp = 0;
		numTicks = 0;
	}*/
	timeStart = std::chrono::system_clock::now();

	//although once in a while shrink allocated memory size
	//thread_local auto prevCleanup = 
	//	std::chrono::system_clock::now() + std::chrono::seconds(rand() % 60);
	thread_local int prevCleanup = 0;
	if (++prevCleanup > 50)
	{
		prevCleanup = 0;
		shrinkAndClear(matchesList, 2);
		shrinkAndClear(scoreTable, 2);
		shrinkAndClear(backtrackTable, 2);
	}
	timeMemory += std::chrono::duration_cast<std::chrono::duration<float>>
						(std::chrono::system_clock::now() - timeStart).count();
	timeStart = std::chrono::system_clock::now();

	for (const auto& curKmerPos : IterKmers(fastaRec.sequence))
	{
		if (_vertexIndex.isRepetitive(curKmerPos.kmer))
		{
			curFilteredPos.push_back(curKmerPos.position);
		}
		if (!_vertexIndex.kmerFreq(curKmerPos.kmer)) continue;

		//FastaRecord::Id prevSeqId = FastaRecord::ID_NONE;
		for (const auto& extReadPos : _vertexIndex.iterKmerPos(curKmerPos.kmer))
		{
			//no trivial matches
			if ((extReadPos.readId == fastaRec.id &&
				extReadPos.position == curKmerPos.position)) continue;

			vecMatches.emplace_back(curKmerPos.position, 
									extReadPos.position,
									extReadPos.readId);
		}
	}
	timeKmerIndexFirst += std::chrono::duration_cast<std::chrono::duration<float>>
							(std::chrono::system_clock::now() - timeStart).count();
	timeStart = std::chrono::system_clock::now();

	std::sort(vecMatches.begin(), vecMatches.end(),
			  [](const KmerMatch& k1, const KmerMatch& k2)
			  {return k1.extId != k2.extId ? k1.extId < k2.extId : 
			  								 k1.curPos < k2.curPos;});

	timeKmerIndexSecond += std::chrono::duration_cast<std::chrono::duration<float>>
								(std::chrono::system_clock::now() - timeStart).count();
	timeStart = std::chrono::system_clock::now();

	const int STAT_WND = 10000;
	std::vector<OverlapRange> divStatWindows(curLen / STAT_WND + 1);

	std::vector<OverlapRange> detectedOverlaps;
	size_t extRangeBegin = 0;
	size_t extRangeEnd = 0;
	while(extRangeEnd < vecMatches.size())
	{
		if (maxOverlaps != 0 &&
			detectedOverlaps.size() >= (size_t)maxOverlaps) break;

		extRangeBegin = extRangeEnd;
		size_t uniqueMatches = 0;
		int32_t prevPos = 0;
		while (extRangeEnd < vecMatches.size() &&
			   vecMatches[extRangeBegin].extId == 
			   vecMatches[extRangeEnd].extId)
		{
			if (vecMatches[extRangeEnd].curPos != prevPos)
			{
				++uniqueMatches;
				prevPos = vecMatches[extRangeEnd].curPos;
			}
			++extRangeEnd;
		}
		if (uniqueMatches < minKmerSruvivalRate * _minOverlap) continue;

		matchesList.assign(vecMatches.begin() + extRangeBegin,
						   vecMatches.begin() + extRangeEnd);
		assert(matchesList.size() > 0 && 
			   matchesList.size() < (size_t)std::numeric_limits<int32_t>::max());

		FastaRecord::Id extId = matchesList.front().extId;
		int32_t extLen = _seqContainer.seqLen(extId);

		//pre-filtering
		int32_t minCur = matchesList.front().curPos;
		int32_t maxCur = matchesList.back().curPos;
		int32_t minExt = std::numeric_limits<int32_t>::max();
		int32_t maxExt = std::numeric_limits<int32_t>::min();
		for (const auto& match : matchesList)
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
					//int32_t gapCost = jumpDiv ? 
					//		kmerSize * jumpDiv + ilog2_32(jumpDiv) : 0;
					int32_t gapCost = (jumpDiv > 100 ? LG_GAP : SM_GAP) * jumpDiv;
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

		//backtracking
		std::vector<OverlapRange> extOverlaps;
		std::vector<int32_t> shifts;
		std::vector<std::pair<int32_t, int32_t>> kmerMatches;
		
		for (int32_t chainStart = backtrackTable.size() - 1; 
			 chainStart > 0; --chainStart)
		{
			if (backtrackTable[chainStart] == -1) continue;

			int32_t chainMaxScore = scoreTable[chainStart];
			int32_t lastMatch = chainStart;
			int32_t firstMatch = 0;

			int32_t chainLength = 0;
			shifts.clear();
			kmerMatches.clear();
			//int32_t totalMatch = kmerSize;
			//int32_t totalGap = 0;
			
			int32_t pos = chainStart;
			while (pos != -1)
			{
				//found a new maximum, shorten the chain end
				if (scoreTable[pos] > chainMaxScore)
				{
					chainMaxScore = scoreTable[pos];
					lastMatch = pos;
					chainLength = 0;
					shifts.clear();
					kmerMatches.clear();
				}

				firstMatch = pos;
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
					int32_t jumpDiv = abs((curNext - curPrev) - 
										  (extNext - extPrev));
					int32_t gapCost = (jumpDiv > 50) ? 2 * jumpDiv : jumpDiv;

					totalMatch += matchScore;
					totalGap += gapCost;
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

				assert(pos >= 0 && pos < (int32_t)backtrackTable.size());
				int32_t newPos = backtrackTable[pos];
				backtrackTable[pos] = -1;
				pos = newPos;
			}

			//Logger::get().debug() << chainStart - firstMatch << " " << lastMatch - firstMatch;

			OverlapRange ovlp(fastaRec.id, matchesList.front().extId,
							  matchesList[firstMatch].curPos, 
							  matchesList[firstMatch].extPos,
							  curLen, extLen);
			ovlp.curEnd = matchesList[lastMatch].curPos + kmerSize - 1;
			ovlp.extEnd = matchesList[lastMatch].extPos + kmerSize - 1;
			ovlp.score = scoreTable[lastMatch] - scoreTable[firstMatch] + 
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

				float normLen = std::max(ovlp.curRange(), 
										 ovlp.extRange()) - filteredPositions;
				float matchRate = (float)chainLength * 
								  _vertexIndex.getSampleRate() / normLen;
				matchRate = std::min(matchRate, 1.0f);
				//float repeatRate = (float)filteredPositions / ovlp.curRange();
				ovlp.seqDivergence = std::log(1 / matchRate) / kmerSize;
				ovlp.seqDivergence += _estimatorBias;

				if(_nuclAlignment)
				{
					auto trimmedOverlaps = 
						checkIdyAndTrim(ovlp, fastaRec.sequence, 
										_seqContainer.getSeq(extId),
										_maxDivergence, _minOverlap,
										/*show alignment*/ false);
					for (auto& trimOvlp : trimmedOverlaps)
					{
						extOverlaps.push_back(trimOvlp);
					}
				}
				else if (ovlp.seqDivergence < _maxDivergence)
				{
					extOverlaps.push_back(ovlp);
				}

				//statistics
				size_t wnd = ovlp.curBegin / STAT_WND;
				if (ovlp.curRange() > divStatWindows[wnd].curRange())
				{
					divStatWindows[wnd] = ovlp;
				}

				//benchmarking divergence
				/*float alnDiff = kswAlign(fastaRec.sequence, ovlp.curBegin, ovlp.curRange(),
										 _seqContainer.getSeq(extId), ovlp.extBegin, 
										  ovlp.extRange(), 1, -2, 2, 1, false);
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
			const OverlapRange* maxOvlp = nullptr;
			for (const auto& ovlp : extOverlaps)
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
			
			for (const auto& ovlp : extOverlaps)
			{
				bool isContained = false;
				for (const auto& prim : primOverlaps)
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
			for (const auto& ovlp : primOverlaps)
			{
				detectedOverlaps.push_back(ovlp);
			}
		}
	}

	timeDp += std::chrono::duration_cast<std::chrono::duration<float>>
				(std::chrono::system_clock::now() - timeStart).count();
	timeStart = std::chrono::system_clock::now();

	for (const auto& ovlp : divStatWindows)
	{
		if (ovlp.curRange() > 0)
		{
			divStats.add(ovlp.seqDivergence);
		}
	}
	return detectedOverlaps;
}

bool OverlapContainer::hasSelfOverlaps(FastaRecord::Id readId)
{
	this->lazySeqOverlaps(readId);
	if (!readId.strand()) readId = readId.rc();
	return _overlapIndex.find(readId).suggestChimeric;
}


std::vector<OverlapRange> 
	OverlapContainer::quickSeqOverlaps(FastaRecord::Id readId, int maxOverlaps)
{
	bool suggestChimeric;
	const FastaRecord& record = _queryContainer.getRecord(readId);
	return _ovlpDetect.getSeqOverlaps(record, suggestChimeric, 
									  _divergenceStats, maxOverlaps);
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
											   _divergenceStats,
											   _ovlpDetect._maxCurOverlaps);
	overlaps.shrink_to_fit();

	std::vector<OverlapRange> revOverlaps;
	revOverlaps.reserve(overlaps.size());
	for (const auto& ovlp : overlaps) revOverlaps.push_back(ovlp.complement());

	_overlapIndex.update_fn(readId,
		[&wrapper, &overlaps, &revOverlaps, &suggestChimeric, this]
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
	for (const auto& seqIt : _overlapIndex.lock_table()) 
	{
		allSeqs.push_back(seqIt.first);
		allSeqs.push_back(seqIt.first.rc());
	}

	int totalOverlaps = 0;
	for (const auto& seq : allSeqs)
	{
		const auto& curOvlps = this->unsafeSeqOverlaps(seq);
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
		for (const auto& ovlp : ovlpsToAdd)
		{
			this->unsafeSeqOverlaps(ovlp.curId).push_back(ovlp);
		}
	}
}


void OverlapContainer::findAllOverlaps()
{
	//Logger::get().info() << "Finding overlaps:";
	std::vector<FastaRecord::Id> allQueries;
	for (const auto& seq : _queryContainer.iterSeqs())
	{
		allQueries.push_back(seq.id);
	}

	std::mutex indexMutex;
	std::function<void(const FastaRecord::Id&)> indexUpdate = 
	[this] (const FastaRecord::Id& seqId)
	{
		this->lazySeqOverlaps(seqId);	//automatically stores overlaps
	};
	processInParallel(allQueries, indexUpdate, 
					  Parameters::get().numThreads, true);
	this->ensureTransitivity(false);

	int numOverlaps = 0;
	for (const auto& seqOvlps : _overlapIndex.lock_table()) 
	{
		numOverlaps += seqOvlps.second.fwdOverlaps->size() * 2;
	}
	Logger::get().debug() << "Found " << numOverlaps << " overlaps";

	this->filterOverlaps();

	numOverlaps = 0;
	for (const auto& seqOvlps : _overlapIndex.lock_table()) 
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
	for (const auto& seq : _queryContainer.iterSeqs())
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
		for (const auto& cluster : clusters)
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


void OverlapContainer::estimateOverlaperParameters()
{
	Logger::get().debug() << "Estimating k-mer identity bias";

	const int NEDEED_OVERLAPS = 1000;
	const int MAX_SEQS = 1000;

	std::vector<FastaRecord::Id> readsToCheck;
	for (size_t i = 0; i < MAX_SEQS; ++i) 
	{
		size_t randId = rand() % _queryContainer.iterSeqs().size();
		readsToCheck.push_back(_queryContainer.iterSeqs()[randId].id);
	}

	std::mutex storageMutex;
	std::vector<float> biases;
	std::vector<float> trueDivergence;
	std::function<void(const FastaRecord::Id& seqId)> computeParallel =
	[this, &storageMutex, &biases, &trueDivergence] (const FastaRecord::Id& seqId)
	{
		auto overlaps = this->quickSeqOverlaps(seqId, /*max ovlps*/ 0);
		for (const auto& ovlp : overlaps)
		{
			float trueDiv = 
				getAlignmentIdy(ovlp, _queryContainer.getSeq(seqId),
								_ovlpDetect._seqContainer.getSeq(ovlp.extId),
								false);

			std::lock_guard<std::mutex> lock(storageMutex);
			biases.push_back(trueDiv - ovlp.seqDivergence);
			trueDivergence.push_back(trueDiv);
			if (biases.size() >= NEDEED_OVERLAPS) return;
		}
	};
	processInParallel(readsToCheck, computeParallel, 
					  Parameters::get().numThreads, false);

	if (!biases.empty())
	{
		_kmerIdyEstimateBias = median(biases);
		_meanTrueOvlpDiv = median(trueDivergence);

		//set the parameters and reset statistics
		_ovlpDetect._estimatorBias = _kmerIdyEstimateBias;
		_divergenceStats.vecSize = 0;
	}
	else
	{
		Logger::get().warning() << "No overlaps found - unable to estimate parameters";
		_meanTrueOvlpDiv = 0.5f;
		_kmerIdyEstimateBias = 0.0f;
	}

	Logger::get().debug() << "Median overlap divergence: " << _meanTrueOvlpDiv;
	Logger::get().debug() << "K-mer estimate bias (true - est): " << _kmerIdyEstimateBias;
}


void OverlapContainer::setRelativeDivergenceThreshold(float relThreshold)
{
	_ovlpDetect._maxDivergence = _meanTrueOvlpDiv + relThreshold;
	Logger::get().debug() << "Max divergence threshold set to " 
		<< _ovlpDetect._maxDivergence;
}

void OverlapContainer::overlapDivergenceStats()
{
	this->overlapDivergenceStats(_divergenceStats, 
								 _ovlpDetect._maxDivergence);
}

void OverlapContainer::overlapDivergenceStats(const OvlpDivStats& stats,
											  float divCutoff)
{
	std::vector<float> ovlpDivergence(stats.divVec.begin(),
									  stats.divVec.begin() + stats.vecSize);
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
	int threshold = divCutoff * mult * 100;
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
			else if (i != threshold)
			{
				histString += ' ';
			}
			else
			{
				histString += '|';
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
	for (const auto& seqIt : _overlapIndex.lock_table()) 
	{
		allSeqs.push_back(seqIt.first);
		allSeqs.push_back(seqIt.first.rc());
	}

	for (const auto& seq : allSeqs)
	{
		std::vector<Interval<const OverlapRange*>> intervals;
		auto& overlaps = this->unsafeSeqOverlaps(seq);
		for (const auto& ovlp : overlaps)
		{
			intervals.emplace_back(ovlp.curBegin, ovlp.curEnd, &ovlp);
		}
		_ovlpTree[seq] = IntervalTree<const OverlapRange*>(intervals);
	}
}

std::vector<Interval<const OverlapRange*>> 
	OverlapContainer::getCoveringOverlaps(FastaRecord::Id seqId, 
								  int32_t start, int32_t end) const
{
	return _ovlpTree.at(seqId).findOverlapping(start, end);
}
