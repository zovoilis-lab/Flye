//(c) 2016 by Authors
//This file is a part of ABruijn program.
//Released under the BSD license (see LICENSE file)

#include <set>
#include <iostream>
#include <cassert>
#include <algorithm>
#include <thread>

#include "overlap.h"
#include "config.h"
#include "parallel.h"


//pre-filtering
bool OverlapDetector::goodStart(int32_t curPos, int32_t extPos, 
								int32_t curLen, int32_t extLen) const
{	
	if (_checkOverhang && 
		std::min(curPos, extPos) > _maxOverhang) return false;

	if (extPos > extLen - _minOverlap ||
	    curPos > curLen - _minOverlap) return false;

	return true;
}


OverlapDetector::JumpRes 
OverlapDetector::jumpTest(int32_t curPrev, int32_t curNext,
						  int32_t extPrev, int32_t extNext) const
{
	if (curNext - curPrev > Constants::maximumJump) return J_END;

	if (0 < curNext - curPrev && curNext - curPrev < _maxJump &&
		0 < extNext - extPrev && extNext - extPrev < _maxJump)
	{
		if (abs((curNext - curPrev) - (extNext - extPrev)) 
			< _maxJump / Constants::closeJumpRate)
		{
			return J_CLOSE;
		}
		if (abs((curNext - curPrev) - (extNext - extPrev)) 
			< _maxJump / Constants::farJumpRate)
		{
			return J_FAR;
		}
	}
	return J_INCONS;
}


//Check if it is a proper overlap
bool OverlapDetector::overlapTest(const OverlapRange& ovlp, int32_t curLen, 
								  int32_t extLen) const
{

	if (ovlp.curRange() < _minOverlap || 
		ovlp.extRange() < _minOverlap) 
	{
		return false;
	}

	float lengthDiff = abs(ovlp.curRange() - ovlp.extRange());
	float meanLength = (ovlp.curRange() + ovlp.extRange()) / 2.0f;
	if (lengthDiff > meanLength / Constants::overlapDivergenceRate)
	{
		return false;
	}

	if (_checkOverhang)
	{
		if (std::min(ovlp.curBegin, ovlp.extBegin) > 
			_maxOverhang) 
		{
			return false;
		}
		if (std::min(curLen - ovlp.curEnd, extLen - ovlp.extEnd) > 
			_maxOverhang)
		{
			return false;
		}
	}

	return true;
}

std::vector<OverlapRange> 
OverlapDetector::getSeqOverlaps(const std::string& sequence,
								FastaRecord::Id idTrivial) const
{
	std::unordered_map<FastaRecord::Id, 
					   std::vector<OverlapRange>> activePaths;
	std::set<size_t> eraseMarks;
	size_t curLen = sequence.length();
	//std::vector<KmerPosition> solidKmersCache;

	//for all kmers in this read
	for (auto curKmerPos : IterKmers(sequence))
	{
		if (!_vertexIndex.isSolid(curKmerPos.kmer)) continue;

		int32_t curPos = curKmerPos.position;
		//solidKmersCache.push_back(curKmerPos);

		//for all other occurences of this kmer (extension candidates)
		for (const auto& extReadPos : _vertexIndex.byKmer(curKmerPos.kmer))
		{
			if (_vertexIndex.isRepetitive(curKmerPos.kmer) &&
				!activePaths.count(extReadPos.readId))
			{
				continue;
			}

			//no trivial matches
			if (extReadPos.readId == idTrivial &&
				extReadPos.position == curPos) continue;

			size_t extLen = _seqContainer.seqLen(extReadPos.readId);
			if (extLen < (size_t)_minOverlap) continue;

			int32_t extPos = extReadPos.position;
			auto& extPaths = activePaths[extReadPos.readId];

			size_t maxCloseId = 0;
			size_t maxFarId = 0;
			int32_t maxCloseLen = 0;
			int32_t maxFarLen = 0;
			bool extendsClose = false;
			bool extendsFar = false;
			eraseMarks.clear();

			//searching for longest possible extension
			for (size_t pathId = 0; pathId < extPaths.size(); ++pathId)
			{
				JumpRes jumpResult = this->jumpTest(extPaths[pathId].curEnd, curPos,
													extPaths[pathId].extEnd, extPos);
				int32_t jumpLength = curPos - extPaths[pathId].curBegin;
				//int32_t jumpLength = extPaths[pathId].score + 1;

				switch (jumpResult)
				{
					case J_END:
						break;
					case J_INCONS:
						break;
					case J_CLOSE:
						eraseMarks.insert(pathId);
						if (jumpLength > maxCloseLen)
						{
							extendsClose = true;
							maxCloseId = pathId;	
							maxCloseLen = jumpLength;
						}
						break;
					case J_FAR:
						if (jumpLength > maxFarLen)
						{
							extendsFar = true;
							maxFarId = pathId;
							maxFarLen = jumpLength;
						}
						break;
				}
			}
			//update the best close extension
			if (extendsClose)
			{
				eraseMarks.erase(maxCloseId);
				extPaths[maxCloseId].curEnd = curPos;
				extPaths[maxCloseId].extEnd = extPos;
				//++extPaths[maxCloseId].score;
			}
			//update the best far extension, keep the old path as a copy
			if (extendsFar)
			{
				extPaths.push_back(extPaths[maxFarId]);
				extPaths.back().curEnd = curPos;
				extPaths.back().extEnd = extPos;
				//++extPaths.back().score;
			}
			//if no extensions possible (or there are no active paths), start a new path
			if (!extendsClose && !extendsFar &&
				this->goodStart(curPos, extPos, curLen, extLen))
			{
				extPaths.emplace_back(FastaRecord::ID_NONE, extReadPos.readId,
									  curPos, extPos);
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
	
	std::vector<OverlapRange> detectedOverlaps;
	for (auto& ap : activePaths)
	{
		size_t extLen = _seqContainer.seqLen(ap.first);
		for (auto& ovlp : ap.second)
		{
			if (this->overlapTest(ovlp, curLen, extLen))
			{
				detectedOverlaps.push_back(ovlp);
			}
		}
	}

	return detectedOverlaps;

}

void OverlapContainer::findAllOverlaps()
{
	Logger::get().info() << "Finding overlaps:";
	std::vector<FastaRecord::Id> allQueries;
	for (auto& hashPair : _queryContainer.getIndex())
	{
		allQueries.push_back(hashPair.first);
	}

	std::mutex indexMutex;
	std::function<void(const FastaRecord::Id&)> indexUpdate = 
	[this, &indexMutex] (const FastaRecord::Id& seqId)
	{
		const std::string& seq = _queryContainer.getSeq(seqId);
		auto overlaps = _ovlpDetect.getSeqOverlaps(seq, seqId);

		indexMutex.lock();
		for (auto& ovlp : overlaps)
		{
			ovlp.curId = seqId;
			_overlapIndex[seqId].push_back(ovlp);
		}
		indexMutex.unlock();
	};

		/*
		if (_overlapMatrix.contains(std::make_tuple(ovlp.curId, 
													ovlp.extId)))
		{
			continue;
		}
		//detected overlap
		_overlapMatrix[std::make_tuple(ovlp.curId, ovlp.extId)] = true;
		_overlapIndex[ovlp.curId].push_back(ovlp);
		//in opposite direction
		ovlp = ovlp.reverse();
		_overlapMatrix[std::make_tuple(ovlp.curId, ovlp.extId)] = true;
		_overlapIndex[ovlp.curId].push_back(ovlp);
		//on reverse strands
		ovlp = ovlp.complement();
		_overlapMatrix[std::make_tuple(ovlp.curId, ovlp.extId)] = true;
		_overlapIndex[ovlp.curId].push_back(ovlp);
		//opposite again
		ovlp = ovlp.reverse();
		_overlapMatrix[std::make_tuple(ovlp.curId, ovlp.extId)] = true;
		_overlapIndex[ovlp.curId].push_back(ovlp);
		*/
	processInParallel(allQueries, indexUpdate, Parameters::numThreads, true);
}


void OverlapContainer::saveOverlaps(const std::string& filename)
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
}

/*
namespace
{
	template<typename T>
	T median(std::vector<T>& vec)
	{
		std::sort(vec.begin(), vec.end());
		//NOTE: there's a bug in libstdc++ nth_element, 
		//that sometimes leads to a segfault
		//std::nth_element(vec.begin(), vec.begin() + vec.size() / 2, 
		//				 vec.end());
		return vec[vec.size() / 2];
	}
}

void OverlapDetector::addOverlapShifts(OverlapRange& ovlp,
									   const std::vector<KmerPosition>& 
									   		solidKmersCache,
									   int32_t curLen, int32_t extLen) const
{
	std::vector<int32_t> ovlpShifts;
	for (auto curKmer : solidKmersCache)
	{
		if (ovlp.curBegin <= curKmer.position &&
			curKmer.position <= ovlp.curEnd)
		{
			for (const auto& extKmer : _vertexIndex.byKmer(curKmer.kmer))
			{
				if (extKmer.readId == ovlp.extId &&
				    ovlp.extBegin <= extKmer.position &&
					extKmer.position <= ovlp.extEnd)
				{
					ovlpShifts.push_back(curKmer.position - extKmer.position);
				}
			}
		}
	}

	ovlp.leftShift = median(ovlpShifts);
	ovlp.rightShift = extLen - curLen + ovlp.leftShift;
}*/
