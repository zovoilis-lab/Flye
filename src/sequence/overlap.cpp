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

	if (_checkOverhang && ovlp.curId != ovlp.extId.rc())
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

namespace
{
	struct DPRecord
	{
		OverlapRange ovlp;
		std::vector<int32_t> shifts;
	};
	
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

std::vector<OverlapRange> 
OverlapDetector::getSeqOverlaps(const FastaRecord& fastaRec, 
								bool uniqueExtensions) const
{
	std::unordered_map<FastaRecord::Id, 
					   std::vector<DPRecord>> activePaths;
	std::set<size_t> eraseMarks;
	size_t curLen = fastaRec.sequence.length();
	//std::vector<KmerPosition> solidKmersCache;

	//for all kmers in this read
	for (auto curKmerPos : IterKmers(fastaRec.sequence))
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
			if (extReadPos.readId == fastaRec.id &&
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
				JumpRes jumpResult = 
					this->jumpTest(extPaths[pathId].ovlp.curEnd, curPos,
								   extPaths[pathId].ovlp.extEnd, extPos);
				int32_t jumpLength = curPos - extPaths[pathId].ovlp.curBegin;
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
				extPaths[maxCloseId].ovlp.curEnd = curPos;
				extPaths[maxCloseId].ovlp.extEnd = extPos;
				extPaths[maxCloseId].shifts.push_back(curPos - extPos);
				//++extPaths[maxCloseId].score;
			}
			//update the best far extension, keep the old path as a copy
			if (extendsFar)
			{
				extPaths.push_back(extPaths[maxFarId]);
				extPaths.back().ovlp.curEnd = curPos;
				extPaths.back().ovlp.extEnd = extPos;
				extPaths.back().shifts.push_back(curPos - extPos);
				//++extPaths.back().score;
			}
			//if no extensions possible (or there are no active paths), start a new path
			if (!extendsClose && !extendsFar &&
				(this->goodStart(curPos, extPos, curLen, extLen) ||
				fastaRec.id == extReadPos.readId.rc()))	//TODO: temporary bypass overhang
			{
				OverlapRange ovlp(fastaRec.id, extReadPos.readId,
								  curPos, extPos);
				extPaths.push_back({ovlp});
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
		OverlapRange maxOverlap;
		bool passedTest = false;
		for (auto& dpRec : ap.second)
		{
			if (this->overlapTest(dpRec.ovlp, curLen, extLen))
			{
				dpRec.ovlp.leftShift = median(dpRec.shifts);
				dpRec.ovlp.rightShift = extLen - curLen + 
										dpRec.ovlp.leftShift;

				if (!uniqueExtensions)
				{
					detectedOverlaps.push_back(dpRec.ovlp);
				}
				else
				{
					passedTest = true;
					if (dpRec.ovlp.curRange() > maxOverlap.curRange())
					{
						maxOverlap = dpRec.ovlp;
					}
				}
			}
		}
		if (uniqueExtensions && passedTest)
		{
			detectedOverlaps.push_back(maxOverlap);
		}
	}

	return detectedOverlaps;

}

const std::vector<OverlapRange>&
OverlapContainer::getSeqOverlaps(FastaRecord::Id readId)
{
	if (!_cached.count(readId))
	{
		_cached.insert(readId);
		_cached.insert(readId.rc());

		const FastaRecord& record = _queryContainer.getIndex().at(readId);
		auto overlaps = _ovlpDetect.getSeqOverlaps(record, _onlyMax);

		auto& fwdOverlaps = _overlapIndex[record.id];
		auto& revOverlaps = _overlapIndex[record.id.rc()];

		std::unordered_set<FastaRecord::Id> extisting;
		if (_onlyMax)
		{
			for (auto& ovlp : fwdOverlaps) extisting.insert(ovlp.extId);
		}

		for (auto& ovlp : overlaps)
		{
			if (_onlyMax && extisting.count(ovlp.extId)) continue;

			int32_t curLen = _queryContainer.seqLen(ovlp.curId);
			int32_t extLen = _queryContainer.seqLen(ovlp.extId);

			auto revOvlp = ovlp.reverse();
			fwdOverlaps.push_back(ovlp);
			revOverlaps.push_back(ovlp.complement(curLen, extLen));
			_overlapIndex[revOvlp.curId].push_back(revOvlp);
			_overlapIndex[revOvlp.curId.rc()]
					.push_back(revOvlp.complement(extLen, curLen));
		}
	}
	return _overlapIndex[readId];
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
		auto& fastaRec = _queryContainer.getIndex().at(seqId);
		auto overlaps = _ovlpDetect.getSeqOverlaps(fastaRec, false);

		indexMutex.lock();
		for (auto& ovlp : overlaps)
		{
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
	processInParallel(allQueries, indexUpdate, 
					  Parameters::get().numThreads, true);
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
