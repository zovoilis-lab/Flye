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
#include "../common/parallel.h"


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
bool OverlapDetector::overlapTest(const OverlapRange& ovlp) const
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

	//for all kmers in this read
	for (auto curKmerPos : IterKmers(fastaRec.sequence))
	{
		if (!_vertexIndex.isSolid(curKmerPos.kmer)) continue;

		int32_t curPos = curKmerPos.position;

		//for all other occurences of this kmer (extension candidates)
		for (const auto& extReadPos : _vertexIndex.byKmer(curKmerPos.kmer))
		{
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
			}
			//update the best far extension, keep the old path as a copy
			if (extendsFar)
			{
				extPaths.push_back(extPaths[maxFarId]);
				extPaths.back().ovlp.curEnd = curPos;
				extPaths.back().ovlp.extEnd = extPos;
				extPaths.back().shifts.push_back(curPos - extPos);
			}
			//if no extensions possible (or there are no active paths), start a new path
			if (!extendsClose && !extendsFar &&
				(this->goodStart(curPos, extPos, curLen, extLen) ||
				fastaRec.id == extReadPos.readId.rc()))	//TODO: temporary bypass overhang
			{
				OverlapRange ovlp(fastaRec.id, extReadPos.readId,
								  curPos, extPos, curLen, extLen);
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
		DPRecord* maxOverlap = nullptr;
		bool passedTest = false;
		for (auto& dpRec : ap.second)
		{
			if (this->overlapTest(dpRec.ovlp))
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
					if (!maxOverlap || 
						dpRec.ovlp.curRange() > maxOverlap->ovlp.curRange())
					{
						maxOverlap = &dpRec;
					}
				}
			}
		}
		if (uniqueExtensions && passedTest)
		{
			maxOverlap->ovlp.leftShift = median(maxOverlap->shifts);
			maxOverlap->ovlp.rightShift = extLen - curLen + 
									maxOverlap->ovlp.leftShift;
			detectedOverlaps.push_back(maxOverlap->ovlp);
		}
	}

	return detectedOverlaps;

}

const std::vector<OverlapRange>&
OverlapContainer::lazySeqOverlaps(FastaRecord::Id readId)
{
	if (!_cached.count(readId))
	{
		const FastaRecord& record = _queryContainer.getIndex().at(readId);
		auto overlaps = _ovlpDetect.getSeqOverlaps(record, _onlyMax);
		this->storeOverlaps(overlaps, readId);	
	}
	return _overlapIndex[readId];
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
		this->storeOverlaps(overlaps, seqId);
		indexMutex.unlock();
	};

	processInParallel(allQueries, indexUpdate, 
					  Parameters::get().numThreads, true);

	this->filterOverlaps();
}


void OverlapContainer::filterOverlaps()
{
	std::vector<FastaRecord::Id> seqIds;
	for (auto seqIndex : _queryContainer.getIndex())
	{
		seqIds.push_back(seqIndex.first);
	}

	Logger::get().debug() << "Filtering overlaps";
	std::atomic<size_t> numIdent(0);
	std::atomic<size_t> numContained(0);
	std::function<void(const FastaRecord::Id& seqId)> filterParallel =
	[&numIdent, &numContained, this] (const FastaRecord::Id& seqId)
	{
		auto& overlaps = _overlapIndex[seqId];
		std::sort(overlaps.begin(), overlaps.end(), 
				  [](const OverlapRange& o1, const OverlapRange& o2)
				  {return o1.curBegin < o2.curBegin;});

		std::unordered_map<int32_t, std::vector<OverlapRange*>> ovlpsByStart;
		for (auto& ovlp : overlaps)
		{
			bool found = false;
			for (auto& otherOvlp : ovlpsByStart[ovlp.curBegin])
			{
				if (ovlp.equals(*otherOvlp))
				{
					found = true;
					break;
				}
			}
			if (!found) ovlpsByStart[ovlp.curBegin].push_back(&ovlp);
		}
		numIdent += overlaps.size() - ovlpsByStart.size();
		overlaps.clear();
		for (auto& pointOvlps : ovlpsByStart)
		{
			for (auto& ovlp : pointOvlps.second) overlaps.push_back(*ovlp);
		}
		/////////////
		
		std::unordered_set<OverlapRange*> contained;
		for (auto& ovlp : overlaps)
		{
			for (auto& otherOvlp : overlaps)
			{
				if (ovlp.containedBy(otherOvlp))
				{
					contained.insert(&ovlp);
					break;
				}
			}
		}

		std::vector<OverlapRange> nonContained;
		for (auto& ovlp : overlaps)
		{
			if (!contained.count(&ovlp)) nonContained.push_back(ovlp);
		}
		numContained += overlaps.size() - nonContained.size();
		overlaps = std::move(nonContained);

	};
	processInParallel(seqIds, filterParallel, 
					  Parameters::get().numThreads, false);

	Logger::get().debug() << "Filtered " << numIdent << " identical and "
		<< numContained << " contained overlaps";
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
