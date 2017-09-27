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
#include "../common/disjoint_set.h"


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
		DPRecord() {}
		DPRecord(const OverlapRange& ovlp): ovlp(ovlp) {}

		OverlapRange ovlp;
		//int score;
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

	struct pairhash 
	{
	public:
		template <typename T, typename U>
		std::size_t operator()(const std::pair<T, U> &x) const
		{
			return std::hash<T>()(x.first) ^ std::hash<U>()(x.second);
		}
	};
}

std::vector<OverlapRange> 
OverlapDetector::getSeqOverlaps(const FastaRecord& fastaRec, 
								bool uniqueExtensions) const
{
	std::unordered_map<FastaRecord::Id, 
					   std::vector<DPRecord>> activePaths;
	std::unordered_map<FastaRecord::Id, 
					   std::vector<DPRecord>> completedPaths;
	std::set<size_t> eraseMarks;
	size_t curLen = fastaRec.sequence.length();

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

			size_t extLen = _seqContainer.seqLen(extReadPos.readId);
			if (extLen < (size_t)_minOverlap) continue;

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

				int32_t jumpLength = curPos - extPaths[pathId].ovlp.curEnd;
				int32_t gapScore = (jumpLength - Constants::gapJump) / 
													   Constants::penaltyWindow;
				if (jumpLength < Constants::gapJump) gapScore = 1;

				int32_t jumpScore = extPaths[pathId].ovlp.score + gapScore;

				switch (jumpResult)
				{
					case J_END:
						eraseMarks.insert(pathId);
						if (this->overlapTest(extPaths[pathId].ovlp))
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
				extPaths[maxCloseId].ovlp.curEnd = curPos;
				extPaths[maxCloseId].ovlp.extEnd = extPos;
				extPaths[maxCloseId].ovlp.score = maxCloseScore;
				extPaths[maxCloseId].shifts.push_back(curPos - extPos);
				if (_keepAlignment)
				{
					extPaths[maxCloseId].ovlp.kmerMatches
										.emplace_back(curPos, extPos);
				}
			}
			//update the best far extension, keep the old path as a copy
			if (extendsFar)
			{
				extPaths.push_back(extPaths[maxFarId]);
				extPaths.back().ovlp.curEnd = curPos;
				extPaths.back().ovlp.extEnd = extPos;
				extPaths.back().ovlp.score = maxFarScore;
				extPaths.back().shifts.push_back(curPos - extPos);
				if (_keepAlignment)
				{
					extPaths.back().ovlp.kmerMatches
										.emplace_back(curPos, extPos);
				}
			}
			//if no extensions possible (or there are no active paths), start a new path
			if (!extendsClose && !extendsFar &&
				(this->goodStart(curPos, extPos, curLen, extLen) ||
				fastaRec.id == extReadPos.readId.rc()))	//TODO: temporary bypass overhang
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
			if (this->overlapTest(dpRec.ovlp))
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
		size_t extLen = _seqContainer.seqLen(ap.first);
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

std::vector<OverlapRange>
OverlapContainer::seqOverlaps(FastaRecord::Id seqId) const
{
	const FastaRecord& record = _queryContainer.getIndex().at(seqId);
	return _ovlpDetect.getSeqOverlaps(record, _onlyMax);
}

const std::vector<OverlapRange>&
OverlapContainer::lazySeqOverlaps(FastaRecord::Id readId)
{
	_indexMutex.lock();
	if (!_cached.count(readId))
	{
		_indexMutex.unlock();
		auto overlaps = this->seqOverlaps(readId);
		_indexMutex.lock();
		this->storeOverlaps(overlaps, readId);
	}
	_indexMutex.unlock();
	return _overlapIndex.at(readId);
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
		auto overlaps = _ovlpDetect.getSeqOverlaps(fastaRec, false);

		indexMutex.lock();
		this->storeOverlaps(overlaps, seqId);
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

				if (curDiff < Constants::maximumJump / 10 && 
					extDiff < Constants::maximumJump / 10) 
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
