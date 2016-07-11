//(c) 2016 by Authors
//This file is a part of ABruijn program.
//Released under the BSD license (see LICENSE file)

#include <set>
#include <iostream>
#include <cassert>
#include <algorithm>
#include <thread>

#include "overlap.h"

void OverlapDetector::findAllOverlaps(size_t numThreads)
{
	Logger::get().info() << "Finding overlaps:";
	_overlapMatrix.clear();
	_overlapMatrix.reserve(_seqContainer.getIndex().size() * _coverage);
	_jobQueue.clear();
	_nextJob = 0;

	for (const auto& seqHash : _seqContainer.getIndex())
	{
		_jobQueue.push_back(seqHash.first);
	}

	std::vector<std::thread> threads(numThreads);
	for (size_t i = 0; i < threads.size(); ++i)
	{
		threads[i] = std::thread(&OverlapDetector::parallelWorker, this);
	}
	for (size_t i = 0; i < threads.size(); ++i)
	{
		threads[i].join();
	}

	_overlapMatrix.clear();
}


void OverlapDetector::parallelWorker()
{
	_fetchMutex.lock();
	while (true)
	{
		if (_nextJob == _jobQueue.size())
		{
			_fetchMutex.unlock();
			return;
		}
		_progress.advance();
		FastaRecord::Id readId = _jobQueue[_nextJob++];
		_overlapIndex[readId];	//empty vector by default

		if (_seqContainer.seqLen(readId) < (size_t)_minimumOverlap ||
			!readId.strand()) 
		{
			continue;
		}

		_fetchMutex.unlock();
		auto detectedOverlaps = this->getReadOverlaps(readId);
		_fetchMutex.lock();
		
		for (auto ovlp : detectedOverlaps)
		{
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
		}
	}
}


void OverlapDetector::saveOverlaps(const std::string filename)
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

void OverlapDetector::loadOverlaps(const std::string filename)
{
	std::ifstream fin(filename);
	if (!fin.is_open())
	{
		throw std::runtime_error("Can't open overlaps file");
	}

	for (auto& seqHash : _seqContainer.getIndex())
	{
		_overlapIndex[seqHash.first];	//empty vector by default
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


//pre-filtering
bool OverlapDetector::goodStart(int32_t curPos, int32_t extPos, 
								int32_t curLen, int32_t extLen) const
{	
	return std::min(curPos, extPos) < _maximumOverhang && 
		   extPos < extLen - _minimumOverlap &&
		   curPos < curLen - _minimumOverlap;
}

OverlapDetector::JumpRes 
OverlapDetector::jumpTest(int32_t curPrev, int32_t curNext,
						  int32_t extPrev, int32_t extNext) const
{
	static const int CLOSE_FRAC = 100;
	static const int FAR_FRAC = 2;
	
	if (curNext - curPrev > _maximumJump) return J_END;

	if (0 < curNext - curPrev && curNext - curPrev < _maximumJump &&
		0 < extNext - extPrev && extNext - extPrev < _maximumJump)
	{
		if (abs((curNext - curPrev) - (extNext - extPrev)) 
			< _maximumJump / CLOSE_FRAC)
		{
			return J_CLOSE;
		}
		if (abs((curNext - curPrev) - (extNext - extPrev)) 
			< _maximumJump / FAR_FRAC)
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
	if (ovlp.curRange() < _minimumOverlap || ovlp.extRange() < _minimumOverlap)
		return false;
	if (abs(ovlp.curRange() - ovlp.extRange()) > _maximumJump)
		return false;
	if (ovlp.curId == ovlp.extId.rc()) return true;	//FIXME: adhoc solution for chimeras
	if (std::min(ovlp.curBegin, ovlp.extBegin) > _maximumOverhang)
		return false;
	if (std::min(curLen - ovlp.curEnd, extLen - ovlp.extEnd) > _maximumOverhang)
		return false;
	return true;
}

//Getting all possible overlaps
//based on the shared kmers (common jump-paths)
std::vector<OverlapRange> 
OverlapDetector::getReadOverlaps(FastaRecord::Id currentReadId) const
{
	std::unordered_map<FastaRecord::Id, 
					   std::vector<OverlapRange>> activePaths;
	std::set<size_t> eraseMarks;
	size_t curLen = _seqContainer.seqLen(currentReadId);
	std::vector<KmerPosition> solidKmersCache;

	//for all kmers in this read
	int curKmerId = 0;
	for (auto curKmerPos : IterSolidKmers(currentReadId))
	{
		++curKmerId;
		int32_t curPos = curKmerPos.position;
		solidKmersCache.push_back(curKmerPos);

		//for all other occurences of this kmer (extension candidates)
		for (const auto& extReadPos : _vertexIndex.byKmer(curKmerPos.kmer))
		{
			int32_t extPos = extReadPos.position;
			size_t extLen = _seqContainer.seqLen(extReadPos.readId);
			//don't want self-overlaps
			if (extReadPos.readId == currentReadId) continue;
			if (extLen < (size_t)_minimumOverlap) continue;

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
			}
			//update the best far extension, keep the old path as a copy
			if (extendsFar)
			{
				extPaths.push_back(extPaths[maxFarId]);
				extPaths.back().curEnd = curPos;
				extPaths.back().extEnd = extPos;
			}
			//if no extensions possible (or there are no active paths), start a new path
			if (!extendsClose && !extendsFar &&
				(this->goodStart(curPos, extPos, curLen, extLen) ||
				currentReadId == extReadPos.readId.rc()))
			{
				extPaths.emplace_back(currentReadId, extReadPos.readId,
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
	//std::vector<OverlapRange> debugOverlaps;
	for (auto& ap : activePaths)
	{
		size_t extLen = _seqContainer.seqLen(ap.first);
		OverlapRange maxOverlap;
		//OverlapRange debugOverlap;
		bool passedTest = false;
		for (auto& ovlp : ap.second)
		{
			if (this->overlapTest(ovlp, curLen, extLen))
			{
				passedTest = true;
				if (maxOverlap.curRange() < ovlp.curRange()) maxOverlap = ovlp;
			}
			//if (debugOverlap.curRange() < ovlp.curRange()) debugOverlap = ovlp;
		}

		if (passedTest)
		{
			this->addOverlapShifts(maxOverlap, solidKmersCache, curLen, extLen);
			detectedOverlaps.push_back(maxOverlap);
		}
		//if (debugOverlap.curRange() > 3000) debugOverlaps.push_back(debugOverlap);
	}
	
	
	/*
	if (!debugOverlaps.empty())
	{
		_logMutex.lock();
		Logger::get().debug() << "Ovlps for " 
					<< _seqContainer.seqName(currentReadId);
		for (auto& ovlp : debugOverlaps)
		{
			auto extLen = _seqContainer.seqLen(ovlp.extId);
			Logger::get().debug() << "\t" 
					<< _seqContainer.getIndex().at(ovlp.extId).description
					<< "\tcs:" << ovlp.curBegin << "\tcl:" << ovlp.curRange()
					<< "\tes:" << ovlp.extBegin << "\tel:" << ovlp.extRange()
					<< "\t" << this->overlapTest(ovlp, curLen, extLen);
		}
		_logMutex.unlock();
	}
	*/

	return detectedOverlaps;
}

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
}
