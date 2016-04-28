//(c) 2016 by Authors
//This file is a part of ABruijn program.
//Released under the BSD license (see LICENSE file)

#include <limits>
#include <cassert>
#include <algorithm>

#include "utility.h"
#include "extender.h"

ContigPath Extender::extendRead(FastaRecord::ReadIdType startRead)
{
	ContigPath contigPath;
	FastaRecord::ReadIdType curRead = startRead;
	contigPath.reads.push_back(curRead);
	_visitedReads.insert(curRead);
	_visitedReads.insert(-curRead);
	bool rightExtension = true;

	DEBUG_PRINT("Start Read: " << 
				_seqContainer.getIndex().at(startRead).description);

	while(true)
	{
		FastaRecord::ReadIdType extRead = this->stepRight(curRead, startRead);

		if (extRead == startRead)	//circular
		{
			DEBUG_PRINT("Circular contig");
			contigPath.circular = true;
			break;
		}

		if (_visitedReads.count(extRead))	//loop
		{
			LOG_PRINT("Looped contig");
			break;
		}

		if (extRead == FastaRecord::ID_NONE)	//dead end
		{
			if (rightExtension)
			{
				DEBUG_PRINT("Changing direction");
				break;
				rightExtension = false;
				extRead = -startRead;
				std::reverse(contigPath.reads.begin(), contigPath.reads.end());
				for (size_t i = 0; i < contigPath.reads.size(); ++i) 
				{
					contigPath.reads[i] = -contigPath.reads[i];
				}
				contigPath.reads.pop_back();
			}
			else
			{
				DEBUG_PRINT("Linear contig");
				break;
			}
		}

		DEBUG_PRINT("Extension: " << 
				    _seqContainer.getIndex().at(extRead).description);

		_visitedReads.insert(extRead);
		_visitedReads.insert(-extRead);
		curRead = extRead;
		contigPath.reads.push_back(curRead);

	}
	LOG_PRINT("Made " << contigPath.reads.size() << " extensions");
	return contigPath;
}

void Extender::assembleContigs()
{
	LOG_PRINT("Extending reads");
	_visitedReads.clear();

	while (true)
	{
		//choose a read for extension
		int maxExtension = 0;
		FastaRecord::ReadIdType startRead = FastaRecord::ID_NONE;
		for (auto& indexPair : _seqContainer.getIndex())
		{	
			if (_visitedReads.count(indexPair.first) ||
				_chimDetector.isChimeric(indexPair.first)) continue;

			/*
			bool overlapsVisited = false;
			for (auto& ovlp : _ovlpDetector.getOverlapIndex().at(indexPair.first))
			{
				if (_visitedReads.count(ovlp.extId))
				{
					overlapsVisited = true;
					break;
				}
			}
			if (overlapsVisited) continue;
			*/

			if (this->countRightExtensions(indexPair.first) > maxExtension)
			{
				maxExtension = this->countRightExtensions(indexPair.first);
				startRead = indexPair.first;
			}
		}
		if (startRead == FastaRecord::ID_NONE) break;

		_contigPaths.push_back(this->extendRead(startRead));
		for (auto& readId : _contigPaths.back().reads)
		{
			for (auto& ovlp : _ovlpDetector.getOverlapIndex().at(readId))
			{
				if (this->branchIndex(ovlp.extId) > 0.5)
				{
					_visitedReads.insert(ovlp.extId);
					_visitedReads.insert(-ovlp.extId);
				}
			}
		}
	}
}

float Extender::branchIndex(FastaRecord::ReadIdType readId)
{
	auto& overlaps = _ovlpDetector.getOverlapIndex().at(readId);
	std::unordered_set<FastaRecord::ReadIdType> extensions;
	for (auto& ovlp : overlaps)
	{
		if (this->isProperRightExtension(ovlp) &&
			!_chimDetector.isChimeric(ovlp.extId))
		{
			extensions.insert(ovlp.extId);
		}
	}

	std::vector<int> ovlpIndices;
	for (auto& ovlp : overlaps)
	{
		if (!extensions.count(ovlp.extId)) continue;

		int ovlpIndex = 0;
		auto& extOverlaps = _ovlpDetector.getOverlapIndex().at(ovlp.extId);
		for (auto& extOvlp : extOverlaps)
		{
			if (extensions.count(extOvlp.extId)) ++ovlpIndex;
		}
		ovlpIndices.push_back(ovlpIndex);
	}
	if (extensions.empty()) return 0.0f;

	float total = 0;
	for (int ovlpIndex : ovlpIndices)
	{
		total += ((float)ovlpIndex + 1) / extensions.size();
	}
	float ovlpIndex = total / ovlpIndices.size();
	return ovlpIndex;
}

//makes one extension to the right
FastaRecord::ReadIdType Extender::stepRight(FastaRecord::ReadIdType readId, 
										    FastaRecord::ReadIdType startReadId)
{
	auto& overlaps = _ovlpDetector.getOverlapIndex().at(readId);

	int32_t maxOverlap = std::numeric_limits<int32_t>::min();
	FastaRecord::ReadIdType bestExtension = FastaRecord::ID_NONE;
	bool reliableExtension = false;

	for (auto& ovlp : overlaps)
	{
		assert(ovlp.curId != ovlp.extId);
		if (!this->isProperRightExtension(ovlp))
			continue;
		if (ovlp.extId == startReadId)
			return startReadId;
		if (_visitedReads.count(ovlp.extId))
			continue;
		
		if (!_chimDetector.isChimeric(ovlp.extId) &&
			this->branchIndex(ovlp.extId) > 0.5 &&
			this->countRightExtensions(ovlp.extId) > 0)
		{
			if (!reliableExtension || maxOverlap < ovlp.curRange())
			{
				maxOverlap = ovlp.curRange();
				bestExtension = ovlp.extId;
			}
			reliableExtension = true;
		}
		else
		{
			if (!reliableExtension && maxOverlap < ovlp.curRange())
			{
				maxOverlap = ovlp.curRange();
				bestExtension = ovlp.extId;
			}
		}
	}

	if (bestExtension != FastaRecord::ID_NONE)
	{
		float ovlpIndex = this->branchIndex(bestExtension);
		if (ovlpIndex < 0.5) 
			DEBUG_PRINT("Making branching extension: " << ovlpIndex);
	}
	if (!reliableExtension) DEBUG_PRINT("Making non-reliable extension!");

	return bestExtension;
}


int Extender::countRightExtensions(FastaRecord::ReadIdType readId)
{
	int count = 0;
	for (auto& ovlp : _ovlpDetector.getOverlapIndex().at(readId))
	{
		if (this->isProperRightExtension(ovlp))
			++count;
	}
	return count;
}


//Checks if read is extended to the right
bool Extender::isProperRightExtension(const OverlapDetector::OverlapRange& ovlp)
{
	int32_t curLen = _seqContainer.getIndex().at(ovlp.curId).sequence.length();
	int32_t extLen = _seqContainer.getIndex().at(ovlp.extId).sequence.length();
	return extLen - ovlp.extEnd > curLen - ovlp.curEnd;
}
