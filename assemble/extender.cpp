//(c) 2016 by Authors
//This file is a part of ABruijn program.
//Released under the BSD license (see LICENSE file)

#include <limits>
#include <cassert>
#include <algorithm>

#include "utility.h"
#include "extender.h"

ContigPath Extender::extendRead(FastaRecord::Id startRead)
{
	ContigPath contigPath;
	FastaRecord::Id curRead = startRead;
	contigPath.reads.push_back(curRead);
	_visitedReads.insert(curRead);
	_visitedReads.insert(curRead.rc());
	bool rightExtension = true;

	DEBUG_PRINT("Start Read: " << 
				_seqContainer.getIndex().at(startRead).description);

	std::unordered_set<FastaRecord::Id> curPathVisited;
	while(true)
	{
		FastaRecord::Id extRead = this->stepRight(curRead, startRead);

		if (extRead == startRead)	//circular
		{
			DEBUG_PRINT("Circular contig");
			contigPath.circular = true;
			break;
		}

		if (curPathVisited.count(extRead)) DEBUG_PRINT("Visited in this path");
		if (_visitedReads.count(extRead)) DEBUG_PRINT("Visited globally");
		if (extRead == FastaRecord::ID_NONE) DEBUG_PRINT("No extension found");

		if (_visitedReads.count(extRead) || extRead == FastaRecord::ID_NONE)
		{
			if (rightExtension)
			{
				DEBUG_PRINT("Changing direction");
				rightExtension = false;
				extRead = startRead.rc();
				std::reverse(contigPath.reads.begin(), contigPath.reads.end());
				for (size_t i = 0; i < contigPath.reads.size(); ++i) 
				{
					contigPath.reads[i] = contigPath.reads[i].rc();
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
		_visitedReads.insert(extRead.rc());
		curPathVisited.insert(extRead);
		curPathVisited.insert(extRead.rc());
		curRead = extRead;
		contigPath.reads.push_back(curRead);
	}

	for (auto& readId : contigPath.reads)
	{
		for (auto& ovlp : _ovlpDetector.getOverlapIndex().at(readId))
		{
			_visitedReads.insert(ovlp.extId);
			_visitedReads.insert(ovlp.extId.rc());
		}
	}

	LOG_PRINT("Made " << contigPath.reads.size() - 1 << " extensions");
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
		FastaRecord::Id startRead = FastaRecord::ID_NONE;
		for (auto& indexPair : _seqContainer.getIndex())
		{	
			if (_visitedReads.count(indexPair.first) ||
				_chimDetector.isChimeric(indexPair.first)) continue;

			if (this->countRightExtensions(indexPair.first) > maxExtension)
			{
				maxExtension = this->countRightExtensions(indexPair.first);
				startRead = indexPair.first;
			}
		}
		if (startRead == FastaRecord::ID_NONE) break;

		ContigPath path = this->extendRead(startRead);
		if (path.reads.size() > 2) _contigPaths.push_back(std::move(path));
	}
	LOG_PRINT("Assembled " << _contigPaths.size() << " contigs");
}

float Extender::branchIndex(FastaRecord::Id readId)
{
	auto& overlaps = _ovlpDetector.getOverlapIndex().at(readId);
	std::unordered_set<FastaRecord::Id> extensions;
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
	if (extensions.empty()) return 1.0f;

	float total = 0;
	for (int ovlpIndex : ovlpIndices)
	{
		total += ((float)ovlpIndex + 1) / extensions.size();
	}
	float ovlpIndex = total / ovlpIndices.size();
	return ovlpIndex;
}

//makes one extension to the right
FastaRecord::Id Extender::stepRight(FastaRecord::Id readId, 
									FastaRecord::Id startReadId)
{
	auto& overlaps = _ovlpDetector.getOverlapIndex().at(readId);
	std::unordered_set<FastaRecord::Id> extensions;

	for (auto& ovlp : overlaps)
	{
		assert(ovlp.curId != ovlp.extId);
		if (this->isProperRightExtension(ovlp) &&
			this->countRightExtensions(ovlp.extId) > 0) 
		{
			extensions.insert(ovlp.extId);
		}
	}

	//rank extension candidates
	std::unordered_map<FastaRecord::Id, 
					   std::tuple<int, int, int>> supportIndex;
	for (auto& extCandidate : extensions)
	{
		int leftSupport = 0;
		int rightSupport = 0;
		int ovlpSize = 0;
		for (auto& ovlp : _ovlpDetector.getOverlapIndex().at(extCandidate))
		{
			if (ovlp.extId == readId) ovlpSize = ovlp.curRange();
			if (!extensions.count(ovlp.extId)) continue;

			//auto revOvlp = ovlp;
			//revOvlp.reverse();

			if (this->isProperRightExtension(ovlp)) ++rightSupport;
			if (this->isProperLeftExtension(ovlp)) ++leftSupport;
		}
		int minSupport = std::min(leftSupport, rightSupport);
		supportIndex[extCandidate] = std::make_tuple(minSupport, ovlpSize,
													 rightSupport);
		DEBUG_PRINT("\t" << _seqContainer.getIndex().at(extCandidate).description
					<< " " << leftSupport << " " << rightSupport << " "
					<< this->branchIndex(extCandidate));
	}

	auto bestSupport = std::make_tuple(0, 0, 0);
	auto bestExtension = FastaRecord::ID_NONE;
	for (auto& extCandidate : extensions)
	{
		if (extCandidate == startReadId) return startReadId;
		//if (this->branchIndex(extCandidate) < 0.5) continue;

		if (supportIndex[extCandidate] > bestSupport)
		{
			bestSupport = supportIndex[extCandidate];
			bestExtension = extCandidate;
		}
	}


	if (bestExtension != FastaRecord::ID_NONE)
	{
		//DEBUG_PRINT(this->branchIndex(bestExtension));
		if (std::get<2>(supportIndex[bestExtension]) == 0)
			DEBUG_PRINT("No right support! " << 
					_seqContainer.getIndex().at(bestExtension).description);
		if (_chimDetector.isChimeric(bestExtension))
			DEBUG_PRINT("Chimeric extension! " << 
					_seqContainer.getIndex().at(bestExtension).description);
		if (this->branchIndex(bestExtension) < 0.5)
			DEBUG_PRINT("Branching extension! " << 
						this->branchIndex(bestExtension) << " " << 
						_seqContainer.getIndex().at(bestExtension).description);
	}
	return bestExtension;
}

int Extender::countRightExtensions(FastaRecord::Id readId)
{
	int count = 0;
	for (auto& ovlp : _ovlpDetector.getOverlapIndex().at(readId))
	{

		//if (!_visitedReads.count(ovlp.extId) &&
		//	this->isProperRightExtension(ovlp)) ++count;
		if (this->isProperRightExtension(ovlp)) ++count;
	}
	return count;
}

//Checks if read is extended to the right
bool Extender::isProperRightExtension(const OverlapRange& ovlp)
{
	return ovlp.rightShift > _maximumJump;
}

//Checks if read is extended to the left
bool Extender::isProperLeftExtension(const OverlapRange& ovlp)
{
	return ovlp.leftShift < _maximumJump;
}
