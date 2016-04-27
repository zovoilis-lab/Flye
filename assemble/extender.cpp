//(c) 2016 by Authors
//This file is a part of ABruijn program.
//Released under the BSD license (see LICENSE file)

#include <limits>
#include <cassert>

#include "utility.h"
#include "extender.h"

void Extender::extendReads()
{
	//TODO: multiple chromosome support
	//TODO: handle looped structures
	LOG_PRINT("Extending reads");
	FastaRecord::ReadIdType startRead = FastaRecord::ID_NONE;

	int maxExtensions = std::numeric_limits<int>::min();
	for (auto& indexPair : _seqContainer.getIndex())
	{
		if (!_chimDetector.isChimeric(indexPair.first) &&
			this->countRightExtensions(indexPair.first) > maxExtensions)
		{
			maxExtensions = this->countRightExtensions(indexPair.first);
			startRead = indexPair.first;
		}
		//if (indexPair.second.description == "+seq_19095/0_27206")
		//	startRead = indexPair.first;
	}
	if (startRead == FastaRecord::ID_NONE) 
		throw std::runtime_error("No non-chimeric extensions found");

	FastaRecord::ReadIdType curRead = startRead;
	ReadPath curPath;
	curPath.push_back(curRead);
	std::unordered_set<FastaRecord::ReadIdType> visited;
	while(true)
	{
		FastaRecord::ReadIdType extRead = this->stepRight(curRead, startRead);

		if (extRead == startRead)	//all good
			break;
		if (extRead == FastaRecord::ID_NONE)
		{
			WARNING_PRINT("No further extension found");
			break;
		}
		if (visited.count(extRead))
		{
			WARNING_PRINT("Looped structure while extending");
			break;
		}
		visited.insert(extRead);
		//DEBUG_PRINT("Extension: " << 
		//		    _seqContainer.getIndex().at(extRead).description);

		curRead = extRead;
		curPath.push_back(curRead);
	}
	_contigPaths.push_back(curPath);
}

float Extender::isBranching(FastaRecord::ReadIdType readId)
{
	auto& overlaps = _ovlpDetector.getOverlapIndex().at(readId);
	std::unordered_set<FastaRecord::ReadIdType> extensions;
	for (auto& ovlp : overlaps)
	{
		if (this->isProperRightExtension(ovlp))
		{
			extensions.insert(ovlp.extId);
		}
	}

	std::vector<int> ovlpIndices;
	for (auto& ovlp : overlaps)
	{
		if (!this->isProperRightExtension(ovlp)) continue;

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
		
		if (!_chimDetector.isChimeric(ovlp.extId) &&
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
	if (!reliableExtension)
	{
		DEBUG_PRINT("Making non-reliable extension!");
	}
	if (bestExtension != FastaRecord::ID_NONE)
	{
		float ovlpIndex = this->isBranching(bestExtension);
		if (ovlpIndex < 0.5) 
			DEBUG_PRINT("Making branching extension: " << ovlpIndex);
	}

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
