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
		throw std::runtime_error("Error: No non-chimeric extensions found, exiting");

	FastaRecord::ReadIdType curRead = startRead;
	ReadPath curPath;
	curPath.push_back(curRead);
	std::unordered_set<FastaRecord::ReadIdType> visited;
	while(true)
	{
		FastaRecord::ReadIdType extRead = this->stepRight(curRead, startRead);
		DEBUG_PRINT("Extension: " << 
				    _seqContainer.getIndex().at(extRead).description);

		if (extRead == FastaRecord::ID_NONE)
			WARNING_PRINT("No further extension found");
		if (extRead == startRead || extRead == FastaRecord::ID_NONE)
			break;

		if (visited.count(extRead))
			throw std::runtime_error("Looped structure while extending");
		visited.insert(extRead);

		curRead = extRead;
		curPath.push_back(curRead);
	}
	_contigPaths.push_back(curPath);
}


//makes one extension to the right
FastaRecord::ReadIdType Extender::stepRight(FastaRecord::ReadIdType readId, 
										    FastaRecord::ReadIdType startReadId)
{
	auto& overlaps = _ovlpDetector.getOverlapIndex().at(readId);

	int32_t maxOverlap = std::numeric_limits<int32_t>::min();
	FastaRecord::ReadIdType bestExtension = FastaRecord::ID_NONE;
	bool reliableExtension = true;

	for (auto& ovlp : overlaps)
	{
		assert(ovlp.curId != ovlp.extId);
		if (!this->isProperRightExtension(ovlp))
			continue;
		if (ovlp.extId == startReadId)
			return startReadId;
		
		if (!_chimDetector.isChimeric(ovlp.extId) &&
			this->countRightExtensions(ovlp.extId) > 0)
		//if (!_chimDetector.isChimeric(ovlp.extId))
		{
			//replace non-reliable with reliable
			if (reliableExtension)
			{
				maxOverlap = ovlp.curRange();
				bestExtension = ovlp.extId;
				reliableExtension = false;
				continue;
			}
		}
		else if (!reliableExtension)
		{
			//don't replace reliable with non-reliable
			continue;
		}

		if (maxOverlap < ovlp.curRange())
		{
			maxOverlap = ovlp.curRange();
			bestExtension = ovlp.extId;
		}
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
