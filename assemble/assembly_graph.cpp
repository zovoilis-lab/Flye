
#include "assembly_graph.h"
#include "config.h"
#include "logger.h"
#include <set>

void AssemblyGraph::construct()
{
	Logger::get().info() << "Overlapping";
	for (auto readId : _seqContainer.getIndex())
	{
		//if (readId.first.strand())
		{
			auto overlaps = this->getReadOverlaps(readId.first);

			for (auto ovlp : overlaps)
			{
				_overlapIndex[ovlp.curId].push_back(ovlp);
			}
		}
	}
}

AssemblyGraph::JumpRes 
AssemblyGraph::jumpTest(int32_t curPrev, int32_t curNext,
						  int32_t extPrev, int32_t extNext) const
{
	if (curNext - curPrev > Constants::maximumJump) return J_END;

	const int JMP = 100;
	if (0 < curNext - curPrev && curNext - curPrev < JMP &&
		0 < extNext - extPrev && extNext - extPrev < JMP)
	{
		if (abs((curNext - curPrev) - (extNext - extPrev)) 
			< Constants::maximumJump / Constants::closeJumpRate)
		{
			return J_CLOSE;
		}
		if (abs((curNext - curPrev) - (extNext - extPrev)) 
			< JMP / Constants::farJumpRate)
		{
			return J_FAR;
		}
	}
	return J_INCONS;
}

void AssemblyGraph::saveOverlaps(const std::string& filename)
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

//Getting all possible overlaps
//based on the shared kmers (common jump-paths)
std::vector<OverlapRange> 
AssemblyGraph::getReadOverlaps(FastaRecord::Id currentReadId) const
{
	std::unordered_map<FastaRecord::Id, 
					   std::vector<OverlapRange>> activePaths;
	std::set<size_t> eraseMarks;
	size_t curLen = _seqContainer.seqLen(currentReadId);
	std::vector<KmerPosition> solidKmersCache;

	//for all kmers in this read
	for (auto curKmerPos : IterSolidKmers(currentReadId))
	{
		int32_t curPos = curKmerPos.position;
		solidKmersCache.push_back(curKmerPos);

		//for all other occurences of this kmer (extension candidates)
		for (const auto& extReadPos : _vertexIndex.byKmer(curKmerPos.kmer))
		{
			if (_vertexIndex.isRepetitive(curKmerPos.kmer) &&
				!activePaths.count(extReadPos.readId))
			{
				continue;
			}

			if (extReadPos.readId == currentReadId &&
				extReadPos.position == curPos) continue;

			//size_t extLen = _seqContainer.seqLen(extReadPos.readId);
			//if (extLen < (size_t)Parameters::minimumOverlap) continue;

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
			if (!extendsClose && !extendsFar)
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
				detectedOverlaps.push_back(ovlp);
			}
			//if (debugOverlap.curRange() < ovlp.curRange()) debugOverlap = ovlp;
		}

		if (passedTest)
		{
			//detectedOverlaps.push_back(maxOverlap);
		}
		//if (debugOverlap.curRange() > 3000) debugOverlaps.push_back(debugOverlap);
	}
	
	//if (!debugOverlaps.empty())
	//{
		//_logMutex.lock();
		Logger::get().debug() << "Ovlps for " 
					<< _seqContainer.seqName(currentReadId);
		for (auto& ovlp : detectedOverlaps)
		{
			auto extLen = _seqContainer.seqLen(ovlp.extId);
			auto curBegin = ovlp.curId.strand() ? ovlp.curBegin : 
							curLen - ovlp.curBegin - ovlp.curRange();
			auto extBegin = ovlp.extId.strand() ? ovlp.extBegin : 
							extLen - ovlp.extBegin - ovlp.extRange();
			Logger::get().debug() << "\t" 
					<< _seqContainer.getIndex().at(ovlp.extId).description
					<< "\tcs:" << curBegin << "\tcl:" << ovlp.curRange()
					<< "\tes:" << extBegin << "\tel:" << ovlp.extRange()
					<< "\t" << this->overlapTest(ovlp, curLen, extLen);
		}
		//_logMutex.unlock();
	//}

	return detectedOverlaps;
}

//Check if it is a proper overlap
bool AssemblyGraph::overlapTest(const OverlapRange& ovlp, int32_t curLen, 
								  int32_t extLen) const
{
	if (ovlp.curRange() < Parameters::minimumOverlap || 
		ovlp.extRange() < Parameters::minimumOverlap) 
	{
		return false;
	}

	float lengthDiff = abs(ovlp.curRange() - ovlp.extRange());
	float meanLength = (ovlp.curRange() + ovlp.extRange()) / 2.0f;
	if (lengthDiff > meanLength / Constants::overlapDivergenceRate)
	{
		return false;
	}

	return true;
}

