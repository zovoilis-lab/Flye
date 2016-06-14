//(c) 2016 by Authors
//This file is a part of ABruijn program.
//Released under the BSD license (see LICENSE file)

#include <limits>
#include <cassert>
#include <algorithm>
#include <queue>
#include <iomanip>

#include "logger.h"
#include "extender.h"

namespace
{
	template <class T>
	T median(const std::vector<T>& vals)
	{
		if (vals.empty()) return T();
		std::vector<T> tmp = vals;
		std::sort(tmp.begin(), tmp.end());
		//NOTE: there's a bug in libstdc++ nth_element, 
		//that sometimes leads to a segfault
		//std::nth_element(tmp.begin(), tmp.begin() + tmp.size() / 2, 
		//				 tmp.end());
		return tmp[tmp.size() / 2];
	}

	template <class T>
	T mean(const std::vector<T>& vals)
	{
		if (vals.empty()) return T();
		T sum = 0;
		for (const T& val : vals) sum += val;
		return sum / vals.size();
	}
	
	template <class T>
	T vecMax(const std::vector<T>& vals)
	{
		if (vals.empty()) return T();
		T max = 0;
		for (const T& val : vals) max = std::max(max, val);
		return max;
	}
}

ContigPath Extender::extendRead(FastaRecord::Id startRead)
{
	ContigPath contigPath;
	FastaRecord::Id curRead = startRead;
	contigPath.reads.push_back(curRead);
	_visitedReads.insert(curRead);
	_visitedReads.insert(curRead.rc());
	bool rightExtension = true;

	Logger::get().debug() << "Start Read: " << 
				_seqContainer.seqName(startRead);

	std::unordered_set<FastaRecord::Id> curPathVisited;

	while(true)
	{
		FastaRecord::Id extRead = this->stepRight(curRead, startRead);

		if (extRead == startRead && rightExtension)	//circular
		{
			Logger::get().debug() << "Circular contig";
			contigPath.circular = true;
			break;
		}

		if (extRead != FastaRecord::ID_NONE) 
		{
			Logger::get().debug() << "Extension: " << 
				    	_seqContainer.seqName(extRead);
			if (curPathVisited.count(extRead)) 
				Logger::get().debug() << "Visited in this path";
			if (_visitedReads.count(extRead)) 
				Logger::get().debug() << "Visited globally";
		}
		else
		{
			Logger::get().debug() << "No extension found"; 
		}

		if (_visitedReads.count(extRead) || extRead == FastaRecord::ID_NONE)
		{
			if (rightExtension)
			{
				Logger::get().debug() << "Changing direction";
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
				Logger::get().debug() << "Linear contig";
				break;
			}
		}

		_visitedReads.insert(extRead);
		_visitedReads.insert(extRead.rc());
		curPathVisited.insert(extRead);
		curPathVisited.insert(extRead.rc());

		curRead = extRead;
		contigPath.reads.push_back(curRead);
	}

	Logger::get().debug() << "Made " << contigPath.reads.size() - 1 
						  << " extensions";
	return contigPath;
}

int Extender::rightMultiplicity(FastaRecord::Id readId)
{
	std::unordered_set<FastaRecord::Id> extensions;
	auto& overlaps = _ovlpDetector.getOverlapIndex().at(readId);
	for (auto& ovlp : overlaps)
	{
		if (this->isProperRightExtension(ovlp)) 
		{
			extensions.insert(ovlp.extId);
		}
	}
	if (extensions.size() < 2) return 1;

	/*
	for (auto& extCandidate : extensions)
	{
		Logger::get().debug() << "Read " << _seqContainer.seqName(extCandidate);
		for (auto& ovlp : _ovlpDetector.getOverlapIndex().at(extCandidate))
		{
			if (extensions.count(ovlp.extId) &&
		   		ovlp.rightShift < 0)
			{
				Logger::get().debug() << "\t" 
							<< _seqContainer.seqName(ovlp.extId);
			}
		}
	}
	*/
	
	std::unordered_map<FastaRecord::Id, int> clusters;
	std::unordered_set<FastaRecord::Id> coveredReads;
	while(true)
	{
		FastaRecord::Id maxUniqueCoveredId = FastaRecord::ID_NONE;
		int maxUniqueCovered = 0;
		for (auto& extCandidate : extensions)
		{
			//ignore reads that are already belong to clusters
			if (coveredReads.count(extCandidate)) continue;

			//count number of reads the candidate extends
			int extUniqueCovered = 0;
			for (auto& ovlp : _ovlpDetector.getOverlapIndex().at(extCandidate))
			{
				if (extensions.count(ovlp.extId) && ovlp.rightShift < 0 && 
					!coveredReads.count(ovlp.extId))
				{
					++extUniqueCovered;
				}
			}

			//keeping maximum
			if (extUniqueCovered > maxUniqueCovered)
			{
				maxUniqueCovered = extUniqueCovered;
				maxUniqueCoveredId = extCandidate;
			}
		}

		if (maxUniqueCoveredId == FastaRecord::ID_NONE) break;

		clusters[maxUniqueCoveredId] = maxUniqueCovered;
		for (auto& ovlp : _ovlpDetector.getOverlapIndex().at(maxUniqueCoveredId))
		{
			if (extensions.count(ovlp.extId) && ovlp.rightShift < 0)
			{
				coveredReads.insert(ovlp.extId);
			}
		}
	}

	int numClusters = 0;
	for (auto& clustHash : clusters)
	{
		if (clustHash.second > 1) ++numClusters;
	}

	/*
	std::string strClusters;
	for (auto& clustHash : clusters) 
		strClusters += std::to_string(clustHash.second) + " ";
	Logger::get().debug() << "Clusters: " << strClusters;
	*/

	return numClusters;
}

void Extender::assembleContigs()
{
	Logger::get().info() << "Extending reads";
	//FIXME: improve contig lengths filter
	std::vector<size_t> contigLengths;
	_visitedReads.clear();

	for (auto& indexPair : _seqContainer.getIndex())
	{	
		if (_visitedReads.count(indexPair.first) ||
			_chimDetector.isChimeric(indexPair.first) ||
			this->countRightExtensions(indexPair.first) == 0 ||
			this->isBranching(indexPair.first)) continue;

		//additionaly, should not overlap with any of visited reads
		bool ovlpsVisited = false;
		for (auto& ovlp : _ovlpDetector.getOverlapIndex().at(indexPair.first))
		{
			if (_visitedReads.count(ovlp.extId)) ovlpsVisited = true;
		}
		if (ovlpsVisited) continue;

		ContigPath path = this->extendRead(indexPair.first);

		if (contigLengths.empty() || 
			vecMax(contigLengths) / 100 < path.reads.size())
		{
			contigLengths.push_back(path.reads.size());

			//marking visited reads
			std::unordered_set<FastaRecord::Id> leftSupported;
			std::unordered_set<FastaRecord::Id> rightSupported;
			for (auto& readId : path.reads)
			{
				for (auto& ovlp : _ovlpDetector.getOverlapIndex().at(readId))
				{
					if (this->isProperRightExtension(ovlp))
					{
						leftSupported.insert(ovlp.extId);
						rightSupported.insert(ovlp.extId.rc());
					}
					if (this->isProperLeftExtension(ovlp))
					{
						rightSupported.insert(ovlp.extId);
						leftSupported.insert(ovlp.extId.rc());
					}
				}
			}
			for (auto readId : leftSupported)
			{
				if (rightSupported.count(readId))
				{
					_visitedReads.insert(readId);
					_visitedReads.insert(readId.rc());
				}
			}
			
			_contigPaths.push_back(std::move(path));
		}
	}

	Logger::get().info() << "Assembled " << _contigPaths.size() << " contigs";
}


float Extender::extensionIndex(FastaRecord::Id readId)
{
	std::unordered_set<FastaRecord::Id> extensions;
	auto& overlaps = _ovlpDetector.getOverlapIndex().at(readId);

	for (auto& ovlp : overlaps)
	{
		if (this->isProperRightExtension(ovlp) &&
			this->countRightExtensions(ovlp.extId) > 0) 
		{
			extensions.insert(ovlp.extId);
		}
	}
	//if (extensions.size() < 2) return 0.0f;
	
	int maxCovered = 0;
	std::unordered_set<FastaRecord::Id> covered;
	for (auto& extCandidate : extensions)
	{
		int extCovered = 0;
		for (auto& ovlp : _ovlpDetector.getOverlapIndex().at(extCandidate))
		{
			if (extensions.count(ovlp.extId) &&
		   		ovlp.rightShift < 0)
			{
				++extCovered;
				covered.insert(ovlp.extId);
			}
		}
		maxCovered = std::max(maxCovered, extCovered);
	}

	//int nonCovered = extensions.size() - covered.size();
	//return (float)maxCovered / (extensions.size() - 1);
	return !covered.empty() ? (float)maxCovered / covered.size() : 0.0f;
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
		if (this->isProperRightExtension(ovlp)) 
		{
			if (ovlp.extId == startReadId) return startReadId;
			if (this->countRightExtensions(ovlp.extId) > 0)
			{
				extensions.insert(ovlp.extId);
			}
		}
	}
	
	//rank extension candidates
	std::unordered_map<FastaRecord::Id, 
					   std::tuple<int, int, int, int>> supportIndex;

	Logger::get().debug() << "Extension index " << this->extensionIndex(readId);
	//Logger::get().debug() << "Multiplicity " << this->rightMultiplicity(readId);
	//this->rightMultiplicity(readId);

	for (auto& extCandidate : extensions)
	{
		int leftSupport = 0;
		int rightSupport = 0;
		int ovlpSize = 0;
		int ovlpShift = 0;
		for (auto& ovlp : _ovlpDetector.getOverlapIndex().at(extCandidate))
		{
			if (ovlp.extId == readId) ovlpSize = ovlp.curRange();
			if (ovlp.extId == readId) ovlpShift = -ovlp.rightShift;
			if (!extensions.count(ovlp.extId)) continue;

			auto revOvlp = ovlp;
			revOvlp.reverse();

			if (this->isProperRightExtension(ovlp)) ++rightSupport;
			if (this->isProperRightExtension(revOvlp)) ++leftSupport;
			//if (this->isProperLeftExtension(ovlp)) ++leftSupport;
		}
		int minSupport = std::min(leftSupport, rightSupport);
		//int endsRepeat = 1 - this->isBranching(extCandidate);
		//endsRepeat = 1;
		if (!this->isBranching(readId))
		{
			supportIndex[extCandidate] = std::make_tuple(1, minSupport, 
													 	 rightSupport, ovlpSize);
		}
		else
		{
			int startsRepeat = 1 - this->isBranching(extCandidate.rc());
			supportIndex[extCandidate] = std::make_tuple(startsRepeat, minSupport, 
														 rightSupport, ovlpSize);
		}
		Logger::get().debug() << "\t" 
				    << _seqContainer.seqName(extCandidate)
					<< "\t" << leftSupport << "\t" << rightSupport << "\t"
					<< std::fixed << std::setprecision(2) 
					<< this->extensionIndex(extCandidate.rc()) << "\t" << ovlpSize
					<< "\t" << ovlpShift;
	}

	auto bestSupport = std::make_tuple(0, 0, 0, 0);
	auto bestExtension = FastaRecord::ID_NONE;
	for (auto& extCandidate : extensions)
	{
		if (extCandidate == startReadId) return startReadId;

		if (supportIndex[extCandidate] > bestSupport)
		{
			bestSupport = supportIndex[extCandidate];
			bestExtension = extCandidate;
		}
	}


	if (bestExtension != FastaRecord::ID_NONE)
	{
		if (_chimDetector.isChimeric(bestExtension))
			Logger::get().debug() << "Chimeric extension! ";
		if (this->isBranching(bestExtension))
			Logger::get().debug() << "Branching extension! ";
	}
	return bestExtension;
}

int Extender::countRightExtensions(FastaRecord::Id readId)
{
	int count = 0;
	for (auto& ovlp : _ovlpDetector.getOverlapIndex().at(readId))
	{
		if (this->isProperRightExtension(ovlp)) ++count;
	}
	return count;
}

bool Extender::isBranching(FastaRecord::Id readId)
{
	//return this->branchIndex(readId) > 2.0f;
	return this->extensionIndex(readId) <= 0.75f;
	//return this->rightMultiplicity(readId) > 1;
}

//Checks if read is extended to the right
bool Extender::isProperRightExtension(const OverlapRange& ovlp)
{
	return !_chimDetector.isChimeric(ovlp.extId) && 
		   ovlp.rightShift > _maximumJump;
}

//Checks if read is extended to the left
bool Extender::isProperLeftExtension(const OverlapRange& ovlp)
{
	return !_chimDetector.isChimeric(ovlp.extId) &&
		   ovlp.leftShift < -_maximumJump;
}
