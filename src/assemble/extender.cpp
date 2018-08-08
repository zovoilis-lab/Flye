//(c) 2016 by Authors
//This file is a part of ABruijn program.
//Released under the BSD license (see LICENSE file)

#include <limits>
#include <algorithm>
#include <iomanip>
#include <stack>

#include "../common/config.h"
#include "../common/logger.h"
#include "../common/parallel.h"
#include "extender.h"


Extender::ExtensionInfo Extender::extendContig(FastaRecord::Id startRead)
{

	std::unordered_set<FastaRecord::Id> currentReads;
	currentReads.insert(startRead);
	currentReads.insert(startRead.rc());

	bool rightExtension = true;
	FastaRecord::Id currentRead = startRead;
	std::vector<int> numExtensions;
	std::vector<int> overlapSizes;
	ExtensionInfo exInfo;
	exInfo.reads.push_back(startRead);
	exInfo.assembledLength = _readsContainer.seqLen(startRead);

	auto leftExtendsStart = [startRead, this](const FastaRecord::Id readId)
	{
		for (auto& ovlp : _ovlpContainer.lazySeqOverlaps(readId))
		{
			if (ovlp.extId == startRead &&
				ovlp.leftShift > (int)Config::get("maximum_jump")) 
			{
				return true;
			}
		}
		return false;
	};

	while(true)
	{
		auto& overlaps = _ovlpContainer.lazySeqOverlaps(currentRead);
		std::vector<OverlapRange> extensions;
		//int innerOverlaps = 0;
		for (auto& ovlp : overlaps)
		{
			//if (_innerReads.contains(ovlp.extId)) ++innerOverlaps;
			if (this->extendsRight(ovlp)) extensions.push_back(ovlp);
		}
		numExtensions.push_back(extensions.size());

		//sort from longes to shortest overlap
		std::sort(extensions.begin(), extensions.end(), 
				  [](const OverlapRange& a, const OverlapRange& b)
					 {return a.curRange() > b.curRange();});

		//checking if read overlaps with one of the used reads
		bool foundExtension = false;

		//getting extension
		int minExtensions = std::max(1, (int)extensions.size() / 
										(int)Config::get("max_coverage_drop_rate"));

		/*Logger::get().debug() << _readsContainer.seqName(currentRead) 
			<< "\t" << overlaps.size();
		for (auto& ovlp : overlaps)
		{
			if (this->extendsRight(ovlp))
			{
				Logger::get().debug() << "\t" << 
					_readsContainer.seqName(ovlp.extId) << "\t" 
					<< "o:" << ovlp.curRange() << "\tchim:"
					<< _chimDetector.isChimeric(ovlp.extId)
					<< "\trep:" << this->isRightRepeat(ovlp.extId)
					<< "\text:" << this->countRightExtensions(ovlp.extId);
			}
		}*/

		OverlapRange* maxExtension = nullptr;
		for (auto& ovlp : extensions)
		{
			if(leftExtendsStart(ovlp.extId)) continue;

			//try to find a good one
			if (!_chimDetector.isChimeric(ovlp.extId) &&
				this->countRightExtensions(ovlp.extId) > minExtensions)
			{
				foundExtension = true;
				maxExtension = &ovlp;
				exInfo.assembledLength += ovlp.rightShift;
				currentRead = ovlp.extId;
				break;
			}

			if (!maxExtension || maxExtension->rightShift < ovlp.rightShift)
			{
				maxExtension = &ovlp;
			}
		}
		//in case of suspicious extension make the farthest jump possible
		if (!foundExtension && maxExtension)
		{
			++exInfo.numSuspicious;
			exInfo.assembledLength += maxExtension->rightShift;
			foundExtension = true;
			currentRead = maxExtension->extId;
		}

		bool overlapsVisited = _innerReads.contains(currentRead);
		if (foundExtension) 
		{
			exInfo.reads.push_back(currentRead);
			overlapSizes.push_back(maxExtension->curRange());
			overlapsVisited |= currentReads.count(currentRead);
		}
		else
		{
			rightExtension ? exInfo.leftTip = true : exInfo.rightTip = true;
		}

		if (!foundExtension || overlapsVisited)
		{
			//Logger::get().debug() << "Not found: " << !foundExtension << 
			//	" overlaps visited: " << overlapsVisited;

			//right extension done, try to extend left from start read
			if (rightExtension && !exInfo.reads.empty())
			{
				exInfo.stepsToTurn = exInfo.reads.size();
				rightExtension = false;
				currentRead = exInfo.reads.front().rc();
				std::reverse(exInfo.reads.begin(), exInfo.reads.end());
				for (size_t i = 0; i < exInfo.reads.size(); ++i) 
				{
					exInfo.reads[i] = exInfo.reads[i].rc();
				}
			}
			//done with extension
			else
			{
				break;
			}
		}

		currentReads.insert(currentRead);
		currentReads.insert(currentRead.rc());
	}

	if (!numExtensions.empty())
	{
		exInfo.meanOverlaps = median(numExtensions);
	}
	if (!overlapSizes.empty())
	{
		exInfo.avgOverlapSize = median(overlapSizes);
		exInfo.minOverlapSize = *std::min_element(overlapSizes.begin(), 
												  overlapSizes.end());
	}

	return exInfo;
}


void Extender::assembleContigs()
{
	static const int MAX_JUMP = Config::get("maximum_jump");
	Logger::get().info() << "Extending reads";
	_chimDetector.estimateGlobalCoverage();
	_ovlpContainer.overlapDivergenceStats();
	_innerReads.clear();
	cuckoohash_map<FastaRecord::Id, size_t> coveredReads;
	
	std::mutex indexMutex;
	auto processRead = [this, &indexMutex, &coveredReads] 
		(FastaRecord::Id startRead)
	{
		if (coveredReads.contains(startRead)) return true;

		int numInnerOvlp = 0;
		for (auto& ovlp : _ovlpContainer.lazySeqOverlaps(startRead))
		{
			if (_innerReads.contains(ovlp.extId)) ++numInnerOvlp;
			if (ovlp.extId == startRead || 
				ovlp.extId.rc() == startRead) return true;
		}
		if (numInnerOvlp > 0) return true;

		int maxStartExt = _chimDetector.getOverlapCoverage() * 10;
		int minStartExt = std::max(1, _chimDetector.getOverlapCoverage() / 50);
		int extLeft = this->countRightExtensions(startRead.rc());
		int extRight = this->countRightExtensions(startRead);
		if (_chimDetector.isChimeric(startRead) ||
			std::max(extLeft, extRight) > maxStartExt ||
			std::min(extLeft, extRight) < minStartExt) return true;
		
		//Good to go!
		ExtensionInfo exInfo = this->extendContig(startRead);
		if (exInfo.reads.size() < 
			(size_t)Config::get("min_reads_in_contig")) return true;

		//Exclusive part - updating the overall assembly
		std::lock_guard<std::mutex> guard(indexMutex);
		
		int innerCount = 0;
		for (auto& readId : exInfo.reads)
		{
			if (_innerReads.contains(readId)) ++innerCount;
		}
		int innerThreshold = std::min((int)Config::get("max_inner_reads"),
									  int((float)Config::get("max_inner_fraction") * 
										  exInfo.reads.size()));
		if (innerCount > innerThreshold)
		{
			Logger::get().debug() << "Discarded contig with "
				<< exInfo.reads.size() << " reads and "
				<< innerCount << " inner overlaps";
			return false;
		}

		Logger::get().debug() << "Assembled contig " 
			<< std::to_string(_readLists.size() + 1)
			<< "\n\tWith " << exInfo.reads.size() << " reads"
			<< "\n\tStart read: " << _readsContainer.seqName(startRead)
			<< "\n\tAt position: " << exInfo.stepsToTurn
			<< "\n\tleftTip: " << exInfo.leftTip 
			<< " rightTip: " << exInfo.rightTip
			<< "\n\tSuspicios: " << exInfo.numSuspicious
			<< "\n\tMean extensions: " << exInfo.meanOverlaps
			<< "\n\tAvg overlap len: " << exInfo.avgOverlapSize
			<< "\n\tMin overlap len: " << exInfo.minOverlapSize
			<< "\n\tInner reads: " << innerCount
			<< "\n\tLength: " << exInfo.assembledLength;
		
		//update inner read index
		std::unordered_set<FastaRecord::Id> rightExtended;
		std::unordered_set<FastaRecord::Id> leftExtended;
		for (auto& readId : exInfo.reads)
		{
			coveredReads.insert(readId, true);
			coveredReads.insert(readId.rc(), true);
			_innerReads.insert(readId, true);
			_innerReads.insert(readId.rc(), true);

			//repetitive read - will bring to many "off-target" reads
			int maxExtensions = exInfo.meanOverlaps * 2;
			if (this->countRightExtensions(readId) > maxExtensions||
				this->countRightExtensions(readId.rc()) > maxExtensions) continue;

			//so each read is covered from the left and right
			for (auto& ovlp : _ovlpContainer.lazySeqOverlaps(readId))
			{
				coveredReads.insert(ovlp.extId, true);
				coveredReads.insert(ovlp.extId.rc(), true);

				//contained inside the other read - probably repetitive
				//if (ovlp.leftShift > MAX_JUMP &&
				//	ovlp.rightShift < -MAX_JUMP) continue;

				if (ovlp.leftShift > MAX_JUMP)
				{
					leftExtended.insert(ovlp.extId);
					rightExtended.insert(ovlp.extId.rc());
				}
				if (ovlp.rightShift < -MAX_JUMP)	
				{
					rightExtended.insert(ovlp.extId);
					leftExtended.insert(ovlp.extId.rc());
				}
			}
		}
		for (auto& read : rightExtended)
		{
			if (leftExtended.count(read)) 
			{
				_innerReads.insert(read, true);
				_innerReads.insert(read.rc(), true);
			}
		}

		Logger::get().debug() << "Inner: " << 
			_innerReads.size() << " covered: " << coveredReads.size()
			<< " total: "<< _readsContainer.iterSeqs().size();
		
		_readLists.push_back(std::move(exInfo));
		return true;
	};

	std::function<void(const FastaRecord::Id&)> threadWorker = 
		[processRead] (const FastaRecord::Id& readId)
	{
		processRead(readId);
	};
	std::vector<FastaRecord::Id> allReads;
	for (auto& seq : _readsContainer.iterSeqs())
	{
		allReads.push_back(seq.id);
	}
	std::random_shuffle(allReads.begin(), allReads.end());
	processInParallel(allReads, threadWorker,
					  Parameters::get().numThreads, true);
	_ovlpContainer.ensureTransitivity(/*only max*/ true);

	bool addSingletons = (bool)Config::get("add_unassembled_reads");
	if (addSingletons)
	{
		std::vector<FastaRecord::Id> sortedByLength;
		for (auto& seq : _readsContainer.iterSeqs())
		{
			if (seq.id.strand() && !_innerReads.contains(seq.id) &&
				_readsContainer.seqLen(seq.id) > Parameters::get().minimumOverlap)
			{
				sortedByLength.push_back(seq.id);
			}
		}
		std::sort(sortedByLength.begin(), sortedByLength.end(),
				  [this](FastaRecord::Id idOne, FastaRecord::Id idTwo)
				  	{return _readsContainer.seqLen(idOne) > 
							_readsContainer.seqLen(idTwo);});

		int singletonsAdded = 0;
		std::unordered_set<FastaRecord::Id> coveredLocal;
		for (auto readId : sortedByLength)
		{
			if (!coveredLocal.count(readId))
			{
				for (auto& ovlp : _ovlpContainer.lazySeqOverlaps(readId))
				{
					if (ovlp.leftShift >= 0 && ovlp.rightShift <= 0)
					{
						coveredLocal.insert(ovlp.extId);
						coveredLocal.insert(ovlp.extId.rc());
					}
				}
				ExtensionInfo path;
				path.singleton = true;
				path.reads.push_back(readId);
				_readLists.push_back(path);
				++singletonsAdded;
			}
		}
		Logger::get().info() << "Added " << singletonsAdded << " singleton reads";
	}

	this->convertToContigs();
	Logger::get().info() << "Assembled " << _contigPaths.size() << " draft contigs";
}

void Extender::convertToContigs()
{
	for (auto& exInfo : _readLists)
	{
		ContigPath path;
		if (!exInfo.singleton)
		{
			path.name = "contig_" + std::to_string(_contigPaths.size() + 1);
		}
		else
		{
			path.name = "read_" + std::to_string(_contigPaths.size() + 1);
		}

		for (size_t i = 0; i < exInfo.reads.size() - 1; ++i)
		{
			bool found = false;
			OverlapRange readsOvlp;

			for (auto& ovlp : _ovlpContainer.lazySeqOverlaps(exInfo.reads[i]))
			{
				if (ovlp.extId == exInfo.reads[i + 1]) 
				{
					readsOvlp = ovlp;
					found = true;
					break;
				}
			}
			if (!found) throw std::runtime_error("Ovlp not found!");

			path.sequences.push_back(_readsContainer.getSeq(exInfo.reads[i]));
			path.overlaps.push_back(readsOvlp);
		}
		path.sequences.push_back(_readsContainer.getSeq(exInfo.reads.back()));
		_contigPaths.push_back(std::move(path));
	}
}

int Extender::countRightExtensions(FastaRecord::Id readId) const
{
	int count = 0;
	for (auto& ovlp : _ovlpContainer.lazySeqOverlaps(readId))
	{
		if (this->extendsRight(ovlp)) ++count;
	}
	return count;
}

/*bool Extender::isRightRepeat(FastaRecord::Id readId) const
{
	int maxExtensions = _chimDetector.getOverlapCoverage() * 10;
	return this->countRightExtensions(readId) >= maxExtensions;
}*/

bool Extender::extendsRight(const OverlapRange& ovlp) const
{
	static const int MAX_JUMP = Config::get("maximum_jump");
	return ovlp.rightShift > MAX_JUMP;
}
