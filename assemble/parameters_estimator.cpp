//(c) 2016 by Authors
//This file is a part of ABruijn program.
//Released under the BSD license (see LICENSE file)

#include "parameters_estimator.h"
#include "utility.h"

int ParametersEstimator::estimateMinKmerCount(int coverage, 
											  int upperCutoff)
{
	const int MIN_CUTOFF = 4;

	int totalReadLen = 0;
	for (auto& seqPair : _seqContainer.getIndex()) 
	{
		totalReadLen += _seqContainer.seqLen(seqPair.first);
	}
	int kmersNeeded = totalReadLen / coverage;
	LOG_PRINT("Genome size estimate: " << kmersNeeded / 2);

	int takenKmers = 0;
	int cutoff = 0;
	int repetitiveKmers = 0;
	int prevDiff = 0;
	for (auto mapPair = _vertexIndex.getKmerHist().rbegin();
		 mapPair != _vertexIndex.getKmerHist().rend(); ++mapPair)
	{
		if (mapPair->first <= upperCutoff)
		{
			takenKmers += mapPair->second;
			if (takenKmers >= kmersNeeded)
			{
				if (abs(takenKmers - kmersNeeded) < prevDiff)
				{
					cutoff = mapPair->first;
				}
				else
				{
					cutoff = mapPair->first + 1;
					takenKmers -= mapPair->second;
				}
				break;
			}
			prevDiff = abs(takenKmers - kmersNeeded);
		}
		else
		{
			repetitiveKmers += mapPair->second;
		}
	}

	if (cutoff < 4)
	{
		WARNING_PRINT("Unable to choose minimum kmer count cutoff. "
					  "Check if the coverage parameter is correct. "
					  "Running with default parameter t = " << MIN_CUTOFF);
		cutoff = 4;
	}
	
	LOG_PRINT("Filtered " << repetitiveKmers << " repetitive kmers");
	LOG_PRINT("Estimated minimum kmer coverage: " << cutoff <<
			  ", " << takenKmers << " unique kmers selected");

	return cutoff;
}
