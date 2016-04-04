#pragma once

#include <vector>

#include "fasta.h"
#include "extender.h"
#include "overlap.h"

class ContigGenerator
{
public:
	ContigGenerator(const Extender& extender, 
					const OverlapDetector& overlapDetector,
					const SequenceContainer& seqContainer):
		_extender(extender), _overlapDetector(overlapDetector),
		_seqContainer(seqContainer) {}
	
	void generateContigs();
	void outputContigs(const std::string& fileName);
	
private:
	const Extender& _extender;
	const OverlapDetector& _overlapDetector;
	const SequenceContainer& _seqContainer;

	std::tuple<int32_t, int32_t> getSwitchPositions(FastaRecord::ReadIdType leftRead, 
													FastaRecord::ReadIdType rightRead,
													int32_t prevSwitch);

	std::vector<std::vector<FastaRecord>> _contigs;
};
