#include "contig.h"

void ContigGenerator::generateContigs()
{
	const auto& readPaths = _extender.getContigPaths();
	for (const Extender::ReadPath& path : readPaths)
	{
	}
}


std::tuple<int32_t, int32_t> 
ContigGenerator::getSwitchPositions(FastaRecord::ReadIdType leftRead, 
									FastaRecord::ReadIdType rightRead,
									int32_t prevSwitch)
{
	return std::make_tuple(0, 0);
}
