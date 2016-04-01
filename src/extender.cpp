#include "extender.h"

std::vector<Extender::ReadPath> Extender::extendPaths()
{
	return std::vector<ReadPath>();
}


void Extender::extendRead(FastaRecord::ReadIdType readId)
{
	auto& overlaps = _ovlpDetector.getOverlapIndex().at(readId);
	for (auto& ovlp : overlaps)
	{
		//check if it a proper extension to the right
	}
}


void Extender::stepRight(FastaRecord::ReadIdType readId)
{
}
