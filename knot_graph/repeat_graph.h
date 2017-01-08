//(c) 2016 by Authors
//This file is a part of ABruijn program.
//Released under the BSD license (see LICENSE file)

#pragma once

#include "sequence_container.h"

struct GluePoint
{
	GluePoint(size_t id = 0, FastaRecord::Id seqId = FastaRecord::ID_NONE,
			  int32_t position = 0):
		pointId(id), seqId(seqId), position(position) {}

	size_t 	pointId;
	FastaRecord::Id seqId;
	int32_t	position;
};


class RepeatGraph
{
public:
	RepeatGraph(const SequenceContainer& asmSeqs):
		_asmSeqs(asmSeqs)
	{}

	void build();
	void outputDot(const std::string& filename);

private:
	const SequenceContainer& _asmSeqs;

	std::unordered_map<FastaRecord::Id, 
					   std::vector<GluePoint>> _gluePoints;
	std::unordered_map<FastaRecord::Id, 
					   std::vector<std::pair<int32_t, int32_t>>> 
						   							_repetitiveRegions;
};
