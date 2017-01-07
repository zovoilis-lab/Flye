//(c) 2016 by Authors
//This file is a part of ABruijn program.
//Released under the BSD license (see LICENSE file)

#pragma once

#include "sequence_container.h"

class RepeatGraph
{
public:
	RepeatGraph(const SequenceContainer& asmSeqs):
		_asmSeqs(asmSeqs)
	{}

	void build();

private:
	const SequenceContainer& _asmSeqs;
};
