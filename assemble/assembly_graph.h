//(c) 2016 by Authors
//This file is a part of ABruijn program.
//Released under the BSD license (see LICENSE file)

#pragma once

#include "sequence_container.h"
#include "vertex_index.h"
#include "overlap.h"


class AssemblyGraph
{
public:
	AssemblyGraph():
		_seqContainer(SequenceContainer::get()),
		_vertexIndex(VertexIndex::get())
	{}

	void construct();
	void saveOverlaps(const std::string& filename);

	typedef std::unordered_map<FastaRecord::Id, 
					   std::vector<OverlapRange>> OverlapIndex;
private:
	enum JumpRes {J_END, J_INCONS, J_CLOSE, J_FAR};

	std::vector<OverlapRange> getReadOverlaps(FastaRecord::Id readId) const;
	JumpRes jumpTest(int32_t currentPrev, int32_t currentNext,
				     int32_t extensionPrev, int32_t extensionNext) const;
	bool    overlapTest(const OverlapRange& ovlp, 
						int32_t curLen, int32_t extLen) const;

	const SequenceContainer& _seqContainer;
	const VertexIndex& _vertexIndex;
	OverlapIndex _overlapIndex;
};
