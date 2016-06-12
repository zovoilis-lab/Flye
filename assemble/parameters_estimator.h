//(c) 2016 by Authors
//This file is a part of ABruijn program.
//Released under the BSD license (see LICENSE file)

#include "vertex_index.h"
#include "sequence_container.h"

class ParametersEstimator
{
public:
	ParametersEstimator():
		_vertexIndex(VertexIndex::get()), 
		_seqContainer(SequenceContainer::get())
	{}

	int estimateMinKmerCount(int coverage, int upperCutoff);
private:

	const VertexIndex& _vertexIndex;
	const SequenceContainer& _seqContainer;
};
