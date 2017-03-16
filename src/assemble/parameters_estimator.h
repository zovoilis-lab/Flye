//(c) 2016 by Authors
//This file is a part of ABruijn program.
//Released under the BSD license (see LICENSE file)

#include "../sequence/vertex_index.h"
#include "../sequence/sequence_container.h"

class ParametersEstimator
{
public:
	ParametersEstimator(const SequenceContainer& seqContainer,
						const VertexIndex& vertexIndex):
		_vertexIndex(vertexIndex), 
		_seqContainer(seqContainer)
	{}

	int estimateMinKmerCount(int coverage, int upperCutoff);
private:

	const VertexIndex& _vertexIndex;
	const SequenceContainer& _seqContainer;
};
