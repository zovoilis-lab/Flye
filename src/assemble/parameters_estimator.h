//(c) 2016 by Authors
//This file is a part of ABruijn program.
//Released under the BSD license (see LICENSE file)

#include "../sequence/vertex_index.h"
#include "../sequence/sequence_container.h"

class ParametersEstimator
{
public:
	ParametersEstimator(const SequenceContainer& seqContainer,
						const VertexIndex& vertexIndex, int coverage):
		_vertexIndex(vertexIndex), 
		_seqContainer(seqContainer),
		_coverage(coverage)
	{}

	void estimateMinKmerCount(int upperCutoff);
	int genomeSizeEstimate();
	int minKmerCount() {return _minKmerCount;}
private:

	const VertexIndex& _vertexIndex;
	const SequenceContainer& _seqContainer;
	const int _coverage;
	int _takenKmers;
	int _minKmerCount;
};
