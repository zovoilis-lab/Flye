//(c) 2016 by Authors
//This file is a part of ABruijn program.
//Released under the BSD license (see LICENSE file)

#pragma once

#include "sequence_container.h"
#include "vertex_index.h"
#include "overlap.h"




struct Edge;
struct Knot
{
	Knot(size_t id): knotId(id) {}

	typedef size_t Id;
	size_t knotId;
	std::vector<Edge*> inEdges;
	std::vector<Edge*> outEdges;
};

struct Edge
{
	Edge(FastaRecord::Id id, int32_t seqBegin, 
		 int32_t seqEnd, Knot::Id knotId):
		 seqId(id), seqBegin(seqBegin), seqEnd(seqEnd), knotId(knotId)
	{}

	FastaRecord::Id seqId;
	int32_t seqBegin;
	int32_t seqEnd;
	Knot::Id knotId;
};

class AssemblyGraph
{
public:
	AssemblyGraph()
	{}

	void construct(OverlapDetector& ovlp);
	void outputDot(const std::string& filename);

	typedef std::unordered_map<FastaRecord::Id, 
					   std::vector<OverlapRange>> OverlapIndex;
private:
	std::vector<Knot> _knots;
	std::unordered_map<FastaRecord::Id, std::vector<Edge>> _edges;
};


