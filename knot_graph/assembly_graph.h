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

	typedef int32_t Id;
	Id knotId;
	std::vector<Edge*> inEdges;
	std::vector<Edge*> outEdges;
};

struct Edge
{
	Edge(FastaRecord::Id id, int32_t seqBegin, 
		 int32_t seqEnd, Knot::Id knotBegin, Knot::Id knotEnd):

		 seqId(id), seqBegin(seqBegin), seqEnd(seqEnd), 
		 knotBegin(knotBegin), knotEnd(knotEnd)
	{}

	FastaRecord::Id seqId;
	int32_t seqBegin;
	int32_t seqEnd;
	Knot::Id knotBegin;
	Knot::Id knotEnd;
};

struct Connection
{
	Edge* inEdge;
	OverlapRange inOverlap;
	Edge* outEdge;
	OverlapRange outOverlap;
};

class AssemblyGraph
{
public:
	AssemblyGraph(const SequenceContainer& seqAssembly, 
				  const SequenceContainer& seqReads):
		_seqAssembly(seqAssembly), _seqReads(seqReads)
	{
		_knots.push_back(Knot(SEQ_BEGIN));
		_knots.push_back(Knot(SEQ_END));
	}

	void construct(const OverlapContainer& ovlp);
	void untangle(const OverlapContainer& ovlp);
	void outputDot(const std::string& filename);
	void outputFasta(const std::string& filename);

	typedef std::unordered_map<FastaRecord::Id, 
					   std::vector<OverlapRange>> OverlapIndex;
private:
	std::vector<Connection> getConnections(const OverlapContainer& ovlp);

	const SequenceContainer& _seqAssembly;
	const SequenceContainer& _seqReads;

	std::vector<Knot> _knots;
	std::unordered_map<FastaRecord::Id, std::list<Edge>> _edges;
	
	const Knot::Id SEQ_BEGIN = 0;
	const Knot::Id SEQ_END = 1;
};


