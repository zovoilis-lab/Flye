//(c) 2016 by Authors
//This file is a part of ABruijn program.
//Released under the BSD license (see LICENSE file)

#pragma once

#include "repeat_graph.h"
#include "repeat_resolver.h"

struct Contig
{
	Contig(const GraphPath& path, FastaRecord::Id id = FastaRecord::ID_NONE,
		   bool circular = false):
		   	path(path), id(id), circular(circular) {}
	GraphPath path;
	FastaRecord::Id id;
	bool circular;
};

class GraphProcessor
{
public:
	GraphProcessor(RepeatGraph& graph, const SequenceContainer& asmSeqs,
				   const SequenceContainer& readSeqs):
		_graph(graph), _asmSeqs(asmSeqs), _readSeqs(readSeqs),
		_tipThreshold(Parameters::get().minimumOverlap) {}

	void condence();
	void trimTips();
	void generateContigs();
	void dumpRepeats(const std::vector<GraphAlignment>& readAlignments,
					 const std::string& outFile);

	void outputContigsGraph(const std::string& filename);
	void outputContigsFasta(const std::string& filename);

private:
	void unrollLoops();
	void condenceEdges();
	void updateEdgesMultiplicity();

	RepeatGraph& _graph;
	const SequenceContainer& _asmSeqs;
	const SequenceContainer& _readSeqs;
	const int _tipThreshold;

	std::vector<Contig> _contigs;
};
