//(c) 2016 by Authors
//This file is a part of ABruijn program.
//Released under the BSD license (see LICENSE file)

#pragma once

#include "repeat_graph.h"
#include "graph_processing.h"

class OutputGenerator
{
public:
	OutputGenerator(RepeatGraph& graph, const SequenceContainer& asmSeqs,
				    const SequenceContainer& readSeqs):
		_graph(graph), _asmSeqs(asmSeqs), _readSeqs(readSeqs) {}

	void generateContigs();
	void dumpRepeats(const std::vector<GraphAlignment>& readAlignments,
					 const std::string& outFile);

	void outputDot(bool contigs, const std::string& filename);
	void outputGfa(bool contigs, const std::string& filename);
	void outputFasta(bool contigs, const std::string& filename);

private:
	void generateContigSequences(std::vector<UnbranchingPath>& paths) const;
	void outputEdgesDot(const std::vector<UnbranchingPath>& paths,
						const std::string& filename);
	void outputEdgesGfa(const std::vector<UnbranchingPath>& paths,
						const std::string& filename);
	void outputEdgesFasta(const std::vector<UnbranchingPath>& paths,
						  const std::string& filename);
	std::vector<UnbranchingPath> edgesPaths() const;

	RepeatGraph& _graph;
	const SequenceContainer& _asmSeqs;
	const SequenceContainer& _readSeqs;

	std::vector<UnbranchingPath> _contigs;
};
