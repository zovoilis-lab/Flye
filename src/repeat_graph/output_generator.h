//(c) 2016 by Authors
//This file is a part of ABruijn program.
//Released under the BSD license (see LICENSE file)

//This module generates multiple types of output given the graph:
//it output edges sequences in FASTA format, the graph structure
//as dot or gfa. Also, output information about the unresolved
//repeats for the subsequent alasysis


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
	void extendContigs(const std::vector<GraphAlignment>& readAln, 
					   const std::string& outFile);

	//if contigs parameter is false, unbranching paths are
	//output as separate edges (useful for debugging)
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
