#include "haplotype_resolver.h"
#include "graph_processing.h"

//This function collapses simply bubbles caused by
//alternative haplotypes / strains. They are defined as follows:
//1. Structure: 1 input, 2 branches, 1 output: -<>-
//2. Size of each branch is shorter than MAX_BUBBLE_LEN below
//3. Total coverage of bubbles branches roughly equasl to input/output coverages
//4. Each branch is shorter than both entrace and exits. We need this to
//   distinguish from the case of two repeats of multiplicity 2
//Note that we are not using any global coverage assumptions here.
int HaplotypeResolver::collapseHeterozygousBulges(bool removeAlternatives)
{
	const float MAX_COV_VAR = 1.5;
	const int MAX_BUBBLE_LEN = Config::get("max_bubble_length");

	GraphProcessor proc(_graph, _asmSeqs);
	auto unbranchingPaths = proc.getUnbranchingPaths();

	std::unordered_set<FastaRecord::Id> toSeparate;
	int numMasked = 0;
	for (auto& path : unbranchingPaths)
	{
		if (path.isLooped()) continue;

		std::vector<UnbranchingPath*> twoPaths;
		for (auto& candEdge : unbranchingPaths)
		{
			if (candEdge.nodeLeft() == path.nodeLeft() &&
				candEdge.nodeRight() == path.nodeRight()) 
			{
				twoPaths.push_back(&candEdge);
			}
		}

		//making sure the structure is ok
		if (twoPaths.size() != 2) continue;
		if (twoPaths[0]->id == twoPaths[1]->id.rc()) continue;
		if (toSeparate.count(twoPaths[0]->id) || 
			toSeparate.count(twoPaths[1]->id)) continue;
		if (twoPaths[0]->nodeLeft()->inEdges.size() != 1 ||
			twoPaths[0]->nodeLeft()->outEdges.size() != 2 ||
			twoPaths[0]->nodeRight()->outEdges.size() != 1 ||
			twoPaths[0]->nodeRight()->inEdges.size() != 2) continue;

		UnbranchingPath* entrancePath = nullptr;
		UnbranchingPath* exitPath = nullptr;
		for (auto& cand : unbranchingPaths)
		{
			if (cand.nodeRight() == 
				twoPaths[0]->nodeLeft()) entrancePath = &cand;
			if (cand.nodeLeft() == twoPaths[0]->nodeRight()) exitPath = &cand;
		}

		//sanity check for maximum bubble size
		if (std::max(twoPaths[0]->length, twoPaths[1]->length) > 
			MAX_BUBBLE_LEN) continue;

		//coverage requirement: sum over two branches roughly equals to
		//exit and entrance coverage or less
		float covSum = twoPaths[0]->meanCoverage + twoPaths[1]->meanCoverage;
		if (covSum > std::min(entrancePath->meanCoverage * MAX_COV_VAR,
							  exitPath->meanCoverage * MAX_COV_VAR)) continue;

		//require bubble branches to be shorter than entrance or exit,
		//to distinguish from the case of two consecutive repeats
		//of multiplicity 2
		if (std::max(twoPaths[0]->length, twoPaths[1]->length) >
			std::max(entrancePath->length, exitPath->length)) continue;

		if (twoPaths[0]->meanCoverage > twoPaths[1]->meanCoverage)
		{
			std::swap(twoPaths[0], twoPaths[1]);
		}

		if (!twoPaths[0]->path.front()->altHaplotype ||
			!twoPaths[1]->path.front()->altHaplotype) ++numMasked;

		for (size_t i = 0; i < 2; ++i)
		{
			for (auto& edge : twoPaths[i]->path)
			{
				edge->altHaplotype = true;
				_graph.complementEdge(edge)->altHaplotype = true;
			}
		}

		if (removeAlternatives)
		{
			toSeparate.insert(twoPaths[0]->id);
			toSeparate.insert(twoPaths[0]->id.rc());
			for (auto& edge : twoPaths[1]->path)
			{
				edge->meanCoverage += twoPaths[0]->meanCoverage;
				_graph.complementEdge(edge)->meanCoverage += twoPaths[0]->meanCoverage;
				edge->altHaplotype = false;
				_graph.complementEdge(edge)->altHaplotype = false;
			}
		}
	}

	if (removeAlternatives)
	{
		for (auto& path : unbranchingPaths)
		{
			if (toSeparate.count(path.id))
			{
				//Logger::get().debug() << "Seperated branch: " << path.edgesStr();

				GraphNode* newLeft = _graph.addNode();
				GraphNode* newRight = _graph.addNode();
				vecRemove(path.nodeLeft()->outEdges, path.path.front());
				vecRemove(path.nodeRight()->inEdges, path.path.back());
				path.nodeLeft() = newLeft;
				path.nodeRight() = newRight;
				newLeft->outEdges.push_back(path.path.front());
				newRight->inEdges.push_back(path.path.back());
			}
		}

		Logger::get().debug() << "[SIMPL] Removed " << toSeparate.size() / 2 
			<< " heterozygous bulges";

		_aligner.updateAlignments();
		return toSeparate.size() / 2;
	}
	else
	{
		Logger::get().debug() << "[SIMPL] Masked " << numMasked
			<< " heterozygous bulges";
		return numMasked;
	}
}

//This function collapses simple loops:
//1. One loop edge with one entrance and one exit
//2. Loop length is shorter than lengths of entrance/exit
//3. Loop coverage is roughly equal or less than coverage of entrance/exit
int HaplotypeResolver::collapseHeterozygousLoops(bool removeAlternatives)
{
	const float COV_MULT = 1.5;

	GraphProcessor proc(_graph, _asmSeqs);
	auto unbranchingPaths = proc.getUnbranchingPaths();

	std::unordered_set<FastaRecord::Id> toUnroll;
	std::unordered_set<FastaRecord::Id> toRemove;
	int numMasked = 0;
	for (auto& loop : unbranchingPaths)
	{
		if (!loop.id.strand()) continue;
		if (!loop.isLooped()) continue;
		if (loop.path.front()->selfComplement) continue;

		GraphNode* node = loop.nodeLeft();
		if (node->inEdges.size() != 2 ||
			node->outEdges.size() != 2) continue;

		UnbranchingPath* entrancePath = nullptr;
		UnbranchingPath* exitPath = nullptr;
		for (auto& cand : unbranchingPaths)
		{
			if (cand.nodeRight() == node &&
				loop.id != cand.id) entrancePath = &cand;
			if (cand.nodeLeft() == node &&
				loop.id != cand.id) exitPath = &cand;
		}

		if (entrancePath->isLooped()) continue;
		if (entrancePath->id == exitPath->id.rc()) continue;

		//loop coverage should be roughly equal or less
		if (loop.meanCoverage > 
				COV_MULT * std::min(entrancePath->meanCoverage, 
									entrancePath->meanCoverage)) continue;

		//loop should not be longer than other branches
		if (loop.length > std::max(entrancePath->length, 
								   exitPath->length)) continue;

		if (!loop.path.front()->altHaplotype) ++numMasked;
		for (auto& edge : loop.path)
		{
			edge->altHaplotype = true;
			_graph.complementEdge(edge)->altHaplotype = true;
		}
		//either remove or unroll loop, depending on the coverage
		if (loop.meanCoverage < 
			(entrancePath->meanCoverage + exitPath->meanCoverage) / 4)
		{
			toRemove.insert(loop.id);
			toRemove.insert(loop.id.rc());
		}
		else
		{
			toUnroll.insert(loop.id);
			toUnroll.insert(loop.id.rc());
		}
	}

	if (removeAlternatives)
	{
		for (auto& path : unbranchingPaths)
		{
			if (toUnroll.count(path.id))
			{
				//Logger::get().debug() << "Unrolled loop: " << path.edgesStr();

				GraphNode* newNode = _graph.addNode();
				size_t id = path.nodeLeft()->inEdges[0] == path.path.back();
				GraphEdge* prevEdge = path.nodeLeft()->inEdges[id];

				vecRemove(path.nodeLeft()->outEdges, path.path.front());
				vecRemove(path.nodeLeft()->inEdges, prevEdge);
				path.nodeLeft() = newNode;
				newNode->outEdges.push_back(path.path.front());
				prevEdge->nodeRight = newNode;
				newNode->inEdges.push_back(prevEdge);
			}
			if (toRemove.count(path.id))
			{
				//Logger::get().debug() << "Removed loop: " << path.edgesStr();

				GraphNode* newLeft = _graph.addNode();
				GraphNode* newRight = _graph.addNode();

				vecRemove(path.nodeLeft()->outEdges, path.path.front());
				vecRemove(path.nodeLeft()->inEdges, path.path.back());
				path.nodeLeft() = newLeft;
				newRight->inEdges.push_back(path.path.back());
				path.nodeRight() = newRight;
				newLeft->outEdges.push_back(path.path.front());
			}
		}

		Logger::get().debug() << "[SIMPL] Removed " << (toRemove.size() + toUnroll.size()) / 2
			<< " heterozygous loops";
		_aligner.updateAlignments();
		return (toRemove.size() + toUnroll.size()) / 2;
	}
	else
	{
		Logger::get().debug() << "[SIMPL] Masked " << numMasked << " heterozygous loops";
		return numMasked;
	}

}
