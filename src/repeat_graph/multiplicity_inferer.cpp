//(c) 2016 by Authors
//This file is a part of ABruijn program.
//Released under the BSD license (see LICENSE file)

#include "multiplicity_inferer.h"
#include "graph_processing.h"
#include "../common/disjoint_set.h"
#include "../common/utils.h"
#include <cmath>


//Estimates the mean coverage and assingns edges multiplicity accordingly
void MultiplicityInferer::estimateCoverage()
{
	const int WINDOW = Config::get("coverage_estimate_window");
	const int SHORT_EDGE = Config::get("unique_edge_length");

	//alternative coverage
	std::unordered_map<GraphEdge*, std::vector<int32_t>> wndCoverage;

	for (auto& edge : _graph.iterEdges())
	{
		size_t numWindows = edge->length() / WINDOW;
		wndCoverage[edge].assign(numWindows, 0);
	}

	for (auto& path : _aligner.getAlignments())
	{
		for (size_t i = 0; i < path.size(); ++i)
		{
			auto& ovlp = path[i].overlap;
			auto& coverage = wndCoverage[path[i].edge];
			for (int pos = ovlp.extBegin / WINDOW + 1; 
			 	 pos < ovlp.extEnd / WINDOW; ++pos)
			{
				if (pos >= 0 && 
					pos < (int)coverage.size())
				{
					++coverage[pos];
				}
			}
		}
	}

	int64_t sumCov = 0;
	int64_t sumLength = 0;
	for (auto& edgeCoverage : wndCoverage)
	{
		if (edgeCoverage.first->length() < SHORT_EDGE) continue;
		for (auto& cov : edgeCoverage.second)
		{
			sumCov += (int64_t)cov;
			++sumLength;
		}
	}
	_meanCoverage = (sumLength != 0) ? sumCov / sumLength : 1;

	Logger::get().info() << "Mean edge coverage: " << _meanCoverage;

	std::vector<int32_t> edgesCoverage;
	for (auto edge : _graph.iterEdges())
	{
		if (wndCoverage[edge].empty()) continue;

		GraphEdge* complEdge = _graph.complementEdge(edge);
		int32_t medianCov = (median(wndCoverage[edge]) + 
						 	 median(wndCoverage[complEdge])) / 2;

		float minMult = (!edge->isTip()) ? 1 : 0;
		int estMult = std::max(minMult, std::round((float)medianCov / 
													_meanCoverage));
		if (estMult == 1)
		{
			edgesCoverage.push_back(medianCov);
		}

		//std::string match = estMult != edge->multiplicity ? "*" : " ";
		std::string covStr;

		Logger::get().debug() << edge->edgeId.signedId() << "\tlen:"
				<< edge->length() << "\tcov:" << medianCov << "\tmult:"
				<< (float)medianCov / _meanCoverage;

		//edge->multiplicity = estMult;
		edge->meanCoverage = medianCov;
	}

	_uniqueCovThreshold = !edgesCoverage.empty() ? q75(edgesCoverage) : 1;
	Logger::get().debug() << "Unique coverage threshold " << _uniqueCovThreshold;
}

//removes edges with low coverage support from the graph
void MultiplicityInferer::removeUnsupportedEdges()
{
	GraphProcessor proc(_graph, _asmSeqs, _readSeqs);
	auto unbranchingPaths = proc.getUnbranchingPaths();

	int32_t coverageThreshold = this->getMeanCoverage() / 
								Config::get("graph_cov_drop_rate");
	Logger::get().debug() << "Read coverage cutoff: " << coverageThreshold;

	std::unordered_set<GraphEdge*> edgesRemove;
	for (auto& path : unbranchingPaths)
	{
		if (!path.id.strand()) continue;

		if (path.meanCoverage <= coverageThreshold)
		{
			for (auto& edge : path.path)
			{
				edgesRemove.insert(edge);
				edgesRemove.insert(_graph.complementEdge(edge));
			}
		}
	}
	for (auto& edge : edgesRemove) _graph.removeEdge(edge);
	Logger::get().debug() << "Removed " << edgesRemove.size() / 2
		<< " unsupported edges";

	_aligner.updateAlignments();
}

void MultiplicityInferer::removeUnsupportedConnections()
{
	std::unordered_map<GraphEdge*, int32_t> rightConnections;
	std::unordered_map<GraphEdge*, int32_t> leftConnections;

	for (auto& readPath : _aligner.getAlignments())
	{
		if (readPath.size() < 2) continue;
		//int overhang = std::max(readPath.front().overlap.curBegin,
		//						readPath.back().overlap.curLen - 
		//							readPath.back().overlap.curEnd);
		//if (overhang > (int)Config::get("maximum_overhang")) continue;

		for (size_t i = 0; i < readPath.size() - 1; ++i)
		{
			if (readPath[i].edge == readPath[i + 1].edge &&
				readPath[i].edge->isLooped()) continue;
			if (readPath[i].edge->edgeId == 
				readPath[i + 1].edge->edgeId.rc()) continue;

			++rightConnections[readPath[i].edge];
			++leftConnections[readPath[i + 1].edge];
			GraphEdge* complLeft = _graph.complementEdge(readPath[i].edge);
			GraphEdge* complRight = _graph.complementEdge(readPath[i + 1].edge);
			++rightConnections[complRight];
			++leftConnections[complLeft];
		}
	}

	auto disconnectRight = [this](GraphEdge* edge)
	{
		GraphNode* newNode = _graph.addNode();
		vecRemove(edge->nodeRight->inEdges, edge);
		edge->nodeRight = newNode;
		edge->nodeRight->inEdges.push_back(edge);
	};
	auto disconnectLeft = [this](GraphEdge* edge)
	{
		GraphNode* newNode = _graph.addNode();
		vecRemove(edge->nodeLeft->outEdges, edge);
		edge->nodeLeft = newNode;
		edge->nodeLeft->outEdges.push_back(edge);
	};

	for (auto& edge : _graph.iterEdges())
	{
		if (!edge->edgeId.strand() || edge->isLooped()) continue;
		GraphEdge* complEdge = _graph.complementEdge(edge);

		int32_t coverageThreshold = edge->meanCoverage / 
								Config::get("graph_cov_drop_rate");

		//Logger::get().debug() << "Adjacencies: " << edge->edgeId.signedId() << " "
		//	<< leftConnections[edge] / 2 << " " << rightConnections[edge] / 2;

		if (!edge->nodeRight->isEnd() &&
			rightConnections[edge] / 2 < coverageThreshold)
		{
			Logger::get().debug() << "Chimeric right: " <<
				edge->edgeId.signedId() << " " << rightConnections[edge] / 2;

			disconnectRight(edge);
			disconnectLeft(complEdge);

			if (edge->selfComplement) continue;	//already discinnected
		}
		if (!edge->nodeLeft->isEnd() &&
			leftConnections[edge] / 2 < coverageThreshold)
		{
			Logger::get().debug() << "Chimeric left: " <<
				edge->edgeId.signedId() << " " << leftConnections[edge] / 2;

			disconnectLeft(edge);
			disconnectRight(complEdge);
		}
	}

	//GraphProcessor proc(_graph, _asmSeqs, _readSeqs);
	//proc.trimTips();
	_aligner.updateAlignments();
}

void MultiplicityInferer::separateHaplotypes()
{
	const float MAX_VARIATION = 0.25;

	GraphProcessor proc(_graph, _asmSeqs, _readSeqs);
	auto unbranchingPaths = proc.getUnbranchingPaths();

	std::unordered_set<FastaRecord::Id> toSeparate;
	for (auto& path : unbranchingPaths)
	{
		if (path.isLooped()) continue;

		std::vector<UnbranchingPath*> parallelPaths;
		for (auto& candEdge : unbranchingPaths)
		{
			if (candEdge.nodeLeft() == path.nodeLeft() &&
				candEdge.nodeRight() == path.nodeRight()) 
			{
				parallelPaths.push_back(&candEdge);
			}
		}

		if (parallelPaths.size() != 2) continue;
		if (parallelPaths[0]->id == parallelPaths[1]->id.rc()) continue;
		if (toSeparate.count(parallelPaths[0]->id) || 
			toSeparate.count(parallelPaths[1]->id)) continue;
		if (parallelPaths[0]->nodeLeft()->inEdges.size() != 1 ||
			parallelPaths[0]->nodeRight()->outEdges.size() != 1) continue;

		GraphEdge* entranceEdge = parallelPaths[0]->nodeLeft()->inEdges.front();
		GraphEdge* exitEdge = parallelPaths[0]->nodeRight()->outEdges.front();

		float covSum = parallelPaths[0]->meanCoverage + 
					   parallelPaths[1]->meanCoverage;
		
		float entranceDiff = fabsf(covSum - entranceEdge->meanCoverage) / covSum;
		float exitDiff = fabsf(covSum - exitEdge->meanCoverage) / covSum;
		if (entranceDiff > MAX_VARIATION || exitDiff > MAX_VARIATION) continue;

		//Logger::get().debug() << "Var: " << entranceDiff << " " << exitDiff;

		if (parallelPaths[0]->meanCoverage < parallelPaths[1]->meanCoverage)
		{
			toSeparate.insert(parallelPaths[0]->id);
			toSeparate.insert(parallelPaths[0]->id.rc());
		}
		else
		{
			toSeparate.insert(parallelPaths[1]->id);
			toSeparate.insert(parallelPaths[1]->id.rc());
		}
		//if (parallelEdges[0]->length() > Constants::trustedEdgeLength || 
		//	parallelEdges[1]->length() > Constants::trustedEdgeLength) continue;
	}

	for (auto& path : unbranchingPaths)
	{
		if (toSeparate.count(path.id))
		{
			GraphNode* newLeft = _graph.addNode();
			GraphNode* newRight = _graph.addNode();
			vecRemove(path.nodeLeft()->outEdges, path.path.front());
			vecRemove(path.nodeRight()->inEdges, path.path.back());
			path.nodeLeft() = newLeft;
			path.nodeRight() = newRight;
			path.nodeLeft()->outEdges.push_back(path.path.front());
			path.nodeRight()->inEdges.push_back(path.path.back());
		}
	}
	Logger::get().debug() << "Separated " << toSeparate.size() / 2 << " haplotypes";
}
