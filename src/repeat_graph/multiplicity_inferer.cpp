//(c) 2016 by Authors
//This file is a part of ABruijn program.
//Released under the BSD license (see LICENSE file)

#include "multiplicity_inferer.h"
#include "graph_processing.h"
#include "../common/disjoint_set.h"
#include "../common/utils.h"


//Estimates the mean coverage and assingns edges multiplicity accordingly
void MultiplicityInferer::estimateCoverage()
{
	const int WINDOW = Config::get("coverage_estimate_window");
	const int SHORT_EDGE = Config::get("unique_edge_length");

	//alternative coverage
	std::unordered_map<GraphEdge*, std::vector<int>> wndCoverage;

	for (auto& edge : _graph.iterEdges())
	{
		int numWindows = edge->length() / WINDOW;
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
			sumCov += cov;
			++sumLength;
		}
	}
	_meanCoverage = (sumLength != 0) ? sumCov / sumLength : 1;

	Logger::get().debug() << "Mean edge coverage: " << _meanCoverage;

	std::vector<int> edgesCoverage;
	for (auto edge : _graph.iterEdges())
	{
		if (wndCoverage[edge].empty()) continue;

		GraphEdge* complEdge = _graph.complementEdge(edge);
		int medianCov = (median(wndCoverage[edge]) + 
						 median(wndCoverage[complEdge])) / 2;

		float minMult = (!edge->isTip()) ? 1 : 0;
		int estMult = std::max(minMult, 
							   roundf((float)medianCov / _meanCoverage));
		if (estMult == 1)
		{
			edgesCoverage.push_back(medianCov);
		}

		//std::string match = estMult != edge->multiplicity ? "*" : " ";
		std::string covStr;

		Logger::get().debug() << edge->edgeId.signedId() << "\t"
				<< edge->length() << "\t" << medianCov << "\t"
				<< (float)medianCov / _meanCoverage;

		//edge->multiplicity = estMult;
		edge->meanCoverage = medianCov;
	}

	_uniqueCovThreshold = q75(edgesCoverage);
	Logger::get().debug() << "Unique coverage threshold " << _uniqueCovThreshold;
}

//removes edges with low coverage support from the graph
void MultiplicityInferer::removeUnsupportedEdges()
{
	GraphProcessor proc(_graph, _asmSeqs, _readSeqs);
	auto unbranchingPaths = proc.getUnbranchingPaths();

	int coverageThreshold = this->getMeanCoverage() / 
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
	std::unordered_map<GraphEdge*, int> rightConnections;
	std::unordered_map<GraphEdge*, int> leftConnections;

	for (auto& readPath : _aligner.getAlignments())
	{
		if (readPath.size() < 2) continue;
		int overhang = std::max(readPath.front().overlap.curBegin,
								readPath.back().overlap.curLen - 
									readPath.back().overlap.curEnd);
		if (overhang > (int)Config::get("maximum_overhang")) continue;

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

	int coverageThreshold = this->getMeanCoverage() / 
							Config::get("graph_cov_drop_rate");
	for (auto& edge : _graph.iterEdges())
	{
		if (!edge->edgeId.strand() || edge->isLooped()) continue;
		GraphEdge* complEdge = _graph.complementEdge(edge);

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

	_aligner.updateAlignments();
}

void MultiplicityInferer::separateHaplotypes()
{
	std::unordered_set<GraphEdge*> toSeparate;
	for (auto& edge : _graph.iterEdges())
	{
		if (edge->isLooped()) continue;

		std::vector<GraphEdge*> parallelEdges;
		for (auto& parEdge : edge->nodeLeft->outEdges)
		{
			if (parEdge->nodeRight == edge->nodeRight) 
			{
				parallelEdges.push_back(parEdge);
			}
		}

		if (parallelEdges.size() != 2) continue;
		if (parallelEdges[0]->edgeId == parallelEdges[1]->edgeId.rc()) continue;
		if (toSeparate.count(parallelEdges[0]) || 
			toSeparate.count(parallelEdges[1])) continue;
		float covSum = parallelEdges[0]->meanCoverage + 
					   parallelEdges[1]->meanCoverage;
		if (covSum / this->getMeanCoverage() > 1.25) continue;

		if (parallelEdges[0]->meanCoverage < parallelEdges[1]->meanCoverage)
		{
			toSeparate.insert(parallelEdges[0]);
			toSeparate.insert(_graph.complementEdge(parallelEdges[0]));
		}
		else
		{
			toSeparate.insert(parallelEdges[1]);
			toSeparate.insert(_graph.complementEdge(parallelEdges[1]));
		}
		//if (parallelEdges[0]->length() > Constants::trustedEdgeLength || 
		//	parallelEdges[1]->length() > Constants::trustedEdgeLength) continue;
	}

	for (auto& edge : toSeparate)
	{
		GraphNode* newLeft = _graph.addNode();
		GraphNode* newRight = _graph.addNode();
		vecRemove(edge->nodeLeft->outEdges, edge);
		vecRemove(edge->nodeRight->inEdges, edge);
		edge->nodeLeft = newLeft;
		edge->nodeRight = newRight;
		edge->nodeLeft->outEdges.push_back(edge);
		edge->nodeRight->inEdges.push_back(edge);
	}
	Logger::get().debug() << "Separated " << toSeparate.size() / 2 << " haplotypes";
}
