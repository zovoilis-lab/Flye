//(c) 2016 by Authors
//This file is a part of ABruijn program.
//Released under the BSD license (see LICENSE file)

#include "multiplicity_inferer.h"
#include "../common/disjoint_set.h"
#include "../common/utils.h"


//Estimates the mean coverage and assingns edges multiplicity accordingly
void MultiplicityInferer::
	estimateCoverage(const std::vector<GraphAlignment>& readAln)
{
	const int WINDOW = Constants::coverageEstimateWindow;
	const int SHORT_EDGE = Constants::trustedEdgeLength;

	//alternative coverage
	std::unordered_map<GraphEdge*, std::vector<int>> wndCoverage;

	for (auto& edge : _graph.iterEdges())
	{
		int numWindows = edge->length() / WINDOW;
		wndCoverage[edge].assign(numWindows, 0);
	}

	for (auto& path : readAln)
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

		std::string match = estMult != edge->multiplicity ? "*" : " ";
		std::string covStr;
		/*for (int cov : altCoverage[edge])
		{
			covStr += std::to_string(cov) + " ";
		}*/
		Logger::get().debug() << match << "\t" << edge->edgeId.signedId() << "\t"
				<< edge->length() << "\t"
				<< edge->multiplicity << "\t" << estMult << "\t"
				<< medianCov << "\t"
				<< (float)medianCov / _meanCoverage;
		//Logger::get().debug() << covStr;

		edge->multiplicity = estMult;
		edge->meanCoverage = medianCov;
	}

	_uniqueCovThreshold = q75(edgesCoverage);
	Logger::get().debug() << "Unique coverage threshold " << _uniqueCovThreshold;
}

//removes edges with low coverage support from the graph
void MultiplicityInferer::removeUnsupportedEdges()
{
	int coverageThreshold = this->getMeanCoverage() / Constants::readCovRate;
	Logger::get().debug() << "Read coverage cutoff: " << coverageThreshold;

	std::unordered_set<GraphEdge*> edgesRemove;
	for (auto& edge : _graph.iterEdges())
	{
		GraphEdge* complEdge = _graph.complementEdge(edge);
		if (edge->meanCoverage <= coverageThreshold)
		{
			edgesRemove.insert(edge);
			edgesRemove.insert(complEdge);
		}
	}
	for (auto& edge : edgesRemove) _graph.removeEdge(edge);
	Logger::get().debug() << "Removed " << edgesRemove.size() 
		<< " unsupported edges";
}
