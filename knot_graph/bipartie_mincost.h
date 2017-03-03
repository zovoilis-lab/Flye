//(c) 2016 by Authors
//This file is a part of ABruijn program.
//Released under the BSD license (see LICENSE file)

#pragma once

#include <vector>

typedef std::vector<std::vector<double>> BipartieTable;

std::vector<int> bipartieMincost(const BipartieTable& cost);
