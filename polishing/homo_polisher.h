//(c) 2016 by Authors
//This file is a part of ABruijn program.
//Released under the BSD license (see LICENSE file)

#pragma once

#include "subs_matrix.h"
#include "bubble.h"

class HomoPolisher
{
public:
	HomoPolisher(const SubstitutionMatrix& subsMatrix,
				 const HopoMatrix& hopoMatrix):
		_subsMatrix(subsMatrix), _hopoMatrix(hopoMatrix)
	{}
	void polishBubble(Bubble& bubble);

private:
	size_t mostLikelyLen(HopoMatrix::State state, 
						 const std::vector<HopoMatrix::Observation>& obs);

	const SubstitutionMatrix& _subsMatrix;
	const HopoMatrix& _hopoMatrix;
};
