/*
 * Copyright (c) 2013-2015:  G-CSC, Goethe University Frankfurt
 * Author: Martin Rupp
 * 
 * This file is part of UG4.
 * 
 * UG4 is free software: you can redistribute it and/or modify it under the
 * terms of the GNU Lesser General Public License version 3 (as published by the
 * Free Software Foundation) with the following additional attribution
 * requirements (according to LGPL/GPL v3 §7):
 * 
 * (1) The following notice must be displayed in the Appropriate Legal Notices
 * of covered and combined works: "Based on UG4 (www.ug4.org/license)".
 * 
 * (2) The following notice must be displayed at a prominent place in the
 * terminal output of covered works: "Based on UG4 (www.ug4.org/license)".
 * 
 * (3) The following bibliography is recommended for citation and must be
 * preserved in all covered files:
 * "Reiter, S., Vogel, A., Heppner, I., Rupp, M., and Wittum, G. A massively
 *   parallel geometric multigrid solver on hierarchically distributed grids.
 *   Computing and visualization in science 16, 4 (2013), 151-164"
 * "Vogel, A., Reiter, S., Rupp, M., Nägel, A., and Wittum, G. UG4 -- a novel
 *   flexible software system for simulating pde based models on high performance
 *   computers. Computing and visualization in science 16, 4 (2013), 165-179"
 * 
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 * GNU Lesser General Public License for more details.
 */

#ifndef NEIGHBORHOOD_H_
#define NEIGHBORHOOD_H_

namespace ug{

template<typename TMatrix>
void GetNeighborhood(TMatrix &A, size_t node, std::vector<size_t> &onlyN1)
{
	onlyN1.clear();
	for (typename TMatrix::const_row_iterator it = A.begin_row(node);
			it != A.end_row(node); ++it)
	{
		if (it.value() != 0.0 && it.index() != node)
			onlyN1.push_back(it.index());
	}
}

// adds the indices of the neighbors of a node to the vector onlyN1, if they have visited flag 0.
// node: neither onlyN1 nor bvisited get reset
template<typename TMatrix>
void GetNeighborhoodRec(TMatrix &A, size_t node,
		std::vector<size_t> &onlyN1, std::vector<bool> &bvisited)
{
	bvisited[node] = true;

	for (typename TMatrix::const_row_iterator it = A.begin_row(node);
			it != A.end_row(node); ++it)
	{
		if (bvisited[it.index()] == false && it.value() != 0.0)
		{
			bvisited[it.index()] = true;
			onlyN1.push_back(it.index());
		}
	}
}
template<typename TMatrix>
void GetNeighborhoodRec(TMatrix &A, size_t node,
		std::vector<size_t> &onlyN1, std::vector<size_t> &onlyN2,
		std::vector<bool> &bvisited)
{
	bvisited[node] = true;

	GetNeighborhoodRec(A, node, onlyN1, bvisited);
	for (size_t i = 0; i < onlyN1.size(); i++)
		GetNeighborhoodRec(A, onlyN1[i], onlyN2, bvisited);
}
inline void CleanupFlagArray(std::vector<bool> &bFlagged,
		const std::vector<size_t> &iFlagged)
{
	for (size_t i = 0; i < iFlagged.size(); i++)
		bFlagged[iFlagged[i]] = false;
}

template<typename TMatrix>
void GetNeighborhood(TMatrix &A, size_t node, std::vector<size_t> &onlyN1,
		std::vector<size_t> &onlyN2, std::vector<bool> &bvisited)
{
	onlyN1.clear();
	onlyN2.clear();
	GetNeighborhoodRec(A, node, onlyN1, onlyN2, bvisited);

	// cleanup
	bvisited[node] = false;
	CleanupFlagArray(bvisited, onlyN1);
	CleanupFlagArray(bvisited, onlyN2);
}

template<typename TMatrix>
void GetNeighborhood(TMatrix &A, size_t node, std::vector<size_t> &onlyN1,
		std::vector<size_t> &onlyN2, std::vector<size_t> &onlyN3, std::vector<bool> &bvisited)
{
	onlyN1.clear();
	onlyN2.clear();
	onlyN3.clear();

	bvisited[node] = true;

	GetNeighborhoodRec(A, node, onlyN1, bvisited);
	for (size_t i = 0; i < onlyN1.size(); i++)
		GetNeighborhoodRec(A, onlyN1[i], onlyN2, bvisited);
	for (size_t i = 0; i < onlyN2.size(); i++)
		GetNeighborhoodRec(A, onlyN2[i], onlyN3, bvisited);
	// cleanup
	bvisited[node] = false;
	CleanupFlagArray(bvisited, onlyN1);
	CleanupFlagArray(bvisited, onlyN2);
	CleanupFlagArray(bvisited, onlyN3);
}

}

#endif /* NEIGHBORHOOD_H_ */
