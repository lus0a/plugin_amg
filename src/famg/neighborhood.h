/*
 * neighborhood.h
 *
 *  Created on: 20.09.2013
 *      Author: mrupp
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

}

#endif /* NEIGHBORHOOD_H_ */
