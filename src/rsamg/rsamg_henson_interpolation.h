/*
 * Copyright (c) 2013:  G-CSC, Goethe University Frankfurt
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

#ifndef RSAMG_HENSON_INTERPOLATION_H_
#define RSAMG_HENSON_INTERPOLATION_H_

#include "rsamg_prolongation.h"

namespace ug{

/**
 *
 * @param P
 * @param A
 * @param i
 * @param coarse
 * @param strongf
 * @param diag
 */
template<typename TProlongationMatrix, typename TMatrix>
inline void SetHensonInterpolation(TProlongationMatrix &P,
		const TMatrix &A, size_t i,
		std::vector<typename TProlongationMatrix::connection > &coarse,
		std::vector<typename TProlongationMatrix::connection > &strongf,
		double diag)
{
	// sort coarse vector w.r.t index
	sort(coarse.begin(), coarse.end(), CompareByIndex<SparseMatrix<double>::connection>);

	typedef typename TMatrix::const_row_iterator const_row_iterator;

	// \sum_{m \in C_i} a_km for all k in N_i
	for(size_t k=0; k<strongf.size(); k++)
	{
		bool found;

		strongf[k].dValue = 0;
		const size_t kindex = strongf[k].iIndex;
		typename TMatrix::value_type sumAkm=0.0;
		for (size_t m=0; m<coarse.size(); m++)
		{
			const_row_iterator connkm = A.get_connection(kindex, coarse[m].iIndex, found);

			if (found) sumAkm += connkm.value();
		}
		strongf[k].dValue = BlockRef(A(i,kindex),0,0) / BlockRef(sumAkm, 0, 0);

		// find kj-connection
		for(size_t j=0; j<coarse.size(); j++)
		{
			const_row_iterator connkj = A.get_connection(kindex, coarse[j].iIndex, found);  // ATTENTION: deg^2*log deg

			if (found) {
				coarse[j].dValue += strongf[k].dValue * BlockRef(connkj.value(),0,0);
			}
		}
	}

	// correct weight by diag ( + sum of weak offdiag)
	for(size_t j=0; j<coarse.size(); j++)
	{coarse[j].dValue /= (-1.0*diag);}

	// set w_ij in matrix row P[i] forall coarse j.
	P.set_matrix_row(i, &coarse[0], coarse.size());
}



/**
 *
 * @param P
 * @param A
 * @param nodes
 * @param strong
 * @param theta
 * @param epsilonTruncation
 */
template<typename TProlongationMatrix, typename Matrix_type>
void CreateHensonProlongation(TProlongationMatrix &P, const Matrix_type &A, AMGNodes &nodes,
		//stdvector<int> &newIndex,
		const cgraph &strong, double theta, double epsilonTruncation=0.0)
{
	AMG_PROFILE_FUNC();
	P.resize(A.num_rows(), A.num_rows());

	std::vector<typename TProlongationMatrix::connection> coarse; coarse.reserve(255);
	std::vector<typename TProlongationMatrix::connection> strongf; coarse.reserve(255);

	// DIRECT INTERPOLATION
	for(size_t i=0; i < A.num_rows(); i++)
	{
		if(!nodes.is_master_or_inner(i)) continue;

		if(nodes[i].is_coarse())
		{
			// a coarse node
			P(i,i) = 1.0;
		}
		else if(A.is_isolated(i))
		{
			//P[i].initWithoutDiag(); // boundary values need not to be prolongated
		}
		else if(nodes[i].is_fine_direct())
		{
			// a non-interpolated fine node. calculate interpolation weights
			coarse.clear();
			strongf.clear();


			// step 1: set w'_ij = a_ij/a_jj for suitable j

			double diag = 0.0;
			for(typename Matrix_type::const_row_iterator conn = A.begin_row(i); conn != A.end_row(i); ++conn)
			{
				const size_t j = conn.index();

#ifdef SKIPDIRICHLET
				if(A.is_isolated(j)) continue;
#endif
				if( (j==i) || (!strong.has_connection(i,j)) ) {
					// skip diagonal and lump weak neighbors to diagonal
					diag += conn.value();
					continue;
				}

				typename TProlongationMatrix::connection c;
				c.iIndex = j;
				c.dValue = amg_offdiag_value(conn.value());

				// only strong neigbors here
				if(!nodes[j].is_coarse()
#ifdef USE_DIRICHLET_AS_INTERPOLATION_NODES
						&& !nodes[j].is_dirichlet()
#endif
				)
				{
					if(!nodes[j].is_dirichlet()) strongf.push_back(c);
					continue;
				}

				// add coarse nodes to interpolatory set
				if(!nodes[j].is_dirichlet()) coarse.push_back(c);
				UG_ASSERT(j >= 0, "not coarse?");
			}


			if(coarse.size() > 0)
				SetHensonInterpolation(P, A, i, coarse, strongf, diag);
			else
			{
				IF_DEBUG(DebugRSAMGProlongation, 2)
				{
				 /*UG_LOG("node " << i << " was marked direct, but no interpolation found. barrier = " << barrier << std::endl);
					UG_LOG("[\n");
					for(typename Matrix_type::const_row_iterator conn = A.begin_row(i); conn != A.end_row(i); ++conn)
					{
						if(conn.index() == i) continue; // skip diagonal
						if(!nodes[conn.index()].is_coarse())
						{	UG_LOG(" " << conn.index() << " is not coarse\n"); }
						else
						{	UG_LOG(" " << conn.index() << " has value " << amg_offdiag_value(conn.value()) << "\n"); }
					}
					UG_LOG("]\n");*/
				}
				//UG_ASSERT(0,i);
				// no suitable interpolating nodes for node i,
				// so this node has to be treated by indirect interpolation
				//nodes.set_unassigned_fine_indirect(i);
				nodes.set_coarse(i);
			}
		}
		/*else
		{
			UG_ASSERT(0, "?");
			//unassigned++;
			//UG_ASSERT(aggressiveCoarsening != 0, "no aggressive Coarsening but node " << i << " is fine and indirect??");
		}*/
	}

	if(nodes.get_unassigned())
		UG_DLOG(DebugRSAMGProlongation, 1, nodes.get_unassigned() << " ? ")
	if(nodes.get_unassigned_indirect_fine())
		UG_DLOG(DebugRSAMGProlongation, 1, "Pass 1: " << nodes.get_unassigned_indirect_fine() << " nodes need indirect interpolation.")
}

} // namespace ug
#endif /* RSAMG_HENSON_INTERPOLATION_H_ */

