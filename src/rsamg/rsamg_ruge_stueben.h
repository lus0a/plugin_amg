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

#ifndef RSAMG_RUGE_STUEBEN_H_
#define RSAMG_RUGE_STUEBEN_H_

#include <cmath>
#include "rsamg_prolongation.h"
#include "lib_algebra/algebra_common/sparse_vector.h"

namespace ug{


/**
 * @brief sets (adds) the R/S Interpolation for a set of connections
 * @param P						Prolongation Matrix to set
 * @param i						node number to interpolate
 * @param conInterpolatory		strong matrix connections A(i, .) to use
 * @param sumNeighbors			sum A(i, j) over all neighbors
 * @param diag					A(i, i)
 * @param alphaScale			scaling for alpha (1.0 oder -1.0 for pos/neg)
 *
 * This function sets
 * sumInterpolatory = sum_{j in conInt} ( A(i, j) )
 * and for all j in conInterpolatory
 * P(i, j) = A(i, j) * sumNeighbors/sumInterpolatory * alphaScale/diag
 * with A(i, j) the matrix induced by conInterpolatory is meant.
 */
template<typename TProlongationMatrix>
void AddRSInterpolation(TProlongationMatrix &P, size_t i,
		std::vector<typename TProlongationMatrix::connection > &conInterpolatory, double sumNeighbors,
		double diag, double alphaScale, double epsilonTruncation)
{
//	PRINT_VECTOR(conNegative, "SetRSInterpolation- (PRE) node " << i);
	double sumInterpolatory=0.0;
	for(size_t j=0; j<conInterpolatory.size(); j++)
	{
		sumInterpolatory+=amg_signed_value(conInterpolatory[j].value());
	}
	// If signed contributions cancel, try the norm sum of the classified
	// group. A group of zero connections still has no valid denominator.
	if(fabs(sumInterpolatory) < 1e-30)
	{
		for(size_t j=0; j<conInterpolatory.size(); j++)
			sumInterpolatory += alphaScale * BlockNorm(conInterpolatory[j].value());
	}
	UG_COND_THROW(sumInterpolatory == 0.0 || !std::isfinite(sumInterpolatory),
			"RSAMG: invalid interpolation denominator at node " << i
			<< ": " << sumInterpolatory);
	UG_COND_THROW(diag == 0.0 || !std::isfinite(diag),
			"RSAMG: invalid interpolation diagonal at node " << i << ": " << diag);
	const double alpha = (sumNeighbors / sumInterpolatory) / diag * alphaScale;
	UG_COND_THROW(!std::isfinite(alpha),
			"RSAMG: non-finite interpolation scale at node " << i);
	for(size_t j=0; j<conInterpolatory.size(); j++)
		conInterpolatory[j].value() *= alpha;

//	PRINT_VECTOR(conNegative, "AddInterpolation (" << alphaScale << ") (POST) node " << i);
	TruncateInterpolation(conInterpolatory, epsilonTruncation);
//	PRINT_VECTOR(conNegative, "AddInterpolation (" << alphaScale << ") (POST-TRUNCATE) node " << i);

	// set w_ij in matrix row P[i] forall j.
	for(size_t j=0; j<conInterpolatory.size(); j++)
	{
		if(conInterpolatory[j].value() != 0.0)
			P(i, conInterpolatory[j].index()) += conInterpolatory[j].value();
	}
}


/**
 * Set the R/S Interpolation
 * @param P					Prolongation Matrix to set
 * @param i					node to interpolate
 * @param conPositive		strong positive connections to be used
 * @param sumPosNeighbors	sum over all positive connections
 * @param conNegative		strong negative connections to be used
 * @param sumNegNeighbors	sum over all negative connections
 * @param diag
 * @param epsilonTruncation
 */
template<typename TProlongationMatrix>
inline void SetRSInterpolation(TProlongationMatrix &P, size_t i,
		std::vector<typename TProlongationMatrix::connection > &conPositive, double sumPosNeighbors,
		std::vector<typename TProlongationMatrix::connection > &conNegative, double sumNegNeighbors,
		double diag, double epsilonTruncation=0.0)
{
	if (conPositive.size() > 0)
		AddRSInterpolation(P, i, conPositive, sumPosNeighbors, diag, 1.0, epsilonTruncation);
	else
		// lump onto diagonal
		diag += sumPosNeighbors;
	if(conNegative.size() > 0)
		AddRSInterpolation(P, i, conNegative, sumNegNeighbors, diag, -1.0, epsilonTruncation);
}

/**
 *
 * @param Birow
 * @param sumPosNeighbors
 * @param sumNegNeighbors
 * @param conPositive
 * @param conNegative
 * @param nodes
 */
template<typename TMatrix, typename TRowType, typename TConnection>
void GetInterpolatoryConnections(size_t i, const TMatrix &A, const cgraph &strong,
		TRowType &Birow, double &sumPosNeighbors, double &sumNegNeighbors,
		std::vector<TConnection> &conPositive, std::vector<TConnection> &conNegative, AMGNodes &nodes)
{
	conPositive.clear();
	conNegative.clear();
	for(typename TRowType::const_iterator conn = Birow.begin(); conn != Birow.end(); ++conn)
	{
		const size_t j = conn.index() ;
		if(j == i) continue; // skip diagonal

#ifdef SKIPDIRICHLET
		if(A.is_isolated(j))
			continue;
#endif

		double connValue = amg_offdiag_value(conn.value());
		// Sparse rows retain structural zeros, including eliminated entries.
		// Such entries cannot provide an interpolatory connection.
		if(connValue == 0.0) continue;
		if(connValue > 0) 	sumPosNeighbors += connValue;
		else				sumNegNeighbors += connValue;


		if (!strong.has_connection(i,j))
			continue;

		// only strong neigbors here
		if(!nodes[j].is_coarse()
#ifdef USE_DIRICHLET_AS_INTERPOLATION_NODES
			&& !nodes[j].is_dirichlet()
#endif
			) continue;

		if(nodes[j].is_dirichlet()) continue;

		// only coarse nodes here
		TConnection c(j, connValue);

		UG_ASSERT(c.index() >= 0, "not coarse?");

		// add to interpolatory set and update sums

		// 	con.push_back(c);
		if(connValue > 0)	{conPositive.push_back(c);} //sumPosInterpolatory += connValue;
		else				{conNegative.push_back(c);} //sumNegInterpolatory += connValue;
	}
}


///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// CreateProlongation:
//-------------------------
/**
 * Calculates Prolongation P with Matrix_type A and coarse/fine markers in
 * nodes[i].is_/is_coarse by direct interpolation
 * \param 	P					Matrix P: here goes the calculated prolongation
 * \param	A					Matrix A: matrix for which to calculate prolongation on next level
 * \param	nodes				fine/coarse marks of the nodes.
 * \param	graph				strong connections graph
 * \param   epsilonTruncation 	if != 0.0, truncates all entries from interpolation with p_ij < epsilonTruncation*max_j p_ij
 */
template<typename TProlongationMatrix, typename TMatrixType>
void CreateRugeStuebenProlongation(TProlongationMatrix &P, const TMatrixType &A,
		AMGNodes &nodes, const cgraph &strong, double epsilonTruncation=0.0)
{
	AMG_PROFILE_FUNC();
	PROGRESS_START(prog, nodes.size(), "CreateRugeStuebenProlongation size " << nodes.size());
	//UG_ASSERT(newIndex.size() == nodes.size(), "");
	P.resize_and_clear(A.num_rows(), A.num_rows());

	std::vector<typename TProlongationMatrix::connection> conPositive, conNegative;
	conPositive.reserve(255);
	conNegative.reserve(255);

	// DIRECT INTERPOLATION

	for(size_t i=0; i < A.num_rows(); i++)
	{
		PROGRESS_UPDATE(prog, i);
		if(!nodes.is_master_or_inner(i))
			continue;
		if(nodes[i].is_coarse())
		{
			// a coarse node
			//UG_ASSERT(newIndex[i] != -1, "coarse node but no new index?");
			//P(i, newIndex[i]) = 1.0;
			P(i,i) = 1.0;
		}
		else if(A.is_isolated(i))
		{
			//P[i].initWithoutDiag(); // boundary values need not to be prolongated
		}
		else if(nodes[i].is_fine_direct())
		{
			// a non-interpolated fine node. calculate interpolation weights

			double diag = amg_diag_value(A(i,i));

			double sumPosNeighbors=0.0, sumNegNeighbors=0.0;

			// Match the graph construction: only scalar RS eliminates positive
			// connections. Systems use the original couplings by block norm.
			if(block_traits<typename TMatrixType::value_type>::static_num_rows == 1)
			{
				// Eliminate connections
				typedef SparseVector<typename TMatrixType::value_type> LocalRow;
				LocalRow BiRow(A.num_cols());
				//StrongOffDiagCheck<typename TMatrixType::value_type> rule(strong);
				PositiveOffDiagCheck<typename TMatrixType::value_type> positive;
				EliminateEntries(positive, A, i, BiRow, true);

				GetInterpolatoryConnections(i, A, strong, BiRow,
						sumPosNeighbors, sumNegNeighbors, conPositive, conNegative, nodes);
			}
			else
			{
				ConstMatrixRow<TMatrixType> AiRow(A, i);
				GetInterpolatoryConnections(i, A, strong, AiRow,
						sumPosNeighbors, sumNegNeighbors, conPositive, conNegative, nodes);
			}

			// step 1: set w'_ij = a_ij/a_jj for suitable j


			if((conPositive.size() + conNegative.size()) > 0)
				SetRSInterpolation(P, i,
						conPositive, sumPosNeighbors,
						conNegative, sumNegNeighbors,
						diag, epsilonTruncation);
			else
			{
				IF_DEBUG(DebugRSAMGProlongation, 2)
				{
					UG_LOG("node " << i << " was marked direct, but no interpolation found." << std::endl);// barrier = " << barrier <<
					UG_LOG("[\n");
					for(typename TMatrixType::const_row_iterator conn = A.begin_row(i); conn != A.end_row(i); ++conn)
					{
						if(conn.index() == i) continue; // skip diagonal
						if(!nodes[conn.index()].is_coarse())
						{	UG_LOG(" " << conn.index() << " is not coarse\n"); }
						else
						{	UG_LOG(" " << conn.index() << " has value " << amg_offdiag_value(conn.value()) << "\n"); }
					}
					UG_LOG("]\n");
				}
				// Retain a node with no usable coarse neighbor on the coarse
				// grid, including its identity injection into the fine grid.
				nodes.set_coarse(i);
				P(i,i) = 1.0;
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
#endif /* RSAMG_RUGE_STUEBEN_H_ */
