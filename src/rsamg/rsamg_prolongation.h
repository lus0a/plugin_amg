/*
 * Copyright (c) 2012-2013:  G-CSC, Goethe University Frankfurt
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


#ifndef __H__UG__LIB_ALGEBRA__AMG_SOLVER__AMG_RS_PROLONGATION_H__
#define __H__UG__LIB_ALGEBRA__AMG_SOLVER__AMG_RS_PROLONGATION_H__

#include <vector>
#include "rsamg_nodeinfo.h"
#include "rsamg_blockvalues.h"
#include "common/progress.h"
#include "rsamg_eliminate.h"
#include "../amg_misc.h"
#include "lib_algebra/algebra_common/unsorted_sparse_vector.h"

//#define USE_DIRICHLET_AS_INTERPOLATION_NODES
#define SKIPDIRICHLET
namespace ug {

/// \addtogroup plugins_amg_rsamg
/// \{

template<typename TConnection>
bool CompareByValue(const TConnection &i, const TConnection &j)
{
	return dabs(i.dValue) > dabs(j.dValue);
}

template<typename TConnection>
bool CompareByIndex(const TConnection &i, const TConnection &j)
{
	return i.iIndex < j.iIndex;
}


/**
 * Truncates entries smaller than (eps x maximum value)
 * elements with weight[j] < eps*max(weight) are set to zero, all other elements are rescaled
 * so that sum ( weightsNew ) = sum( weightsOld )
 * @param weight weight vector. truncated elements are set to 0.0
 * @param eps
 * @return
 */
template<typename TConnection>
bool TruncateInterpolation(std::vector<TConnection> &weight, double eps)
{
	// find maximum
	double sum = 0.0, dmax = 0.0;
	for(size_t j=0; j<weight.size(); j++)
	{
		sum += weight[j].dValue;
		dmax = std::max(weight[j].dValue, dmax);
	}

	// determine nodes to be removed
	double threshold=eps*dmax;
	double remove = 0.0;
	for(size_t j=0; j<weight.size(); j++)
	{
		if (weight[j].dValue <threshold)
		{
			remove += weight[j].dValue;
			weight[j].dValue= 0.0;
		}
	}

	// remove nodes
	if (fabs(remove) < 1e-8) return false;
	// rescale so that sum ( weightsNew ) = sum( weightsOld )
	const double rescale = 1.0 + remove/(sum-remove);
	for(size_t j=0; j<weight.size(); j++)
		weight[j].dValue *= rescale;

	return true;
};


/**
 *
 * @param P
 * @param i
 * @param con
 * @param sumPosNeighbors
 * @param sumNegNeighbors
 * @param sumPosInterpolatory
 * @param sumNegInterpolatory
 * @param diag
 */
inline void SetRSInterpolation_Sparsify(SparseMatrix<double> &P, size_t i,
		std::vector<SparseMatrix<double>::connection > &con,
		double sumPosNeighbors, double sumNegNeighbors,
		double sumPosInterpolatory, double sumNegInterpolatory,
		double diag)
{
	if(con.size()>=1000)
	{
		sort(con.begin(), con.end(), CompareByValue<SparseMatrix<double>::connection>);
		con.resize(6);
		for(size_t j=6; j<con.size(); j++)
		{
			if( con[i].dValue > 0)
				sumPosInterpolatory-=con[i].dValue;
			else
				sumNegInterpolatory-=con[i].dValue;
		}
	}

	// calculate alpha_i
	double alphaPos = - (sumPosNeighbors / sumPosInterpolatory) / diag;
	double alphaNeg = - (sumNegNeighbors / sumNegInterpolatory) / diag;
	// set w_ij = alpha * w'_ij = alpha * a_ij/a_jj.

	for(size_t j=0; j<con.size(); j++)
	{
		if(con[j].dValue < 0)
			con[j].dValue *= alphaNeg;
		else
			con[j].dValue *= alphaPos;
	}

	// set w_ij in matrix row P[i] forall j.
	P.set_matrix_row(i, &con[0], con.size());
}

/**
 *
 * @param A
 * @param i
 * @param minConnValue
 * @param maxConnValue
 * @param diag
 */
template<typename TMatrix>
void GetNeighborValues(const TMatrix &A, size_t i, double &minConnValue, double &maxConnValue, double &diag)
{
	diag = 0; maxConnValue=0; minConnValue = 1e12;
	for(typename TMatrix::const_row_iterator conn = A.begin_row(i); conn != A.end_row(i); ++conn)
	{
		if(conn.index() == i)
		{
			diag = amg_diag_value(conn.value());
			continue; // skip diagonal
		}
#ifdef SKIPDIRICHLET
		if(A.is_isolated(conn.index()))
			continue;
#endif
		double connValue = amg_offdiag_value(conn.value());
		if(minConnValue > connValue)
			minConnValue = connValue;
		else if(maxConnValue < connValue)
			maxConnValue = connValue;
	}
}

/**
 *
 * @param A
 * @param nodes
 * @param i
 * @param minConnValue
 * @param maxConnValue
 * @param minCoarseConnValue
 * @param maxCoarseConnValue
 * @param sumNegNeighbors
 * @param sumPosNeighbors
 * @param diag
 */
template<typename TMatrix, typename TNodes>
void GetNeighborValues(const TMatrix &A, TNodes &nodes, size_t i, double &minConnValue, double &maxConnValue,
		double &minCoarseConnValue, double &maxCoarseConnValue, double &sumNegNeighbors, double &sumPosNeighbors, double &diag)
{
	PROFILE_FUNC();
	maxConnValue = 0; minConnValue=1e12;
	maxCoarseConnValue = 0; minCoarseConnValue = 1e12;
	sumNegNeighbors = 0; sumPosNeighbors=0;
	diag=0;
	for(typename TMatrix::const_row_iterator conn = A.begin_row(i); conn != A.end_row(i); ++conn)
	{
		if(conn.index() == i)
		{
			diag = amg_diag_value(conn.value());
			continue;
		}
		double connValue = amg_offdiag_value(conn.value());

		if(connValue > 0)	sumPosNeighbors += connValue;
		else				sumNegNeighbors += connValue;

		if(maxConnValue < connValue) maxConnValue = connValue;
		if(minConnValue > connValue) minConnValue = connValue;

		if(nodes[conn.index()].is_coarse() || nodes[conn.index()].is_dirichlet())
		{
			if(maxCoarseConnValue > connValue)	maxCoarseConnValue = connValue;
			if(minCoarseConnValue < connValue)	minCoarseConnValue = connValue;
		}
	}
}


///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// CreateIndirectProlongation:
//-------------------------
/**
 * Assume Prolongation of all normal fine nodes is already computed, it calculates the Interpolation of
 * fineIndirect nodes with Matrix_type A and Coarse/FineIndirect markers in nodes[i].is_coarse/is_fine_indirect
 *
 * [AMGKS99] 7.2.2 Multi-pass-interpolation
 * Using A1 or A2 aggressive coarsening, it can be shown that we need at most 4 passes.
 * \note CreateIndirectProlongation is starts with pass 2, since pass 1 is CreateRugeStuebenProlongation.
 *
 * \param 	P					Matrix P: here goes the calculated prolongation
 * \param	A					Matrix A: matrix for which to calculate prolongation on next level
 * \param	nodes
 * \param	theta
 */
template<typename TProlongationMatrix, typename Matrix_type>
void CreateIndirectProlongation(TProlongationMatrix &P, const Matrix_type &A,
		AMGNodes &nodes, cgraph &strong)
{
	AMG_PROFILE_FUNC();
	std::vector<typename TProlongationMatrix::connection> conPositive, conNegative;
	conPositive.reserve(255);
	conNegative.reserve(255);


	std::vector<int> nrOfPaths;
	nrOfPaths.reserve(255);

	typedef UnsortedSparseVector<typename TProlongationMatrix::value_type> RowType;
	typedef typename RowType::iterator RowIterator;
	RowType row(A.num_cols());

	//P.print();
	// INDIRECT INTERPOLATION

	size_t oldUnassigned = -1;
	int pass=2;
	while(nodes.get_unassigned_indirect_fine())
	{
#ifdef AMG_PRINT_INDIRECT
		UG_DLOG(LIB_ALG_AMG, 1, std::endl);
#endif
		UG_DLOG(LIB_ALG_AMG, 1, "Pass " << pass << ": ");
		for(size_t i=0; i<A.num_rows() && nodes.get_unassigned_indirect_fine() > 0; i++)
		{
			if(!nodes[i].is_unassigned_fine_indirect() || A.is_isolated(i))
				continue;

			// calculate min offdiag-entry
			double diag, maxConnValue, minConnValue;
			GetNeighborValues(A, i, minConnValue, maxConnValue, diag);

			double sumPosNeighbors=0.0, sumNegNeighbors=0.0;
			row.clear();

			//double barrier = theta*std::max(maxConnValue, dabs(minConnValue));

			//cout << "indirect interpolating node " << i << endl;

			// look at neighbors of node i, try to interpolate indirectly through them
			for(typename Matrix_type::const_row_iterator conn = A.begin_row(i); conn != A.end_row(i); ++conn)
			{
				size_t indexN = conn.index();
				if(indexN == i) continue; // skip diagonal

				double connValue = amg_offdiag_value(conn.value());
				if(connValue > 0) 	sumPosNeighbors += connValue;
				else				sumNegNeighbors += connValue;

				if (!strong.has_connection(i,indexN))
						continue;

				// we dont want fine nodes which were indirectly interpolated in THIS pass
				// (no gauss-seidel-style)
				if(nodes[indexN].is_fine_indirect_level(pass))
					continue;
				// all interpolate neighbors are now from pass (pass-1)


				UG_ASSERT(!nodes[indexN].is_coarse(), "Node " << i << " indirect, but neighbor " <<  indexN << " coarse?");

				// now we look at the interpolation values of this neighbor (in P !!!)
				// and set w'(i,indexNN) += A(i,indexN) * P(indexN,indexNN)
				for(typename TProlongationMatrix::row_iterator conn2 = P.begin_row(indexN); conn2 != P.end_row(indexN); ++conn2)
				{
					size_t indexNN = conn2.index();
					AddMult(row(indexNN), connValue, conn2.value());
				}
			}

			if(row.size() == 0)
				continue;

			conPositive.clear();
			conNegative.clear();
			for(RowIterator it = row.begin(); it != row.end(); ++it)
			{
				if(it->value() > 0)
					conPositive.push_back(*it);
				else
					conNegative.push_back(*it);
			}

			nodes.set_fine_indirect_level(i, pass);

#ifdef AMG_PRINT_INDIRECT
			cout << i << " ";
#endif
			//cout << endl;

			SetRSInterpolation(P, i,
					conPositive, sumPosNeighbors,
					conNegative, sumNegNeighbors,
					diag);
		}


		if(nodes.get_unassigned_indirect_fine() == oldUnassigned)
		{
			UG_LOG(std::endl << "unassigned nodes left: " << std::endl);
			for(size_t i=0; i<A.num_rows(); i++)
			{
				if(nodes[i].is_unassigned_fine_indirect())
				 UG_LOG(i << " ");
			}
			UG_LOG("\n");
		}
		UG_ASSERT(nodes.get_unassigned_indirect_fine() != oldUnassigned, "Pass " << pass <<
				": Indirect Interpolation hangs at " << nodes.get_unassigned_indirect_fine() <<
				" unassigned indirect fine nodes.");

#ifdef AMG_PRINT_INDIRECT
		UG_DLOG(LIB_ALG_AMG, 1, "calculated, ");
#endif
		UG_DLOG(LIB_ALG_AMG, 1, nodes.get_unassigned_indirect_fine() << " left. ");
		pass++;
		oldUnassigned = nodes.get_unassigned_indirect_fine();
		break;
	}

	UG_ASSERT(nodes.get_unassigned_indirect_fine() == 0, "number of unassigned nodes is still " <<
			nodes.get_unassigned_indirect_fine());

	P.defragment();
	//P.print();
}
// end group plugins_amg_rsamg
/// \}

} // namespace ug

#include "rsamg_ruge_stueben.h"
#include "rsamg_henson_interpolation.h"

#endif /* __H__UG__LIB_ALGEBRA__AMG_SOLVER__AMG_RS_PROLONGATION_H__ */
