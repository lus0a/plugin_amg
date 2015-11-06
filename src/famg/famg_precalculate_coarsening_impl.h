/*
 * Copyright (c) 2012-2015:  G-CSC, Goethe University Frankfurt
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


#ifndef __H__LIB_ALGEBRA__AMG__FAMG_PRECALCULATE_COARSENING_IMPL_H__
#define __H__LIB_ALGEBRA__AMG__FAMG_PRECALCULATE_COARSENING_IMPL_H__

#include "common/util/histogramm.h"

namespace ug
{

/// \addtogroup plugins_amg_famg
/// \{

template<typename algebra_type>
void FAMGLevelCalculator<algebra_type>::calculate_all_possible_parent_pairs()
{
	AMG_PROFILE_FUNC();
	try{

	Stopwatch SW;;
	if(m_bTiming) { UG_LOG("\ncreating possible parent list... "); SW.start(); }



	possible_parents.clear();
	possible_parents.resize(A.num_rows());
	prolongation_calculated.resize(A.num_rows(), false);

	UG_LOG(std::scientific);

	for(size_t i=0; i<A.num_rows(); i++)
		if(A.is_isolated(i))
		{
			rating.set_dirichlet(i);
			continue;
		}

	PROGRESS_START(prog, A.num_rows(), "FAMG::Calculate all prossible parent pairs, " << A.num_rows() << " nodes, Level " << level);
	for(size_t i=0; i<A.num_rows(); i++)
	{
		if(rating.i_must_assign(i) == false) continue;
		if(rating[i].is_valid_rating() == false) continue; // could be set coarse from others
		double F = calculator.get_possible_parents(i, possible_parents[i], rating);
		if(i < fvalues.size())
			fvalues[i] = F;
		prolongation_calculated[i] = true;

		PROGRESS_UPDATE(prog, i);
	}


	//UG_LOG(HistogrammString(fvalues));
	//UG_LOG(DistributionPercentage(fvalues));
	if(m_bTiming) UG_LOG("took " << SW.ms() << " ms");

	}UG_CATCH_THROW_FUNC();
}


template<typename algebra_type>
void FAMGLevelCalculator<algebra_type>::precalculate_coarsening()
{
	AMG_PROFILE_FUNC();
	try{

	Stopwatch SW;
	if(m_bTiming) { UG_LOG("\ncoarsening... "); SW.start(); }

	size_t N = rating.size();
	for(size_t j=0; j<N; j++)
	{
		if(rating[j].is_valid_rating() && rating.i_must_assign(j) && A_OL2.is_isolated(j))
			rating.set_fine(j);
	}

	PROGRESS_START(prog, N, "FAMG::precalculate coarsening");

	std::vector<bool> mustAssign(N, false);
	while(heap.height() != 0)
	{
		PROGRESS_UPDATE(prog, N-heap.height());
		// get node i with best rating
		size_t i = heap.remove_max();
		neighT &n = possible_parents[i][0];


		IF_DEBUG(DBG_AMG_COARSENING,4)
		{
			UG_LOG("\n\n\nSelect next node...\n");
			UG_LOG("node " << i << " has rating " << rating[i] << ". now gets fine. parents: ");
			for(size_t j=0; j < n.parents.size(); j++)	UG_LOG(n.parents[j].from << " ");
			UG_LOG("\nUpdate Neighbors of " << i << "\n");
		}

		// node i then gets fine, parent nodes get coarse, and neighbors of those updated.

		// node i gets fine. update neighbors.
		rating.set_fine(i);
		UpdateNeighbors(SymmNeighGraph, i, possible_parents, rating, heap);

//		if(fvalues.size() > i) fvalues[i] = n.F;

		UG_DLOG(DBG_AMG_COARSENING, 4, "Set coarse parents:\n");
		// get parent pair, set as coarse (if not already done), update neighbors.

		for(size_t j=0; j < n.parents.size(); j++)
		{
			size_t node = n.parents[j].from;

			if(rating[node].is_coarse()) { UG_DLOG(DBG_AMG_COARSENING, 4, "\nnode [" << node << "] is already coarse\n"); }
			else { UG_DLOG(DBG_AMG_COARSENING, 4, "\nnode [" << node << "] has rating " << rating[node] << ". now gets coarse.\nUpdate Neighbors of [" << node << "]\n"); }

			if(!rating[node].is_coarse())
			{
				heap.remove(node);
				if(!rating[node].is_dirichlet())
				{
					rating.set_coarse(node);
					UpdateNeighbors(SymmNeighGraph, node, possible_parents, rating, heap);
				}
			}
			if(!rating[node].is_dirichlet())
				mustAssign[i] = true;
				//PoldIndices(i, node) = n.parents[j].value;
		}
	}

	for(size_t i=0;i<N; i++)
	{
		if(mustAssign[i]==false) continue;

//		rating.update_rating(i, possible_parents[i], true);

		neighT &n = possible_parents[i][0];
		if(fvalues.size() > i) fvalues[i] = n.F;
		for(size_t j=0; j < n.parents.size(); j++)
		{
			size_t node = n.parents[j].from;
			PoldIndices(i, node) = n.parents[j].value;
		}
	}

	if(m_bTiming) UG_LOG("took " << SW.ms() << " ms.");

	}UG_CATCH_THROW_FUNC();
}

/*
void get_best_parent_pair(size_t &i, neighborstruct2 &n)
{
	i = heap.remove_max();
	n = possible_parents[i][0];
}

void get_best_parent_pair_seqpar(size_t &i, neighborstruct2 &n)
{
	i = heap.get_max();
	neighborstruct2 &n0 = possible_parents[i][0];

	AlgebraID id = pcl::argmin(n.F, PN.local_to_global(i));
	neighborstruct2<AlgebraID> par;
	if(id.master_proc() == pcl::ProcRank())
	{
		for(size_t j=0; j < n0.parents.size(); j++)
		{
			par.parents[j].from = PN.local_to_global( n0.parents[j].from );
			par.parents[j].value= n.parents[j].value;
		}
		n = n0;
		heap.remove_max();
	}
	pcl::broadcast(par);

	if(id.master_proc() != pcl::ProcRank())
	{
		bool b;
		i = PN.get_local_index_if_available(id, b );
		for(size_t j=0; j < n0.parents.size(); j++)
		{
			n.parents[j].from = PN.get_local_index_if_available( par.parents[j].from, b);
			n.parents[j].value= par.parents[j].value;
		}
	}
}

template<typename algebra_type>
void FAMGLevelCalculator<algebra_type>::seqpar_precalculate_coarsening()
{
	AMG_PROFILE_FUNC();
	Stopwatch SW;
	UG_DLOG(LIB_ALG_AMG, 1, std::endl << "seqpar coarsening... "); if(bTiming) SW.start();

	size_t N = rating.size();
	for(size_t j=0; j<N; j++)
	{
		if(rating[j].is_valid_rating() && rating.i_must_assign(j) && A_OL2.is_isolated(j))
			rating.set_fine(j);
	}

	while(heap.height() != 0)
	{
		// get node i with best rating
		size_t i = heap.remove_max();
		neighborstruct2 &n = possible_parents[i][0];


		UG_DLOG(LIB_ALG_AMG, 4, "\n\n\nSelect next node...\n");
		UG_DLOG(LIB_ALG_AMG, 4, "node " << i << " has rating " << rating[i] << ". now gets fine. parents: ");
		for(size_t j=0; j < n.parents.size(); j++)	UG_DLOG(LIB_ALG_AMG, 2, n.parents[j].from << " ");
		UG_DLOG(LIB_ALG_AMG, 4, "\nUpdate Neighbors of " << i << "\n");

		// node i then gets fine, parent nodes get coarse, and neighbors of those updated.

		// node i gets fine. update neighbors.
		if(i != -1)
		{
			rating.set_fine(i);
			UpdateNeighbors(SymmNeighGraph, i, possible_parents, rating, heap);

			if(i < fvalues.size()) fvalues[i] = n.F;
		}



		UG_DLOG(LIB_ALG_AMG, 4, "Set coarse parents:\n");
		// get parent pair, set as coarse (if not already done), update neighbors.

		for(size_t j=0; j < n.parents.size(); j++)
		{
			size_t node = n.parents[j].from;
			if(node == -1) continue;

			if(rating[node].is_coarse()) { UG_DLOG(LIB_ALG_AMG, 4, "\nnode " << node << " is already coarse\n"); }
			else { UG_DLOG(LIB_ALG_AMG, 4, "\nnode " << node << " has rating " << rating[node] << ". now gets coarse.\nUpdate Neighbors of " << node << "\n"); }

			if(!rating[node].is_coarse())
			{
				heap.remove(node);
				if(!rating[node].is_dirichlet())
				{
					rating.set_coarse(node);
					UpdateNeighbors(SymmNeighGraph, node, possible_parents, rating, heap);
				}
			}
			if(!rating[node].is_dirichlet() && i != -1)
				PoldIndices(i, node) = n.parents[j].value;
		}
	}
	if(bTiming) UG_DLOG(LIB_ALG_AMG, 1, "took " << SW.ms() << " ms.");
}
*/

// end group plugins_amg_famg
/// \}

}

#endif
