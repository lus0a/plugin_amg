/*
 * Copyright (c) 2009-2014:  G-CSC, Goethe University Frankfurt
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


#ifndef __H__UG__LIB_ALGEBRA__SEMIAMG_SOLVER__SEMIAMG_IMPL_H__
#define __H__UG__LIB_ALGEBRA__SEMIAMG_SOLVER__SEMIAMG_IMPL_H__

#include "boost/tuple/tuple.hpp"

#include "semi_amg.h"
#include "lib_algebra/common/graph/graph.h"



#include "common/stopwatch.h"

#include "../amg_debug_helper.h"
#include "../amg_misc.h"
#include "../amg_profiling.h"
#include "../amg_parallel.h"

#include "../rsamg/rsamg_nodeinfo.h"
#include "../rsamg/strong_connections.h"

#include "common/error.h"

#ifdef UG_PARALLEL
#include "pcl/pcl.h"
#include "lib_algebra/parallelization/parallelization_util.h"
#include "lib_algebra/parallelization/parallel_matrix_overlap_impl.h"
#include "../rsamg/rsamg_parallel_coarsening.h"
#endif



namespace ug{

/// \addtogroup plugins_amg_rsamg
/// \{


//!
//! amg constructor
template<typename TAlgebra>
SemiAMG<TAlgebra>::SemiAMG()
: AMGBase<TAlgebra>(), m_tli(), m_dTheta(0.25), m_bCreateVirtualNodes(false)
{
#ifdef UG_PARALLEL

#endif
}



template<typename algebra_type>
size_t SemiAMG<algebra_type>::user_defined_coarsening(const cgraph &strong, size_t level, AMGNodes &nodes)
{
	cAMG_helper &helper= m_amghelper;

	const size_t nnodes = helper.positions[level].size();
	const size_t gnodes = strong.size();

	THROW_IF_NOT_EQUAL(nnodes, gnodes);
	size_t nfine = 0;
	for(size_t i=0; i<nnodes; i++)
	{
		MathVector<3> &coord =  helper.positions[level][i];
		size_t flevel =m_tli.template get_node_level<3>(coord);
		if ( flevel <=level)
		{
			nodes.set_fine(i);
			nfine++;
		}
		else
		{
			nodes.external_set_coarse(i);
		}
	}
	UG_LOG("\nUser defined C/F-splitting on level " << level <<": " << nfine << " fine nodes out of "<< nnodes << "gnodes " << gnodes);

	return nnodes;
}

inline static int num_ijk_coord_matches(const MathVector<3> &coordi, const MathVector<3> &coordj)
{
	return ((int) (coordi[0]==coordj[0])) + ((int) (coordi[1]==coordj[1])) +((int) (coordi[2]==coordj[2]));
}

template<typename algebra_type>
void SemiAMG<algebra_type>::user_defined_prolongation(const cgraph &S, const matrix_type &A, prolongation_matrix_type &P, size_t level, AMGNodes &nodes)
{

	cAMG_helper &helper= m_amghelper;


	AMG_PROFILE_FUNC();
//	UG_SET_DEBUG_LEVEL(LIB_ALG_AMG, m_famg.iDebugLevelExternalCoarsening);

	Stopwatch SW;
	UG_DLOG(LIB_ALG_AMG, 0, std::endl << "Computing prolongation... ");
	if(bTiming) SW.start();

	const size_t N = A.num_rows();

	PROGRESS_START(prog, 2*N, "SemiAMG::user_defined_prolongation " << N);
	P.resize_and_clear(N, N);


	// init virtual coarse
	TVirtualCoarseMap & vmap = vcoarse[level];
	size_t & vindex = vcount[level];
	vcoarse[level] = TVirtualCoarseMap();
	vcount[level] = 0;


	typedef boost::tuple<int, int, double> pentry_type;
	std::vector<pentry_type> poffd;
	std::vector<cAMG_helper::position_type> vPositions;

	// first pass
	const stdvector<MathVector<3> > & lpositions = helper.positions[level];
	const size_t nnodes = lpositions.size();
	size_t nskipped = 0;
	for(size_t i=0; i<nnodes; i++)
	{
		PROGRESS_UPDATE(prog, i);

		// map coarse grid points by identity
		if (nodes[i].is_coarse())
		  {
			 //	std::cout << "Ignoring " << i << std::endl;
		     P(i, i) = 1.0; continue;
		  }

		// skip off-proc nodes
		if(!nodes.is_master_or_inner(i))
		{
			//std::cout << "Skipping " << i << std::endl;
			nskipped++; continue;
		}

		// skip isolated (Dirichlet) nodes
		if (S.is_isolated(i))
		{
				//std::cout << "Skipping " << i << std::endl;
				nskipped++; continue;
		}
		  
		// continue for fine nodes here
		{
		  const MathVector<3> &icoord =  lpositions[i];
		  MathVector<3> nbcoords[27];
		  double nbweights[27];
		  const int nparents = m_tli.get_parents(icoord, level, nbcoords, nbweights);
		  UG_ASSERT(nparents<=8, "Huhh: Too many neighbors???");
		 

		  // map every (coarse) nb to best matching (fine/coarse) equivalent on graph
		  MathVector<3> bestcoord[27];
		  size_t bestindex[27];
		  size_t bestmatch[27];

		  for (int nb=0; nb<nparents; ++nb)
		  {
			  // default match: node i
			  bestcoord[nb] = icoord;
			  bestindex[nb] = -1;
			  bestmatch[nb] = -1; // num_ijk_coord_matches(nbcoords[nb], icoord);

			  int maxjmatch = bestmatch[nb];

			  // search graph for (unique) best coordinate match
			  // we found a coarse grid node, if all coordinates match
			  for(typename cgraph::const_row_iterator iconn = S.begin_row(i); iconn != S.end_row(i); ++iconn)
			  {
				  // find coarse neighbors
				  const size_t j = *iconn;
				  if (!nodes[j].is_coarse()) continue;

				  // continue for coarse nodes
				  const MathVector<3> &jcoord =  lpositions[j];
				  const int jmatch = num_ijk_coord_matches(nbcoords[nb], jcoord);

				  if (jmatch > maxjmatch)
				  {
					  maxjmatch = jmatch;

					  bestcoord[nb] = jcoord;
					  bestindex[nb] = j;
					  bestmatch[nb] = jmatch;
				  }

			  }

			  if (bestmatch[nb] == 3) continue;

			  for(typename cgraph::const_row_iterator iconn = S.begin_row(i); iconn != S.end_row(i); ++iconn)
			  {
				  // find (indirect) neighbors
				  const size_t j = *iconn;

				  if (nodes[j].is_coarse()) continue;

				  for(typename cgraph::const_row_iterator jconn = S.begin_row(j); jconn != S.end_row(j); ++jconn)
				  {
					  const size_t k = *jconn;
					  if (!nodes[k].is_coarse()) continue;
					  if (k==i) continue;

					  // continue for coarse nodes
					  const MathVector<3> &kcoord =  lpositions[k];
					  const int kmatch = num_ijk_coord_matches(nbcoords[nb], kcoord);

					  if (kmatch > maxjmatch)
					  {
					  		maxjmatch = kmatch;

					  		bestcoord[nb] = kcoord;
					  		bestindex[nb] = k;
					  		bestmatch[nb] = kmatch;


					  		if (bestmatch[nb]==3)
					  		{
					  			std::cout << i << "," << k << ":" << kcoord << std::endl;
					  		}
					  }

				  }
			  }




		  }

		  // iterate coarse parents
		  for (int nb=0; nb<nparents; ++nb)
		  {
			  // if (bestindex[nb]==i) continue;
			  UG_ASSERT(bestindex[nb]>=0, "Index should be positive");

			  if (bestmatch[nb]==3)
			  {
				  // all coordinates match
				  // thus, we found existing coarse nodes
				  // and store full matches
				  if (!S.is_isolated(bestindex[nb]))
					  P(i, bestindex[nb]) = nbweights[nb];
			  }
			  else if (m_bCreateVirtualNodes)
			  {
				  // no coarse node found in matrix graph,
				  // thus we create a virtual one
				  std::pair<TVirtualCoarseMap::iterator, bool> search
				  	  = vmap.insert(TVirtualCoarseMapPair(nbcoords[nb], vindex));
				  if (search.second)
				  {
					  // added a new node
					  std::cout << "VADD: " << i << " "<< nbcoords[nb] << " "<< vindex << std::endl;
					  vPositions.push_back(nbcoords[nb]);
					  vindex++;
				  }

				  // add entry to matrix entry list
				  poffd.push_back(pentry_type(i, search.first->second, nbweights[nb]));
			  }
		  }
		}
	} // for all nodes i

	if (m_bCreateVirtualNodes)
	{
		// DEBUG: all entries have been filled
		for (TVirtualCoarseMap::const_iterator it = vmap.begin(); it!=vmap.end(); ++it)
		{
			std::cout << "VCoarse: " << it->first << "->" << it->second <<std::endl;
		}

		/////////////////////////////////////////
		// extend by virtual coarse nodes
		/////////////////////////////////////////
		const int num_vcoarse = vcount[level];
		const int num_extended = N+vcount[level];

		// extend nodes vector & parallel nodes
		// TODO: extend layouts
		nodes.extend(num_vcoarse);

		// extend position data in helper
		std::vector<cAMG_helper::position_type> &mypositions = helper.positions[level];
		THROW_IF_NOT_EQUAL(vPositions.size(), (size_t)num_vcoarse);
		mypositions.insert(mypositions.end(), vPositions.begin(), vPositions.end());
		THROW_IF_NOT_EQUAL(mypositions.size(), (size_t)num_extended);

		// extend matrix
		P.resize_and_keep_values(num_extended, num_extended);
		for (std::vector<pentry_type>::iterator it= poffd.begin();
				it!= poffd.end(); ++it)
		{
			std::cout << "P: " <<  boost::get<0>(*it) << ", " <<  N+boost::get<1>(*it)  << ", " <<  boost::get<2>(*it)  <<std::endl;
			P(boost::get<0>(*it), N+boost::get<1>(*it)) =  boost::get<2>(*it);
		}

	}

	if(bTiming) UG_DLOG(LIB_ALG_AMG, 0, "took " << SW.ms() << " ms");
}



template<typename algebra_type>
void SemiAMG<algebra_type>::user_defined_restriction(const cgraph &S, const matrix_type &A, const prolongation_matrix_type &P, prolongation_matrix_type &R, size_t level, AMGNodes &nodes)
{

	UG_DLOG(LIB_ALG_AMG, 1, std::endl << "restriction... ");


		PROFILE_BEGIN(AMGSetAsTransposeOf)
		// construct restriction R = I_{h -> 2h}
		R.set_as_transpose_of(P);
		PROFILE_END();
		// R is already defragmented

	#ifdef UG_PARALLEL
		R.set_storage_type(PST_CONSISTENT);
	#endif



}

template<typename TAlgebra>
void SemiAMG<TAlgebra>::tostring() const
{
	AMGBase<TAlgebra>::tostring();

	UG_LOG("Semi AMG (i.e. user defined) Preconditioner:\n");


#ifdef UG_PARALLEL
	//if(m_pParallelCoarsening)
	//{	UG_LOG("Parallel Coarsening: " << m_pParallelCoarsening->tostring() << "\n"); }
	//else
	//{	UG_LOG("Parallel Coarsening not set.\n"); }
#endif

}




///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// createAMGLevel:
//-------------------------
/**
 * create AMG matrix R, P, and AH = R A P
 * \param AH
 * \param R
 * \param A
 * \param P
 * \param level
 */
template<typename TAlgebra>
void SemiAMG<TAlgebra>::create_AMG_level(matrix_type &AH, prolongation_matrix_type &R, const matrix_type &A,
		prolongation_matrix_type &P, size_t level, AMGLayouts &nextLevelLayouts)
{
	AMG_PROFILE_FUNC();

#ifdef UG_PARALLEL
	//UG_ASSERT(pcl::NumProcs()==1, "not implemented for procs > 1");

	matrix_type AOL1;
	matrix_type &A_ = const_cast<matrix_type&> (A);
	ParallelNodes PN(A_.layouts(), A_.num_rows());

	stdvector<IndexLayout> vMasterLayouts, vSlaveLayouts;
	MakeFullRowsMatrix(A, AOL1, PN, vMasterLayouts, vSlaveLayouts /*,m_pParallelCoarsening,*/);
	
#else
	const matrix_type &AOL1 = A;
#endif
	prolongation_matrix_type PoldIndices;

	const size_t N = A.num_rows();
	const size_t NOL1 = AOL1.num_rows();
#ifdef UG_PARALLEL
	AMGNodes nodes(NOL1, PN);
#else
	AMGNodes nodes(NOL1);
#endif
	//UG_LOG("NOL1 = " << NOL1 << "\n");

	bool bTiming=true;
	UG_DLOG(LIB_ALG_AMG, 0, "Creating level " << level << ". (" << N << " nodes, " << NOL1-N << " overlapping )" << std::endl << std::fixed);
	//AMG_LOG(1, "Creating level " << level << ". (" << N << " nodes, " << NOL1-N << " overlapping )" << std::endl << std::fixed);
	Stopwatch SW;


	/////////////////////////////////////////
	// build graph
	/////////////////////////////////////////
	cgraph graphS(NOL1);
	//CreateStrongConnectionGraph(AOL1, graphS, m_dTheta);
	CreateFullConnectionGraph(AOL1, graphS, m_dTheta);
	SetIsolatedCoarse(graphS, nodes);

	/////////////////////////////////////////
	// Coarsen
	/////////////////////////////////////////
	UG_DLOG(LIB_ALG_AMG, 0, std::endl << "Coarsening... ");
	if(bTiming) SW.start();

	const size_t nnodes = this->user_defined_coarsening(graphS, level, nodes);

	(void) nnodes; // remove "unused" warning.

	//UG_ASSERT(N == nnodes, "Something wrong!");

	UG_DLOG(LIB_ALG_AMG, 0, nodes.get_nr_of_coarse() << " coarse nodes. ");

	if(bTiming) UG_DLOG(LIB_ALG_AMG, 1, "took " << SW.ms() << " ms");

	/////////////////////////////////////////
	// construct prolongation P = I_{2h->h}
	/////////////////////////////////////////
	UG_DLOG(LIB_ALG_AMG, 0, std::endl << "prolongation... ");

	if(bTiming) SW.start();
	this->user_defined_prolongation(graphS, AOL1, PoldIndices, level, nodes);
#ifdef UG_PARALLEL
	//	this->communicate_prolongation_additive_to_consistent(PN, PoldIndices, true);
#endif




	/////////////////////////////////////////
	/// create layouts
	/////////////////////////////////////////
	{
		stdvector<int> newIndex;
#ifdef UG_PARALLEL
	this->parallel_process_prolongation(PoldIndices, P, 0.0, level, nodes,
					    //	PN, false, nextLevelLayouts, newIndex);
					    PN, true, nextLevelLayouts, newIndex);
#else
	this->serial_process_prolongation(PoldIndices, P, 0.0, level, nodes, newIndex);
#endif
	}

	/////////////////////////////////////////
	// construct restriction R = I_{h->2h}
	/////////////////////////////////////////
	if(bTiming) SW.start();
	this->user_defined_restriction(graphS, AOL1, P, R, level, nodes);
	if(bTiming) UG_DLOG(LIB_ALG_AMG, 1, "took " << SW.ms() << " ms");

	/////////////////////////////////////////
	// create Galerkin product
	/////////////////////////////////////////

	UG_DLOG(LIB_ALG_AMG, 1, "\ngalerkin product... ");
	if(bTiming) SW.start();

	// AH = R A P
	PROFILE_BEGIN(AMGCreateAsMultiplyOf)
		CreateAsMultiplyOf(AH, R, A, P, m_dGalerkinTruncation);
	PROFILE_END();

	if(bTiming) UG_DLOG(LIB_ALG_AMG, 1, "took " << SW.ms() << " ms");

	// defragment
	if(bTiming) SW.start();
	AH.defragment();
	if(bTiming) UG_DLOG(LIB_ALG_AMG, 1, std::endl << "Defragment ... took " << SW.ms() << " ms");

#ifdef UG_PARALLEL
	AH.set_storage_type(PST_ADDITIVE);
#endif


	UG_DLOG(LIB_ALG_AMG, 1, "\n");

	if(m_writeMatrices)
	{
		AMGBase<TAlgebra>::write_debug_matrix_markers(level, nodes);
		this->write_debug_matrices(AH, R, A, P, level);
	}
}



// end group plugins_amg_rsamg
/// \}

} // namespace ug

#endif //  __H__UG__LIB_ALGEBRA_SEMIAMG_SOLVER__SEMIAMG_IMPL_H__
