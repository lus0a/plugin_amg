/**
 *
 * \author Martin Rupp
 *
 * \date 03.12.2009
 *
 * implementation file for amg
 *
 * Goethe-Center for Scientific Computing 2009-2010.
 */

#ifndef __H__UG__LIB_ALGEBRA__SEMIAMG_SOLVER__SEMIAMG_IMPL_H__
#define __H__UG__LIB_ALGEBRA__SEMIAMG_SOLVER__SEMIAMG_IMPL_H__


#include "semi_amg.h"
#include "lib_algebra/common/graph/graph.h"



#include "common/stopwatch.h"

#include "../amg_debug_helper.h"
#include "../amg_misc.h"
#include "../amg_profiling.h"
#include "../amg_parallel.h"

#include "../rsamg/rsamg_nodeinfo.h"
#include "../rsamg/strong_connections.h"

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
: AMGBase<TAlgebra>(), m_tli(), m_dTheta(0.25)
{
#ifdef UG_PARALLEL

#endif
}



template<typename algebra_type>
size_t SemiAMG<algebra_type>::user_defined_coarsening(const cgraph &strong, size_t level, AMGNodes &nodes)
{
	cAMG_helper &helper= m_amghelper;

	const size_t nnodes = helper.positions[level].size();
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
	UG_LOG("\nUser defined C/F-splitting on level " << level <<": " << nfine << " fine nodes out of "<< nnodes);

	return nnodes;
}


template<typename algebra_type>
void SemiAMG<algebra_type>::user_defined_prolongation(const cgraph &strong, const matrix_type &A, prolongation_matrix_type &P, size_t level, AMGNodes &nodes)
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

	// first pass
	std::set<size_t> coarse;
	const size_t nnodes = helper.positions[level].size();
	size_t nskipped = 0;
	for(size_t i=0; i<nnodes; i++)
	{
		PROGRESS_UPDATE(prog, i);

		// map coarse grid points by identity
		if (nodes[i].is_coarse())
		  {
		    P(i, i) = 1.0; continue;
		    //	std::cout << "Ignoring " << i << std::endl;
		  }

		// continue for fine nodes here

		// skip off-proc nodes
		if(!nodes.is_master_or_inner(i))
		{
		  //std::cout << "Skipping " << i << std::endl;
			nskipped++;
			continue;
		}

		coarse.clear();

		if (0) {

		  
		std::cout << "Setting " << i;
		for(typename matrix_type::const_row_iterator conn = A.begin_row(i); conn != A.end_row(i); ++conn)
		{
			size_t j = conn.index();
			if (nodes[j].is_coarse() && strong.has_connection(i,j))
			{ coarse.insert(j);}
		}

		// nodes without neighbors are coarse
		const size_t ncoarse=coarse.size();
		if (ncoarse == 0)
		{
			// treat in second step
			continue;
		}

		// determine interpolation weights
		MathVector<3> &coord =  helper.positions[level][i];
		size_t flevels[3];
		m_tli.template get_node_levels<3>(coord, flevels);

		double weight = 1.0;
		for (size_t j=0; j<3; j++)
		{
			if (flevels[j] == level) weight*=0.5;
		}

		
		// set weight
		for(std::set<size_t>::iterator citer=coarse.begin(); citer!=coarse.end(); citer++)
		{
			std::cout << ":" << *citer << std::endl;
			P(i, *citer) = 1.0/ncoarse;//weight;
		}
		std::cout <<std::endl;




		} else {

		  //std::cout << "SettingXCP for  " << i;
		  MathVector<3> &coord =  helper.positions[level][i];
		  MathVector<3> nbcoords[27];
		  double nbweights[27];
		  const int nparents = m_tli.get_parents(coord, level, nbcoords, nbweights);
		  UG_ASSERT(nparents<=8, "Huhh: Too many neighbors???");
		 
		  double sumw = 0.0;
		  for(typename matrix_type::const_row_iterator conn = A.begin_row(i); conn != A.end_row(i); ++conn)
		    {
		      	const size_t j = conn.index();
			if (! nodes[j].is_coarse()) continue;
			
			MathVector<3> &coordj =  helper.positions[level][j];
			for (int jj=0;jj<nparents; ++jj)
			  {
			    if ((nbcoords[jj][0]==coordj[0]) &&
				(nbcoords[jj][1]==coordj[1]) &&
				(nbcoords[jj][2]==coordj[2]) )
			      {
				P(i, j) = nbweights[jj];
				//		std::cout << ":" << j  << "=" << nbweights[jj]<< std::endl;
				sumw += nbweights[jj];
				break;
			      }
			  }

			

		    }
		  //std::cout <<"=" << sumw << std::endl;
		  //UG_ASSERT(sumw<=1.0, "Huhh: Large sum???");
		}



	}
/*
	// second pass (eliminate)
	for(size_t i=0; i<nnodes; i++)
	{
		if (nodes[i].is_coarse()) continue;

		// continue for fine nodes here
		coarse.clear();

		for(typename matrix_type::const_row_iterator conn = A.begin_row(i); conn != A.end_row(i); ++conn)
		{
			size_t j = conn.index();
			if (nodes[j].is_coarse() && strong.has_connection(i,j))
			{ coarse.insert(j);}
		}

		// nodes without neighbors are caorse
		const size_t ncoarse=coarse.size();
		if (ncoarse == 0)
		{

			nodes.external_set_coarse(i);
			P(i, i) = 1.0;
			continue;
		}
	}
	*/
	// second pass (indirect interpolation)
	/*for(size_t i=0; i<N; i++)
	{
		// indirect interpolation
		PROGRESS_UPDATE(prog, N+i);

		coarse.clear();
		for(typename matrix_type::const_row_iterator conni = A.begin_row(i);
				conni != A.end_row(i); ++conni)
		{
			const size_t j = conni.index();

			for(typename matrix_type::const_row_iterator connj = A.begin_row(j);
							connj != A.end_row(j); ++connj)
			{
				const size_t k = connj.index();
				if (nodes[k].is_coarse()) { coarse.push_back(k);}
			}
		}

		const size_t ncoarse=coarse.size();
		if (ncoarse == 0)
		{ nodes.set_coarse(i); continue; }

		for(size_t k=0; k<ncoarse; k++)
		{ P(i, coarse[k]) = 1.0/ncoarse; }

	}
*/
	if(bTiming) UG_DLOG(LIB_ALG_AMG, 0, "took " << SW.ms() << " ms");
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
		prolongation_matrix_type &PnewIndices, size_t level, AMGLayouts &nextLevelLayouts)
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
	CreateStrongConnectionGraph(AOL1, graphS, m_dTheta);
	SetIsolatedCoarse(graphS, nodes);

	/////////////////////////////////////////
	// Coarsen
	/////////////////////////////////////////
	UG_DLOG(LIB_ALG_AMG, 0, std::endl << "Coarsening... ");
	if(bTiming) SW.start();

	const size_t nnodes = this->user_defined_coarsening(graphS, level, nodes);

	(void) nnodes; // remove unused warning.

	//UG_ASSERT(N == nnodes, "Something wrong!");

	UG_DLOG(LIB_ALG_AMG, 0, nodes.get_nr_of_coarse() << " coarse nodes. ");

	if(bTiming) UG_DLOG(LIB_ALG_AMG, 1, "took " << SW.ms() << " ms");

	/////////////////////////////////////////
	// Agressive Coarsening
	/////////////////////////////////////////



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
	this->parallel_process_prolongation(PoldIndices, PnewIndices, 0.0, level, nodes,
					    //	PN, false, nextLevelLayouts, newIndex);
					    PN, true, nextLevelLayouts, newIndex);
#else
	this->serial_process_prolongation(PoldIndices, PnewIndices, 0.0, level, nodes, newIndex);
#endif
	}

	/////////////////////////////////////////
	// construct restriction R = I_{h->2h}
	/////////////////////////////////////////

	UG_DLOG(LIB_ALG_AMG, 1, std::endl << "restriction... ");
	if(bTiming) SW.start();

	PROFILE_BEGIN(AMGSetAsTransposeOf)
	// construct restriction R = I_{h -> 2h}
	R.set_as_transpose_of(PnewIndices);
	PROFILE_END();
	// R is already defragmented

#ifdef UG_PARALLEL
	R.set_storage_type(PST_CONSISTENT);
#endif

	if(bTiming) UG_DLOG(LIB_ALG_AMG, 1, "took " << SW.ms() << " ms");

	/////////////////////////////////////////
	// create Galerkin product
	/////////////////////////////////////////

	UG_DLOG(LIB_ALG_AMG, 1, "\ngalerkin product... ");
	if(bTiming) SW.start();

	// AH = R A P
	PROFILE_BEGIN(AMGCreateAsMultiplyOf)
		CreateAsMultiplyOf(AH, R, A, PnewIndices, m_dGalerkinTruncation);
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
		this->write_debug_matrices(AH, R, A, PnewIndices, level);
	}
}



// end group plugins_amg_rsamg
/// \}

} // namespace ug

#endif //  __H__UG__LIB_ALGEBRA_SEMIAMG_SOLVER__SEMIAMG_IMPL_H__
