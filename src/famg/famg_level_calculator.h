/**
 * \file famg_level_calculator.h
 *
 * \author Martin Rupp
 *
 * \date 7.12.2010
 *
 * implementation file for famg
 *
 * Goethe-Center for Scientific Computing 2010-2012.
 *
 */


#ifndef __H__LIB_ALGEBRA__AMG__FAMG_LEVEL_CALCULATOR_H__
#define __H__LIB_ALGEBRA__AMG__FAMG_LEVEL_CALCULATOR_H__

//#include "sparsematrix_util.h"

#include "common/assert.h"
#include "common/log.h"

#include "famg.h"
#include "lib_algebra/lib_algebra.h"
#include "lib_algebra/cpu_algebra_types.h"
#include "lib_algebra/common/stl_debug.h"
#include "../amg_debug_helper.h"

#include "common/stopwatch.h"
#include "common/assert.h"
#include "lib_algebra/common/heap/maxheap.h"
#include "famg.h"
#include <set>
#include <string>
#include <stack>
#include "../amg_misc.h"

#ifdef UG_PARALLEL
#include "lib_algebra/parallelization/parallel_matrix_overlap_impl.h"
#include "lib_algebra/parallelization/communication_scheme.h"
#endif

#include "lib_algebra/common/connection_viewer_output.h"
#include "common/util/file_util.h"

/// \addtogroup plugins_amg_famg

//const double damping = 0.5;

#if 1

//#define AMG_WRITE_GRAPH


//#define AMG_PRINT_INDIRECT

#define FAMG_PRINT_GRAPH

#define FAMG_PRINT_COARSENING
#define FAMG_PRINT_P
#define FAMG_PRINT_R
#define FAMG_PRINT_AH

#define FAMG_PRINT_POSSIBLE_PARENTS
#endif

//#define GRAPH_WITH_LOCAL_INVERSE
#include "famg_structs.h"
#include "famg_interpolation_calculator.h"
#include "famg_nodeinfo.h"

namespace ug{

/// \addtogroup plugins_amg_famg
/// \{




///	algebra blocks are static.
template <> struct block_traits<MathVector<3> >
{
	enum{
		is_static = 1
	};
};


// CreateSymmConnectivityGraph:
//------------------------------
/**
 * Creates a graph which has a connection from i to j if A_{ij} != 0 or A_{ji} != 0
 * \param A						the matrix A
 * \param SymmNeighGraph		later used to determine which neighbors' rating needs to be updated
 */
template<typename matrix_type>
static void CreateSymmConnectivityGraph(const matrix_type &A, cgraph &SymmNeighGraph)
{
	AMG_PROFILE_FUNC();

	for(size_t i=0; i<A.num_rows(); i++)
		for(typename matrix_type::const_row_iterator conn = A.begin_row(i); conn != A.end_row(i); ++conn)
		{
			if(i == conn.index() || conn.value() == 0.0) continue;
			SymmNeighGraph.set_connection(i, conn.index());
			SymmNeighGraph.set_connection(conn.index(), i);
		}

#ifdef FAMG_PRINT_GRAPH

	IF_DEBUG(DBG_FAMG_CON_GRAPH, 3)
	{
		UG_LOG("\nSymmNeighGraph:\n");
		SymmNeighGraph.print();
		UG_LOG("\n\n");
	}
#endif
}

template<typename TAlgebra>
class FAMGLevelCalculator
{
public:
//	Algebra type
	typedef TAlgebra algebra_type;

//	Vector type
	typedef typename TAlgebra::vector_type vector_type;

//	Matrix type
	typedef typename TAlgebra::matrix_type matrix_type;

	typedef typename matrix_type::value_type value_type;

	//typedef typename CPUAlgebra::matrix_type prolongation_matrix_type;
	typedef typename TAlgebra::matrix_type prolongation_matrix_type;

private:
	// refernce to famg object. this will be done differently in the future
	// (constructor looks ugly)
	FAMG<TAlgebra> &m_famg;

	// the coarse matrices to be filled
	matrix_type &AH;

	// the fine grid matrix (const!)
	const matrix_type &A;

	// the fine grid matrix (const!)
	const matrix_type &Aorig;

	// interpolation/restriction matrices to be filled
	prolongation_matrix_type &R, &PnewIndices;

	int m_verbose;

	size_t level; // current amg level

#ifdef UG_PARALLEL
	ParallelNodes PN;
	//! layout for sending/receiving coarsening information
	//! note: this is not a classical master/slave distribution!
	IndexLayout OLCoarseningReceiveLayout, OLCoarseningSendLayout;

	//! processes with higher/lower color to distribute coarsening information
	std::vector<int> processesWithLowerColor, processesWithHigherColor;
	int m_myColor;

	// overlap 2 matrix in parallel case
	matrix_type A_OL2;
	IndexLayout OL1MasterLayout, OL1SlaveLayout;
	IndexLayout OL2MasterLayout, OL2SlaveLayout;
	std::vector<size_t> overlapSize;
	std::vector<IndexLayout> masterLayouts;
	std::vector<IndexLayout> slaveLayouts;
#else
	const matrix_type &A_OL2;
#endif
	FAMGNodes<prolongation_matrix_type> rating;
	AMGLayouts &nextLevelLayouts;
	stdvector< vector_type > &m_testvectors;

	// the calculator which calculates our interpolation weighs
	FAMGInterpolationCalculator<TAlgebra> calculator;

	bool m_bTiming; //< if true, timing will be outputed

	// this structure holds coarse/fine information as well as
	// master/slave. it tries to replace the missing "node" object here


	//! list of possible interpolating parents for each node
	typedef typename prolongation_matrix_type::value_type smallmat_ptype;
	typedef neighborstruct_var<smallmat_ptype> neighT;
	stdvector<stdvector< neighT > > possible_parents;

	stdvector<bool> prolongation_calculated;

	//! heap used for sorting of ratings
	maxheap<FAMGNode> heap;

	//! used to determine which neighbors' rating needs to be updated:
	cgraph SymmNeighGraph;

	stdvector<double> fvalues;
	// our testvectors


	prolongation_matrix_type PoldIndices;

	void create_restriction(stdvector<int> &newIndex)
	{
		AMG_PROFILE_FUNC();
		Stopwatch SW;
		if(m_bTiming) {UG_LOG("\nrestriction... "); SW.start();}
		// construct restriction R = I_{h -> 2h}

		if(m_famg.m_bInjectionRestriction==false)
			R.set_as_transpose_of(PnewIndices);
		else
		{
			R.resize_and_clear(PnewIndices.num_cols(), PnewIndices.num_rows());
			for(size_t i=0; i<A.num_rows(); i++)
			{
				if(newIndex[i] >= 0)
				{
					UG_ASSERT(rating[i].is_coarse(), i);
					R(newIndex[i], i) = 1.0;
				}
				else
				{
					UG_ASSERT(!rating[i].is_coarse(), i);
				}

			}
		}


		// R is already defragmented
		if(m_bTiming) { UG_LOG("took " << SW.ms() << " ms");}

	#ifdef UG_PARALLEL
		R.set_storage_type(PST_CONSISTENT);
	#endif

	#ifdef FAMG_PRINT_R
		IF_DEBUG(LIB_ALG_AMG, 4)
		{
			UG_DLOG(LIB_ALG_AMG, 0, std::endl << "Restriction level " << level << std::endl);
			R.print();
		}
	#endif
	}

	void create_galerkin_product()
	{
		AMG_PROFILE_FUNC();
		Stopwatch SW;
		if(m_bTiming) { UG_LOG("Galerkin Product..."); SW.start(); }
		// AH = R A P

		/* matrix_type AH2;
		matrix_type RoldIndices;
		RoldIndices.set_as_transpose_of(PoldIndices);
		CreateAsMultiplyOf2(AH2, RoldIndices, A, PoldIndices); */

		// R.print("R");
		// A.print("A");
		// PnewIndices.print("P");
		CreateAsMultiplyOf(AH, R, Aorig, PnewIndices, m_famg.get_galerkin_truncation());
		// AH.print();
		if(m_bTiming) {UG_LOG("took " << SW.ms() << " ms");}

		// finalize
		if(m_bTiming) { UG_LOG("\nFinalizing.."); SW.start(); }
		AH.defragment();
		if(m_bTiming) {UG_LOG("took " << SW.ms() << " ms");}

		if(m_famg.m_bChecks)
		{
			bool b = CheckDiagonalInvertible(AH);
			if(b==false)
			UG_LOG("ERROR: FAMG level " << level << ": Diag(AH) not invertible!!!");
		}

		PROFILE_END();

	#ifdef FAMG_PRINT_AH
		IF_DEBUG(LIB_ALG_AMG, 3)
		{
			UG_DLOG(LIB_ALG_AMG, 0, "AH level " << level << std::endl);
			AH.print();
		}
	#endif
	}
public:
	// todo: clean up this mess of constructor
	FAMGLevelCalculator(FAMG<TAlgebra> &f,
			matrix_type &_AH, prolongation_matrix_type &_R,
			const matrix_type &_A, const matrix_type &_Aorig,
			prolongation_matrix_type &_P, size_t _level,
			stdvector< vector_type > &testvectors, stdvector<double> &omega, AMGLayouts &_nextLevelLayouts)
	: m_famg(f), AH(_AH), A(_A), Aorig(_Aorig), R(_R), PnewIndices(_P), level(_level),
#ifdef UG_PARALLEL
			PN(A.layouts(), A.num_rows()),
			rating(PoldIndices, _level, m_famg.m_amghelper, PN),
#else
			A_OL2(A),
			rating(PoldIndices, _level, m_famg.m_amghelper),
#endif
			nextLevelLayouts(_nextLevelLayouts),
			m_testvectors(testvectors),
			calculator(A, A_OL2, m_famg, testvectors, omega, fvalues)
	{
		m_verbose = f.m_verbose;
	}




public:
	void do_calculation()
	{
		AMG_PROFILE_FUNC();
		m_bTiming=false;
		Stopwatch SW, SWwhole; SWwhole.start();

		// 1. Overlap calculation
		//---------------------------
#ifdef UG_PARALLEL
		create_OL2_matrix();
#else
		rating.create(A_OL2.num_rows());
#endif
		size_t N = A_OL2.num_rows();

//		if(m_famg.m_bWriteFValues)
			fvalues.resize(N, 0.0);
//		else 						fvalues.clear();

		// 2. global Testvector calculation (damping)
		//-----------------------------------------------

		smooth_testvectors();
		get_testvectors_on_OL2();
		calculator.init();


		// 3. create heap, P, SymmNeighGraph
		//-------------------------------------

		bool bExternalCoarsening = m_famg.m_bExternalCoarsening;
		bool bUsePrecalculate = m_famg.m_bUsePrecalculate;
		if(bExternalCoarsening == false && m_famg.m_fileCoarsening != m_famg.READ_COARSENING)
		{
			heap.create(rating.nodes);

			if(m_bTiming) { UG_LOG("\nCreate P, SymmNeighGraph... "); SW.start(); }
			PoldIndices.resize_and_clear(N, N);
			// get neighboring information
			SymmNeighGraph.resize(N);
			CreateSymmConnectivityGraph(A_OL2, SymmNeighGraph);
			if(m_bTiming) { UG_LOG("took " << SW.ms() << " ms"); }

			if(bUsePrecalculate)
				// get possible parent nodes
				calculate_all_possible_parent_pairs();

#ifdef UG_PARALLEL
			// 4. coloring in parallel
			//-------------------------------------
			color_process_graph();
			receive_coarsening_from_processes_with_lower_color();
#endif
			rating.calculate_unassigned();

			// 5. do coarsening
			//-------------------------
			if(bUsePrecalculate)
			{
				// calculate ratings (not precalculateable because of coarse/uninterpolateable)
				if(m_bTiming) { UG_LOG("\ncalculate ratings... "); SW.start(); }
				GetRatings(possible_parents, rating, heap);
				if(m_bTiming) { UG_LOG("took " << SW.ms() << " ms");}
				//heap.print();

				// do coarsening
				precalculate_coarsening();
			}
			else
				on_demand_coarsening();


			// 7. do aggressive coarsening
			//-----------------------------------------------
			// this has to be BEFORE communication of coarse/fine markers
			//

			if(m_bTiming) {UG_LOG("\nsecond coarsening... "); SW.start();}
			if(m_famg.get_aggressive_coarsening() == true)
				get_aggressive_coarsening_interpolation();
			else
				set_uninterpolateable_as_coarse();
			if(m_bTiming) { UG_LOG("took " << SW.ms() << " ms"); }



			#ifdef UG_PARALLEL
			// 6. send coarsening data to processes with higher color
			//----------------------------------------------------------
				send_coarsening_data_to_processes_with_higher_color();
				// we could use this information to do better interpolation on the border
				// send_coarsening_data_to_processes_with_lower_color_and_receive_coarsening_data_from_processes_with_higher_color

				// 7. calculate interpolation of nodes with uncalculated_fine
				//calculate_uncalculated_fine_nodes();
			#endif

			// 8.

			// resize P, output statistics & debug stuff

			UG_DLOG(LIB_ALG_AMG, 1, std::endl << N - rating.get_unassigned() << " nodes assigned, " << rating.get_nr_of_coarse() << " coarse, "
					<< N - rating.get_unassigned() - rating.get_nr_of_coarse() << " fine, " << rating.get_unassigned() << " unassigned.");

			// set P to real size
		}
		else // = !(bExternalCoarsening == false && m_famg.m_fileCoarsening != m_famg.READ_COARSENING)
		{
			if(m_famg.m_writeMatrices) m_famg.write_debug_matrix_markers(level, rating);
			PoldIndices.resize_and_clear(N, N);

			UG_DLOG(DBG_AMG_COARSENING, 1, "XC\n");
			if(m_famg.m_fileCoarsening == m_famg.READ_COARSENING)
				read_coarsening();
			else
				rs_amg_external_coarsening();

			UG_DLOG(DBG_AMG_COARSENING, 1, "external_coarsening_calculate_prolongation...\n");

			external_coarsening_calculate_prolongation();

			rs_amg_external_postcoarsening();
			set_uninterpolateable_as_coarse();
		}

		if(m_famg.m_fileCoarsening == m_famg.WRITE_COARSENING)
			write_coarsening();

		set_uninterpolateable_as_coarse();
		make_coarsening_consistent();

//		IF_DEBUG(LIB_ALG_AMG, 3)
//		{	rating.print(); 	}

//		Write_LE_Error<prolongation_matrix_type, vector_type>
//			(PoldIndices, m_famg.m_amghelper.positions[level], m_famg.m_writeMatrixPath, "FAMG_tv", level);
		{
			stdvector<int> newIndex;
			create_new_indices(newIndex);
			//for(size_t i=0; i<rating.get_nr_of_coarse(); i++) { UG_ASSERT(amghelper.parentIndex[level+1][i] != -1, i << " == -1???"); UG_LOG(i << " = " << amghelper.parentIndex[level+1][i] << "\n"); }

			//UG_DLOG(LIB_ALG_AMG, 0, std::endl << rating.get_nr_of_coarse() << " coarse.");

			// 8. construct restriction R = I_{h->2h}
			//-----------------------------------------
			create_restriction(newIndex);
		}


		for(size_t i=0; i<A.num_rows(); i++)
		{
			if(PoldIndices.num_connections(i) > 1 && (rating[i].is_fine() || rating[i].is_fine_direct()) == false)
			{	UG_LOG("level " << level << " node " << i << "\n");		}
		}
		// 9. create Galerkin product AH = R A P
		//-----------------------------------------

		create_galerkin_product();
		// 10.
		//-----------------------------------------
		// [ debug output
		if(m_famg.m_writeMatrices)
		{
			m_famg.write_debug_matrix_markers(level, rating);
			m_famg.write_debug_matrices(AH, R, A, PnewIndices, level);
		}
		// ]

		calculate_next_testvectors();

#ifdef UG_PARALLEL
		IF_DEBUG(LIB_ALG_AMG, 3)
		{
			PRINTLAYOUT(A_OL2.layouts()->proc_comm(), A_OL2.layouts()->comm(), AH.layouts()->master(), AH.layouts()->slave());
		}
#endif

		if(m_famg.m_bWriteFValues)
		{
			std::string filename = m_famg.GetAMGFilenamePath(A, level, "fvalues", level, ".vec", true);
			ConnectionViewer::WriteVector(filename, fvalues, &m_famg.m_amghelper.positions[level][0], m_famg.m_amghelper.dimension);
		}
	}

private:

	void create_new_indices(stdvector<int> &newIndex)
	{
		AMG_PROFILE_FUNC();
#ifdef UG_PARALLEL
		PoldIndices.set_layouts(A.layouts());
		m_famg.parallel_process_prolongation(PoldIndices, PnewIndices,
				m_famg.m_dProlongationTruncation, level, rating, PN, false, nextLevelLayouts, newIndex);
		if(m_famg.m_bChecks)
			TESTLAYOUTS((&nextLevelLayouts));
#else
		m_famg.serial_process_prolongation(PoldIndices, PnewIndices, m_famg.m_dProlongationTruncation, level, rating, newIndex);
#endif

	}

	//	get_aggressive_coarsening_interpolation
	//! tries for all nodes which are uninterpolateable an indirect interpolation
	//! \sa set_uninterpolateable_as_coarse
	void get_aggressive_coarsening_interpolation()
	{
		// ?????
		AMG_PROFILE_FUNC();

		for(size_t i=0; i<A.num_rows(); i++)
		{
			//UG_ASSERT(rating[i].is_valid_rating() == false, "node " << i << " has valid rating, but has not been processed yet?");
			if(rating[i].is_uninterpolateable() == false || rating.i_must_assign(i) == false) continue;

			// we have problems with nodes on the interface:
			// suppose i is a master node
			// 1. we would need the prolongation of the neighbors of i, which are slave0 nodes
			// 2. we could then indirectly interpolate from slave1 and master2 nodes
			// 3. on the slave side, we do not have the slave2 nodes, so we would have to add them somehow
			if(rating.is_inner_node(i))
				calculator.get_all_neighbors_interpolation(i, PoldIndices, rating, true);
			else
				rating.set_coarse(i);
		}
	}


	//	set_uninterpolateable_as_coarse
	//! sets all nodes which are uninterpolateable as coarse
	//! \sa get_aggressive_coarsening_interpolation
	void set_uninterpolateable_as_coarse()
	{
		AMG_PROFILE_FUNC();

		for(size_t i=0; i<rating.size(); i++)
		{
			if(rating[i].is_uninterpolateable() == false || rating.i_must_assign(i) == false) continue;
			//if(rating[i].is_fine() || rating.i_must_assign(i) == false) continue;
			rating.set_coarse(i);
		}
		for(size_t i=0; i<rating.size(); i++)
		{
			if(rating[i].is_coarse())
				PoldIndices(i,i)=1.0;
		}
	}

	void make_coarsening_consistent();

	void on_demand_coarsening();
	void precalculate_coarsening();
	void calculate_all_possible_parent_pairs();

#ifdef UG_PARALLEL
	void color_process_graph();
	void receive_coarsening_from_processes_with_lower_color();
	void send_coarsening_data_to_processes_with_higher_color();
	void create_OL2_matrix();
	void calculate_uncalculated_fine_nodes();
	void rs_amg_coarsening_new(const cgraph &graphS, const cgraph &graphST, nodeinfo_pq_type &PQ, AMGNodes &nodes);
	void rs_amg_coarsening_LE_griebel(const cgraph &graphS, const cgraph &graphST, nodeinfo_pq_type &PQ, AMGNodes &nodes);
#endif
	void external_coarsening_calculate_prolongation();

	void rs_amg_external_coarsening();
	void write_XC_coarsening(cgraph &graphS, AMGNodes &nodes, std::string addStr);
	void rs_amg_external_postcoarsening();
	void get_external_coarsening_graph(cgraph &graphS);

	void smooth_testvectors();
	void get_testvectors_on_OL2();
	void calculate_next_testvectors();

public:
	void write_coarsening(), read_coarsening();

	void onlyTV()
	{
#ifdef UG_PARALLEL
		create_OL2_matrix();
#endif
		smooth_testvectors();
		get_testvectors_on_OL2();
		calculate_next_testvectors();
	}
};

// end group plugins_amg_famg
/// \}

} // namespace ug

#ifdef UG_PARALLEL
#include "famg_parallel_coarsening_impl.h"
#endif

#include "famg_on_demand_coarsening_impl.h"
#include "famg_precalculate_coarsening_impl.h"
#include "famg_testvectors.h"
#include "famg_external_coarsening.h"
#endif //  __H__LIB_ALGEBRA__AMG__FAMG_IMPL_H__
