/**
 * \file amg_base_impl.h
 * \author Martin Rupp
 *
 * \date 01.12.2010
 *
 * implementation file for base amg functionality
 *
 * Goethe-Center for Scientific Computing 2010.
 */

#ifndef __H__UG__LIB_ALGEBRA__AMG_SOLVER__AMG_BASE_IMPL_H__
#define __H__UG__LIB_ALGEBRA__AMG_SOLVER__AMG_BASE_IMPL_H__

#include <fstream>
#include "common/stopwatch.h"
#include "amg_debug.h"
#ifdef UG_PARALLEL
#include "lib_algebra/parallelization/collect_matrix.h"
#include "lib_algebra/parallelization/parallel_vector.h"
#include "lib_algebra/parallelization/parallel_storage_type.h"
#include "lib_algebra/parallelization/consistency_check.h"
#include "lib_algebra/parallelization/serialization.h"
#include "pcl/pcl_tostring.h"
#endif
#include "common/util/file_util.h"
#include "common/progress.h"
#include "common/util/string_concat.h"
#include "additional_math.h"
#include "lib_algebra/common/connection_viewer_output.h"

#include "lib_algebra/operator/preconditioner/gauss_seidel.h"
#include "lib_algebra/operator/linear_solver/lu.h"


//#define UG_DLOGN(tag, level, msg)			UG_DLOG(tag, level, msg << "\n");

std::string GetProcFilename(std::string path, std::string name, std::string extension);

namespace ug{

/// \addtogroup plugins_amg
/// \{

#ifdef UG_PARALLEL
template<typename TMatrix, typename TVector>
void SetParallelVectorAsMatrix(TVector &v, TMatrix &m, ParallelStorageType t)
{
//	v.set_communicator(m.communicator());
//	v.set_process_communicator(m.layouts()->proc_comm());
//	v.set_layouts(m.layouts()->master(), m.layouts()->slave());
	v.set_layouts(m.layouts());
	v.set_storage_type(t);
}
#endif



template<typename TAlgebra>
void AMGBase<TAlgebra>::calculate_level_information(size_t level, double createAMGlevelTiming)
{
	PROFILE_FUNC();
	AMGLevel &L = *levels[level];
	matrix_operator_type &A = *L.pAgglomeratedA;

	if(m_levelInformation.size() < level+1)
		m_levelInformation.resize(level+1);
	LevelInformation &li = m_levelInformation[level];


	if(level > 0)
		li.calculate(A, level, createAMGlevelTiming, &m_levelInformation[level-1]);
	else
		li.calculate(A, 0, 0.0, NULL);


	AMG_LOG(1, li.tostring() << "\n");
	AMG_LOG(1, " level took " << li.m_dCreationTimeMS << " ms" <<
			std::endl << std::endl);
}

template<typename TAlgebra>
void AMGBase<TAlgebra>::init_smoother(size_t level)
{
	AMGLevel &L = *levels[level];
	L.presmoother = m_presmoother->clone();
	L.presmoother->init(L.pAgglomeratedA);

	if(m_presmoother == m_postsmoother)
		L.postsmoother = L.presmoother;
	else
	{
		L.postsmoother = m_postsmoother->clone();
		L.postsmoother->init(L.pAgglomeratedA);
	}
}

template<typename TAlgebra>
bool AMGBase<TAlgebra>::init_smoothers()
{
	try{

		size_t N0 = levels[0]->pA->num_rows();
		PROGRESS_START(prog, N0, "Creating AMG - Init Smoothers");

		for(size_t level=0; level < my_multigrid_levels(); level++)
		{
			AMGLevel &L = *levels[level];
			PROGRESS_UPDATE(prog, N0-L.pA->num_rows());
			init_smoother(level);
		}
		return true;
	}UG_CATCH_THROW("AMGBase::" << __FUNCTION__ << " failed")

	return false;
}


///////////////////////////////////////////////////////////////////////////////////////////////////////////////////

template<typename TAlgebra>
bool AMGBase<TAlgebra>::preprocess(SmartPtr<MatrixOperator<matrix_type, vector_type> > pOp)
{
	try{
	matrix_operator_type& mat = *pOp;

	AMG_PROFILE_FUNC();

	if(m_bOneInit && m_bInited)
		return true;

#ifdef UG_PARALLEL
	if(mat.num_rows() == 0 && mat.layouts()->proc_comm().empty())
		return true;
#endif

	cleanup();

	update_positions();
	m_amghelper.parentIndex = &m_parentIndex;

	if(m_writeMatrices && m_dbgDimension == 0)
	{
		m_writeMatrices = false;
		UG_ERR_LOG("ERROR: AMG set to write Matrices, but no position data given. Not writing grid.\n");
		UG_LOG("ERROR: AMG set to write Matrices, but no position data given. Not writing grid.\n");
	}

	// init m_amghelper for grid printing
	if(needs_positions())
	{
		UG_COND_THROW(m_dbgDimension == 0, "AMG method needs positions");
		if(m_dbgPositions.size() > 0)
		{
			m_amghelper.positions.resize(1);
			m_amghelper.positions[0] = m_dbgPositions;
		}
		m_amghelper.dimension = m_dbgDimension;
	}

	UG_COND_THROW(m_writeMatrices && m_dbgDimension == 0, "you need to provide position information for writing matrices");

	if(m_basesolver.invalid())
	{
		UG_LOGN("amg_base::preprocess(): No base solver selected. Call set_base_solver(basesolver) to set a base solver.");
		return false;
	}
	if(m_presmoother.invalid())
	{
		UG_LOGN("amg_base::preprocess(): No PreSmoother selected. Call set_presmoother(presmoother) to set a PreSmoother.");
		return false;
	}
	if(m_postsmoother.invalid())
	{
		UG_LOGN("amg_base::preprocess(): No PostSmoother selected. Call set_postsmoother(postsmoother) to set a PostSmoother.");
		return false;
	}

#ifdef UG_PARALLEL
	if(pcl::NumProcs() > 1)
	{
		UG_COND_THROW(!m_presmoother->supports_parallel(), "presmoother has to support parallel!");
		UG_COND_THROW(!m_postsmoother->supports_parallel(), "postsmoother has to support parallel!");
	}
#endif

	UG_DLOG(LIB_ALG_AMG, 1, "Starting AMG Setup." << std::endl << std::endl);

	Stopwatch SWwhole;
	SWwhole.start();

	size_t level=0;

	mat.defragment();

	//UG_LOG("======= WARNING: MODIFIED A by A = D^{-1} A =======\n");
	/*for(size_t r=0; r<mat.num_rows(); r++)
	{
		value_type m = mat(r, r);
		for(typename matrix_type::row_iterator it = mat.begin_row(r); it != mat.end_row(r); ++it)
			it.value() /= m;
	}*/

	double createAMGlevelTiming=0;

#ifdef UG_PARALLEL
	m_agglomerateLevel = -1;
#endif

	levels.clear();
	AMGLevel *pL = new AMGLevel(0);
	pL->pA = pOp;
#ifdef UG_PARALLEL
//	pL->com = mat.layouts()->comm();
	pL->spLayouts->comm() = mat.layouts()->comm();
	pL->bHasBeenMerged = false;
#endif
	pL->spYConvCheck = spYConvCheck;
	levels.push_back(make_sp(pL));

	PROGRESS_START(prog, pOp->num_rows(), "Creating AMG");
	for(; ; level++)
	{

		AMGLevel &L = *levels[level];
		PROGRESS_UPDATE(prog, pOp->num_rows()-L.pA->num_rows());

		L.pAgglomeratedA = L.pA;
		calculate_level_information(level, createAMGlevelTiming);

		if(m_levelInformation[level].get_nr_of_nodes() < m_maxNodesForBase
				|| m_levelInformation[level].get_fill_in() > m_dMaxFillBeforeBase
				|| level >= m_maxLevels-1)
		{
			create_direct_solver(level);
			//L.pAgglomeratedA = L.pA;
			L.pAgglomeratedA = L.collectedA;
			break;
		}


#ifdef UG_PARALLEL
		//PRINTPC(L.pA->process_communicator());
		agglomerate(level);
		// todo: recalculate level information after agglomeration
		/*if(L.bLevelHasMergers && m_agglomerateLevel == level
				&& GetLogAssistant().get_output_process() != -1)
		{
			GetLogAssistant().set_output_process(L.agglomeratedPC.get_proc_id(0));
			UG_LOG("New output proc is " << L.agglomeratedPC.get_proc_id(0) << "\n");
		}*/
			//calculate_level_information(level, createAMGlevelTiming);
		precalc_level(level);
		if(m_agglomerateLevel == level)
		{
			AMG_LOG(0, "Processor " << pcl::ProcRank() << " got agglomerated at level " << level << ".\n");
			//if(GetLogAssistant().get_output_process() != -1)
				//GetLogAssistant().set_output_process(pcl::ProcRank()+1);
			break;
		}
#else
		precalc_level(level);
		L.pAgglomeratedA = L.pA;
#endif

		if(m_writeMatrices)
		{
			AMGWriteToFile(*L.pAgglomeratedA, level, level,
					GetAMGFilenamePath(*L.pAgglomeratedA, level, "A", level, ".mat", true), m_amghelper);
			//PRINTPC(L.pAgglomeratedA->get_process_communicator());
		}

		//smoothem_R[level].init(*m_A[level]);

		levels.push_back(make_sp(new AMGLevel(level+1)));
		AMGLevel &nextL = *levels[level+1];
		nextL.spYConvCheck = spYConvCheck->clone();

		//UG_LOG(ToString(nextL.spLayouts->proc_comm()) << "\n" << ToString(L.spAgglomeratedLayouts->proc_comm()));

		init_smoother(level);

		nextL.pA = make_sp(new matrix_operator_type);
		matrix_operator_type &AH = *nextL.pA;
		prolongation_matrix_type &R = L.R;
		prolongation_matrix_type &P = L.P;

#ifdef UG_PARALLEL
		nextL.spLayouts->proc_comm() = L.spAgglomeratedLayouts->proc_comm();
		nextL.spLayouts->comm() = L.spAgglomeratedLayouts->comm();
		AH.set_storage_type(PST_ADDITIVE);
		R.set_storage_type(PST_ADDITIVE);
		P.set_storage_type(PST_CONSISTENT);

		AH.set_layouts(nextL.spLayouts);
#else
		int tt;	nextL.spLayouts = &tt;
#endif

		Stopwatch SW; SW.start();

		AMG_LOG(1, "===== Creating level " << level << ". (" << L.pAgglomeratedA->num_rows() << " nodes =====\n");

		create_AMG_level(AH, R, *L.pAgglomeratedA, P, level, *nextL.spLayouts);

		SW.stop();
		createAMGlevelTiming = SW.ms();

		// create vectors for AMG multigrid
		/////////////////////////////////////////

		L.corr.resize(L.pAgglomeratedA->num_rows());
		L.cH.resize(AH.num_rows());
		L.dH.resize(AH.num_rows());

	#ifdef UG_PARALLEL
	#ifdef UG_DEBUG
		ConsistencyCheck(levels[level]->is_fine, L.pAgglomeratedA->layouts()->comm(),
				L.pAgglomeratedA->layouts()->proc_comm(),
				L.pAgglomeratedA->layouts()->master(), L.pAgglomeratedA->layouts()->slave(), "fine marks");
//			PrintLayoutIfBroken(AH.get_process_communicator(), nextL.com, nextL.masterLayout, nextL.slaveLayout);
	#endif

		AH.set_layouts(nextL.spLayouts);

		L.corr.set_layouts(L.pAgglomeratedA->layouts());
		L.cH.set_layouts(AH.layouts());
		L.dH.set_layouts(AH.layouts());

		L.corr.set_storage_type(PST_CONSISTENT);
		L.cH.set_storage_type(PST_CONSISTENT);
		L.dH.set_storage_type(PST_ADDITIVE);

		/*UG_LOG("nextLevel Layouts:\n");
				PrintLayout(nextL.pA->process_communicator(), nextL.pA->communicator(), nextL.masterLayout, nextL.slaveLayout);*/
	#endif

		// todo: set size for variable sized blockvectors
		/*for(size_t i=0; i<N; i++)
			if(nodes[i].isCoarse())
			{
				int rows = GetRows(A.begin_row(i).value());
				UG_ASSERT(newIndex[i] >= 0, "");
				SetSize((*m_vec1[level+1])[newIndex[i]], rows);
				SetSize((*m_vec2[level+1])[newIndex[i]], rows);
			}*/
	}

	m_usedLevels = level+1;
#ifdef UG_PARALLEL
	const pcl::ProcessCommunicator &pc = mat.layouts()->proc_comm(); //!
	m_totalUsedLevels = pc.allreduce(m_usedLevels, PCL_RO_MAX);
	AMG_LOG(1, "total used levels: " << m_totalUsedLevels << "\n");
	LevelInformation::parallel_update(pc, m_levelInformation, m_totalUsedLevels);
	for(level=m_usedLevels;level<m_totalUsedLevels; level++)
	{
		LevelInformation &li = m_levelInformation[level];
		AMG_LOG(1, li.tostring() << "\n");
		AMG_LOG(1, " level took " << li.m_dCreationTimeMS << " ms" << std::endl << std::endl);
	}
#else
	m_totalUsedLevels = m_usedLevels;
#endif

	AMG_LOGN(1, "finished. now init fsmoothing.");
	init_fsmoothing();


	AMG_LOG(1, "AMG Setup finished. Total used Levels: " << m_totalUsedLevels << ". ");
	m_dTimingWholeSetupMS = SWwhole.ms();
	UG_DLOGN(LIB_ALG_AMG, 1, "AMG Setup took " << m_dTimingWholeSetupMS << " ms.");

	// calc complexities
	double totalNNZs=0;
	double totalNrOfNodes=0;
	for(size_t i=0; i<m_totalUsedLevels; i++)
	{
		LevelInformation &li = m_levelInformation[i];
		totalNNZs += li.get_nnz();
		totalNrOfNodes += li.get_nr_of_nodes();
	}

	m_dOperatorComplexity = totalNNZs / levels[0]->pA->total_num_connections();
	m_dGridComplexity = totalNrOfNodes / levels[0]->pA->num_rows();

	AMG_LOG(1, "Operator Complexity: " << m_dOperatorComplexity << " grid complexity: "
			<< m_dGridComplexity << std::endl << std::endl);

	m_bInited = true;

	m_vLvlConvRates.clear();
	m_vLvlConvRates.resize(m_totalUsedLevels);

	}UG_CATCH_THROW("AMGBase::" << __FUNCTION__ << " failed")
	catch(...)
	{
		UG_THROW("unknown error");
	}
	return true;
}




template<typename TAlgebra>
void AMGBase<TAlgebra>::create_direct_solver(size_t level)
{
	try{
	AMG_PROFILE_BEGIN(amg_createDirectSolver);

	UG_ASSERT(block_traits< typename vector_type::value_type >::is_static, "dynamic not yet implemented");

	AMGLevel &L = *levels[level];
	matrix_operator_type &A = *L.pA; (void) A;
	size_t static_nrUnknowns = block_traits< typename vector_type::value_type >::static_size;
	AMG_LOG(2, "Creating direct solver on level " << level << " (" << A.num_rows() << " nodes, total "
			<< A.num_rows()*static_nrUnknowns << " unknowns)" << std::endl << "Using Direct Solver on Matrix "
			<< A.num_rows()*static_nrUnknowns << "x" << A.num_rows()*static_nrUnknowns << ". ");
	(void) static_nrUnknowns;
	Stopwatch SW; SW.start();

#ifdef UG_PARALLEL
	const pcl::ProcessCommunicator &pc = A.layouts()->proc_comm(); //!
	if(m_basesolver->supports_parallel() == false && pc.size()>1)
	{
		L.bHasBeenMerged = true;
		L.bLevelHasMergers=true;
		// create basesolver
		CollectMatrixOnOneProc(A, *L.collectedA, L.agglomerateMasterLayout, agglomerateSlaveLayout);
		if(pcl::ProcRank() == pc.get_proc_id(0))
		{
			if(m_writeMatrices)
			{
				m_amghelper.receive_agglomerate_positions(level, A.layouts()->comm(), L.agglomerateMasterLayout, L.collectedA->num_rows());

				std::string name = get_output_directory(level) + "collectedA.mat";
				std::vector<MathVector<3> > vec = m_amghelper.positions[level];
				ConnectionViewer::WriteMatrix(name.c_str(), *L.collectedA, &vec[0], m_amghelper.dimension);
			}

			L.spAgglomeratedLayouts->proc_comm() = pc.create_sub_communicator(true);
			// setting empty layouts to collectedA
			L.collectedA->set_layouts(SmartPtr<AlgebraLayouts>(new AlgebraLayouts));
			m_basesolver->init(L.collectedA);

			L.collC.resize(L.collectedA->num_rows());
			L.collD.resize(L.collectedA->num_rows());
			SetParallelVectorAsMatrix(L.collC, *L.collectedA, PST_ADDITIVE);
			SetParallelVectorAsMatrix(L.collD, *L.collectedA, PST_ADDITIVE);
		}
		else
		{
			if(m_writeMatrices)
				m_amghelper.send_agglomerate_positions(level, A.layouts()->comm(), agglomerateSlaveLayout);

			m_agglomerateLevel = level;
			L.spAgglomeratedLayouts->proc_comm() = pc.create_sub_communicator(false);
		}
	}
	else
	{
		L.bLevelHasMergers=false;
		m_basesolver->init(L.pA);
	}
#else
	m_basesolver->init(L.pA);
#endif

	m_dTimingCoarseSolverSetupMS = SW.ms();
	UG_DLOG(LIB_ALG_AMG, 1, "Coarse Solver Setup took " << m_dTimingCoarseSolverSetupMS << "ms." << std::endl);

	}UG_CATCH_THROW("AMGBase::" << __FUNCTION__ << " failed");
}

//!
//! amg constructor
template<typename TAlgebra>
AMGBase<TAlgebra>::AMGBase() :
	m_numPreSmooth(3),
	m_numPostSmooth(3),

	m_maxLevels(100),
	m_usedLevels(0),

	m_pPositionProvider2d(NULL),
	m_pPositionProvider3d(NULL)
{

	m_presmoother = m_postsmoother = make_sp(new SymmetricGaussSeidel<TAlgebra>);
	m_basesolver = make_sp(new LU<TAlgebra>);
	set_Y_cycle(4, 1e-30, 0.1);

	m_vec4 = NULL;
	m_bInited = false;

	m_maxNodesForBase = 1000;
	m_minNodesOnOneProcessor = 100;
	m_preferredNodesOnOneProcessor = 10000;
	m_dMaxFillBeforeBase = 0.5;

	m_writeMatrices = false;
	m_bFSmoothing = true;

	m_dbgDimension = 0;
	iteration_glboal=0;
	m_checkLevelPostIterations = 10;
	m_iNrOfPreiterationsCheck = 8;
	m_dPreiterationsMimumDefect = 1e-8;
	m_bOneInit = false;
	m_verbose = 0;

	m_dGalerkinTruncation = 1e-12;
	m_bDebugWriteAllVectors = false;
}


template<typename TAlgebra>
void AMGBase<TAlgebra>::cleanup()
{
	m_usedLevels = 0;
	m_bInited=false;

	if(m_vec4 != NULL) delete m_vec4; m_vec4 = NULL;
}
//!
//! amg destructor
template<typename TAlgebra>
AMGBase<TAlgebra>::~AMGBase()
{
	cleanup();
}

template<typename TAlgebra>
void AMGBase<TAlgebra>::init_fsmoothing()
{
	AMG_PROFILE_FUNC();
	// last level is basesolver level
	for(size_t k=0; k<levels.size()-1; k++)
	{
		//UG_LOG("---- level " << k << " ---- ");
		matrix_operator_type &mat = *levels[k]->pAgglomeratedA;


#ifdef UG_PARALLEL
		ParallelVector<Vector< typename matrix_operator_type::value_type > > m_diag;
//		m_diag.set_layouts(mat.layouts()->master(), mat.layouts()->slave());
//		m_diag.set_communicator(mat.communicator());
		m_diag.set_layouts(mat.layouts());
#else
		std::vector<typename matrix_operator_type::value_type > m_diag;
#endif
		size_t size = mat.num_rows();
		levels[k]->m_diagInv.resize(size);
		m_diag.resize(size);


		// copy diagonal
		for(size_t i = 0; i < m_diag.size(); ++i)
			m_diag[i] = mat(i, i);

#ifdef UG_PARALLEL
		//	make diagonal consistent
		m_diag.set_storage_type(PST_ADDITIVE);
		m_diag.change_storage_type(PST_CONSISTENT);
#endif

		for(size_t i = 0; i < m_diag.size(); ++i)
			GetInverse(levels[k]->m_diagInv[i], m_diag[i]);
	}
}


template<typename TAlgebra>
bool AMGBase<TAlgebra>::f_smoothing(vector_type &corr, vector_type &d, size_t level)
{
	AMG_PROFILE_FUNC();
	AMGLevel &L = *levels[level];
	stdvector<bool> &is_fine = L.is_fine;
	matrix_operator_type &A = *L.pAgglomeratedA;


	//UG_ASSERT(L.m_diagInv.size() == is_fine.size(), "fine markers do not match in size on level " << level);
	for(size_t i=0; i<L.m_diagInv.size(); i++)
	{
		if(is_fine[i])
			MatMult(corr[i], 1.0, L.m_diagInv[i], d[i]);
		else
			corr[i] = 0.0;
	}

#ifdef UG_PARALLEL
	corr.set_storage_type(PST_ADDITIVE);
	corr.change_storage_type(PST_CONSISTENT);
#endif

	A.matmul_minus(d, corr);
	return true;
}

template<typename TAlgebra>
bool AMGBase<TAlgebra>::solve_on_base(vector_type &c, vector_type &d, size_t level)
{
	AMG_PROFILE_FUNC();
	m_basesolver->apply_return_defect(c, d);
	return true;
}

template<typename TAlgebra>
void AMGBase<TAlgebra>::calc_level_conv_rates(SmartPtr<IConvergenceCheck<vector_type> > pCheck, size_t level)
{
	if(is_fixed_conv_check()) return;

	LevelConvergenceRates &s = m_vLvlConvRates[level];

	size_t steps = pCheck->step();
	s.avgReduction = pow(s.avgReduction, s.totalSteps/((double)s.totalSteps+steps))
		* pow(pCheck->avg_rate(), 1.0/((double)s.totalSteps+steps));

	s.totalSteps += steps;
	s.calls ++;
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// add_correction_and_update_defect:
//------------------------------------

template<typename TAlgebra>
bool AMGBase<TAlgebra>::add_correction_and_update_defect2(vector_type &c, vector_type &d,
		matrix_operator_type &A, size_t level)
{
	THROW_IF_NOT_EQUAL_3(c.size(), d.size(), A.num_rows());

	AMG_PROFILE_FUNC();
	//PRINTPC(A.get_process_communicator());
	// c = consistent, d = additiv

	if(level == m_usedLevels-1)
	{
		solve_on_base(c, d, level);
		return true;
	}

	AMGLevel &L = *levels[level];


#ifdef UG_PARALLEL
	if(!d.has_storage_type(PST_ADDITIVE) || !c.has_storage_type(PST_CONSISTENT))
	{
		UG_LOG("ERROR: In 'amg::check':Inadequate storage format of Vectors.\n");
		return false;
	}
#endif
	vector_type &corr = levels[level]->corr;
	corr.set(0.0);

#ifdef UG_PARALLEL
	corr.set_storage_type(PST_CONSISTENT);
	c.set_storage_type(PST_CONSISTENT);
#endif


	if(m_bDebugWriteAllVectors)
	{
		writevec(std::string("AMG_") + ToString(iteration_glboal++) + "aa_d_L", d, level);
		writevec(std::string("AMG_") + ToString(iteration_glboal++) + "aa_c_L", c, level);
	}

	// presmooth
	for(size_t i=0; i < m_numPreSmooth; i++)
	{
		corr.set(0.0);
		L.presmoother->apply_update_defect(corr, d);
        c += corr;
	}
#ifdef UG_PARALLEL
	UG_ASSERT(c.has_storage_type(PST_CONSISTENT), "" );
#endif
	// pre f-smoothing
	if(m_bFSmoothing)
	{
		corr.set(0.0);
		f_smoothing(corr, d, level);
		c+=corr;
	}
#ifdef UG_PARALLEL
	UG_ASSERT(c.has_storage_type(PST_CONSISTENT), "" );
#endif

	if(m_bDebugWriteAllVectors)
	{
		writevec(std::string("AMG_") + ToString(iteration_glboal++) + "pa_d_L", d, level);
		writevec(std::string("AMG_") + ToString(iteration_glboal++) + "pa_c_L", c, level);
	}

	vector_type &cH = levels[level]->cH;
	vector_type &dH = levels[level]->dH;
	cH.set(0.0);

#ifdef UG_PARALLEL
	cH.set_storage_type(PST_CONSISTENT);
#endif

	dH.set(0.0);
	// restrict defect
	// R = additive, d additive -> consistent. dH is then additive.
	// dH = R*d;
	L.R.apply(dH, d);
#ifdef UG_PARALLEL
	UG_ASSERT(dH.has_storage_type(PST_ADDITIVE), "" );
#endif


	if(m_bDebugWriteAllVectors)
		writevec(std::string("AMG_") + ToString(iteration_glboal++) + "_dH_L", dH, level+1);

	// apply lmgc on coarser nodes

	if(level+1 == m_usedLevels-1)
		add_correction_and_update_defect(cH, dH, level+1);
	else
	{
		// all V/W/Y-Cycles are handled here
		// V-Cycle: spYConvCheck = FixedConvergenceCheck(1)
		// W-Cycle: spYConvCheck = FixedConvergenceCheck(2)
		// Y-Cycle: spYConvCheck = StdConvergenceCheck(20, 1e-24, 0.01)
		//SmartPtr<IConvergenceCheck<vector_type> > pCheck = levels[level]->spYConvCheck;
	  SmartPtr<IConvergenceCheck<vector_type> > pCheck(new FixedConvergenceCheck<vector_type>(1));
		//FixedConvergenceCheck<vector_type> pCheck(m_iYCycle);
		pCheck->start(dH);

		while(!pCheck->iteration_ended())
		{
			add_correction_and_update_defect(cH, dH, level+1);
			pCheck->update(dH);
		}

		calc_level_conv_rates(pCheck, level+1);
	}

	corr.set(0.0);
	//cH.set(0.0);
	// interpolate correction. P = consistent, cH = consistent. -> corr consistent
	// corr = R*cH
	L.P.apply(corr, cH);
#ifdef UG_PARALLEL
	UG_ASSERT(corr.has_storage_type(PST_CONSISTENT), "" );
#endif

	// add coarse grid correction to level correction
	// c += corr;

	c += corr;

		// update defect
	// d = d - Ah*corr
	A.matmul_minus(d, corr);


	if(m_bDebugWriteAllVectors)
	{
		writevec(std::string("AMG_") + ToString(iteration_glboal++) + "cc_d_L", d, level);
		writevec(std::string("AMG_") + ToString(iteration_glboal++) + "cc_c_L", c, level);
	}


	// post f-smoothing
	if(m_bFSmoothing)
	{
		corr.set(0.0);
		f_smoothing(corr, d, level);
		c+=corr;
	}

	// postsmooth
	for(size_t i=0; i < m_numPostSmooth; i++)
	{
		corr.set(0.0);
		L.postsmoother->apply_update_defect(corr, d);
		c += corr;
	}


	if(m_bDebugWriteAllVectors)
	{
		writevec(std::string("AMG_") + ToString(iteration_glboal++) + "pp_d_L", d, level);
		writevec(std::string("AMG_") + ToString(iteration_glboal++) + "pp_c_L", c, level);
	}

	return true;
}



template<typename TAlgebra>
bool AMGBase<TAlgebra>::get_correction(vector_type &c, const vector_type &const_d)
{
	AMG_PROFILE_FUNC();
	THROW_IF_NOT_EQUAL_3(c.size(), const_d.size(), levels[0]->pA->num_rows());

	m_vLvlConvRates[0].calls ++;
	if(m_vec4 == NULL)
		m_vec4 = new vector_type;

	m_vec4->resize(const_d.size());
#ifdef UG_PARALLEL
	m_vec4->set_layouts(c.layouts());
	m_vec4->set_storage_type(PST_ADDITIVE);
#endif
	vector_type &d = *m_vec4;


	d = const_d;

	/*if(m_usedLevels == 1)
	{
		solve_on_base(c, d, 0);
		return true;
	}
	else*/
	{
		c.set(0.0);
		return add_correction_and_update_defect(c, d);
	}
}



template<typename TAlgebra>
std::string AMGBase<TAlgebra>::config_string() const
{
	std::stringstream ss;
	ss << "AMGBase( " << m_cycleType;
	ss << ", max levels = " << m_maxLevels;
	ss << ", F-Smoothing = " << (m_bFSmoothing ? "[ON]" : "[OFF]") << ")\n";
	ss << " Agglomeration: minNodesOnOneProcessor = " << m_minNodesOnOneProcessor << ", preferredNodesOnOneProcessor = " <<
			m_preferredNodesOnOneProcessor << "\n";

	if(m_presmoother == m_postsmoother)
		ss << " Smoother (" << m_numPreSmooth << "x pre, " << m_numPostSmooth << "x post): "
			<< ConfigShift(m_presmoother->config_string());
	else
	{
		ss << " Presmoother (" << m_numPreSmooth << "x): " << ConfigShift(m_presmoother->config_string());
		ss << " Postsmoother ( " << m_numPostSmooth << "x): " << ConfigShift(m_postsmoother->config_string());
	}

	ss << "\n";
	ss << " Basesolver ( maxNodesForBase = " << m_maxNodesForBase << ", maxFillInBeforeBase = " << m_dMaxFillBeforeBase*100 << "%):";
	ss << ConfigShift(m_basesolver->config_string()) << "\n";

	ss << " Galerkin Truncation: " << m_dGalerkinTruncation << "\n";

	if(m_writeMatrices)
		ss << " Writing Matrices to \"" << m_writeMatrixPath << "\"\n";

	if(m_bInited)
	{
		ss << " AMG is initialized.\n";
		ss << " Used Levels: " << m_usedLevels << "\n";
	}
	return ss.str();
}


template<typename TAlgebra>
void AMGBase<TAlgebra>::update_positions()
{
	AMG_PROFILE_FUNC();
	UG_COND_THROW(m_pPositionProvider2d.valid() && m_pPositionProvider3d.valid(), "specify EITHER positions 2d or 3d");
	if(m_pPositionProvider2d == SPNULL && m_pPositionProvider3d == SPNULL)
		return;
	if(m_pPositionProvider2d.valid())
	{
		std::vector<MathVector<2> > pos;
		m_pPositionProvider2d->get_positions(pos);
		set_positions(&pos[0], pos.size());
	}
	if(m_pPositionProvider3d.valid())
	{
		std::vector<MathVector<3> > pos;
		m_pPositionProvider3d->get_positions(pos);
		set_positions(&pos[0], pos.size());
	}
}



template<typename TAlgebra>
std::string AMGBase<TAlgebra>::GetAMGFilenamePath(const matrix_type &m, int levelDir, std::string name, int level, std::string extension, bool bWriteHeader)
{
#ifdef UG_PARALLEL
	return ConnectionViewer::GetParallelName(get_output_directory(levelDir) + name + "_L" + ToString(level) + extension, m.layouts()->proc_comm(), bWriteHeader);
#else
	return get_output_directory(levelDir) + name + "_L" + ToString(level) + extension;
#endif
}

template<typename TAlgebra>
std::string AMGBase<TAlgebra>::GetAMGFilename(const matrix_type &m, std::string name, int level, std::string extension, bool bWriteHeader)
{
#ifdef UG_PARALLEL
	return ConnectionViewer::GetParallelName(name + "_L" + ToString(level) + extension, m.layouts()->proc_comm(), bWriteHeader);
#else
	return name + "_L" + ToString(level) + extension;
#endif
}

template<typename TAlgebra>
template<typename TNodeType>
void AMGBase<TAlgebra>::write_debug_matrix_markers
	(size_t level, const TNodeType &nodes)
{
	try{
	matrix_type &A = *levels[level]->pA;
	// todo: replace some day with something like nodes.get_mark_count(), nodes.mark_name(i), nodes.is_marked(i)


	std::string path = get_output_directory(level);

	AMG_PROFILE_FUNC();
	std::fstream ffine((path + GetAMGFilename(A, "_fine", level, ".marks")).c_str(), std::ios::out);
	ffine << "1 0 0 1 0\n";
	std::fstream ffine2((path + GetAMGFilename(A, "_aggfine", level, ".marks")).c_str(), std::ios::out);
	ffine2 << "1 0.2 1 1 0\n";
	std::fstream fcoarse((path + GetAMGFilename(A, "_coarse", level, ".marks")).c_str(), std::ios::out);
	fcoarse << "0 0 1 1 2\n";
	std::fstream fother((path + GetAMGFilename(A, "_other", level, ".marks")).c_str(), std::ios::out);
	fother << "1 1 0 1 0\n";
	std::fstream fdirichlet((path + GetAMGFilename(A, "_dirichlet", level, ".marks")).c_str(), std::ios::out);
	fdirichlet << "0 1 1 1 0\n";
	for(size_t i=0; i < nodes.size(); i++)
	{
		//int o = m_amghelper.GetOriginalIndex(level, i);
		int o=i;
		if(nodes[i].is_aggressive_fine()) ffine2 << o << "\n";
		else if(nodes[i].is_fine_direct()) ffine << o << "\n";
		else if(nodes[i].is_coarse()) fcoarse << o << "\n";
		else if(nodes[i].is_dirichlet()) fdirichlet << o << "\n";
		else fother << o << "\n";
	}
	}UG_CATCH_THROW("AMGBase::write_debug_matrix_markers failed, level = " << level);
}




template<typename TAlgebra>
template<typename TMatrix>
void AMGBase<TAlgebra>::
	write_debug_matrix(TMatrix &mat, size_t fromlevel, size_t tolevel, const char *name)
{
	try{
	if(m_writeMatrices)
	{
		AMG_PROFILE_FUNC();

		int level = std::min(fromlevel, tolevel);

		std::string filename = GetAMGFilenamePath(mat, level, name, fromlevel, ".mat", true);
		AMGWriteToFile(mat, fromlevel, tolevel, filename, m_amghelper);
		ConnectionViewer::AddMarkers(filename, GetAMGFilename(mat, "_fine", level, ".marks"));
		ConnectionViewer::AddMarkers(filename, GetAMGFilename(mat, "_coarse", level, ".marks"));
		ConnectionViewer::AddMarkers(filename, GetAMGFilename(mat, "_other", level, ".marks"));
		ConnectionViewer::AddMarkers(filename, GetAMGFilename(mat, "_dirichlet", level, ".marks"));
	}
	}UG_CATCH_THROW("AMGBase::write_debug_matrix failed, fromlevel = " << fromlevel << " tolevel = " << tolevel <<
			" name = " << name);
}

/*template<typename TAlgebra>
void AMGBase<TAlgebra>::save(size_t level, const vector_type &u, const vector_t &b)
{
	AMGLevel &L = *levels[level];
	matrix_operator_type &A = *L.pA;

	std::fstream f;
	SerializeParallelData(f, A);
	SerializeUniquePart(f, A);
	SerializeUniquePart(f, u);
	SerializeUniquePart(f, b);
	bool b = m_amghelper.has_positions();
	if(b)
		Serialize(f, m_amghelper.positions[level]]);

}

template<typename TAlgebra>
void AMGBase<TAlgebra>::load(size_t level, const vector_type &u, const vector_t &b)
{
	AMGLevel &L = *levels[level];
	matrix_operator_type &A = *L.pA;

	std::fstream f;
	SerializeParallelData(f, A);
	SerializeUniquePart(f, A);
	SerializeUniquePart(f, u);
	SerializeUniquePart(f, b);
	bool b = m_amghelper.has_positions();
	if(b)
		Serialize(f, m_amghelper.positions[level]]);

}
*/


template<typename TAlgebra>
void AMGBase<TAlgebra>::write_debug_matrices(matrix_type &AH, prolongation_matrix_type &R, const matrix_type &A,
		prolongation_matrix_type &P, size_t level)
{
	try{
	AMG_PROFILE_FUNC();
	AMG_LOG(1, "writing 4 matrices to " << m_writeMatrixPath << ": A ");

	//if(A.num_rows() < AMG_WRITE_MATRICES_MAX)
	{
		write_debug_matrix(A, level, level, "A");		AMG_LOG(1, "P ");
		write_debug_matrix(P, level+1, level, "P");	AMG_LOG(1, "R ");
		write_debug_matrix(R, level, level+1, "R");	AMG_LOG(1, "AH ");
	}

	AMGWriteToFile(AH, level+1, level+1,
			GetAMGFilenamePath(AH, level, "AH", level+1, ".mat", true), m_amghelper);

	AMG_LOG(1, " done.\n");
	}UG_CATCH_THROW("AMGBase::write_debug_matrices failed, level = " << level);
}

template<typename TAlgebra>
bool AMGBase<TAlgebra>::injection(vector_type &vH, const vector_type &v, size_t level)
{
	size_t nH = levels[level]->P.num_cols(); (void) nH;
	size_t n = levels[level]->P.num_rows(); (void) n;
	UG_COND_THROW(nH != vH.size() || n > v.size(), nH << " != " << vH.size() << " or " << n << " > " << v.size());

	for(size_t j=0; j < nH; j++)
		vH[j] = v[m_parentIndex[level+1][j]];
	return true;
}


template<typename TAlgebra>
std::string
AMGBase<TAlgebra>::get_output_directory(int level)
{
	std::string path = m_writeMatrixPath << "/level" << level << "/";
	CreateDirectory(path.c_str(), 0777);
	return path;
}



template<typename TAlgebra>
double
AMGBase<TAlgebra>::Dnorm(vector_type &v, size_t level)
{
	return DNorm(*levels[level]->pA, v);
}

template<typename TAlgebra>
double
AMGBase<TAlgebra>::Anorm(vector_type &v, size_t level)
{
	return ANorm(*levels[level]->pA, v);
}

/*template<typename TAlgebra>
void
AMGBase<TAlgebra>::calculate_smoothing_property(vector_type &t)
{
	// <ASt, St> <= <Ae, e> - \sigma <D^{-1}Ae, Ae}
	// dh. sigma = min_e (<ASt, St> - <Ae, e>) / <D^{-1}Ae, Ae}
}*/

// end group plugins_amg
/// \}

} // namespace ug

#include "amg_agglomeration.h"
#include "amg_checks.h"
#include "amg_communicate_prolongation.h"
#endif //  __H__UG__LIB_ALGEBRA__AMG_SOLVER__AMG_BASE_IMPL_H__

