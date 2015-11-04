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
#include "famg/debug_ids.h"

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
	AMGLevelType &L = *levels[level];
	matrix_operator_type &A = *L.pAgglomeratedA;

	if(m_levelInformation.size() < level+1)
		m_levelInformation.resize(level+1);
	AMGLevelInformation &li = m_levelInformation[level];


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
	PROFILE_FUNC();
	try{
	if(isMergingSlave(level) || m_baseLevel == level) return;

	AMGLevelType &L = *levels[level];
	L.presmoother = m_presmoother->clone();
	L.presmoother->init(L.pAgglomeratedA);

	if(m_presmoother == m_postsmoother)
		L.postsmoother = L.presmoother;
	else
	{
		L.postsmoother = m_postsmoother->clone();
		L.postsmoother->init(L.pAgglomeratedA);
	}

	} UG_CATCH_THROW(PRETTY_FUNCTION << " level " << level << " failed. ")
}

template<typename TAlgebra>
bool AMGBase<TAlgebra>::init_smoothers()
{
	try{
		if(levels.size() == 0) return false;

		size_t N0 = levels[0]->pA->num_rows();
		PROGRESS_START(prog, N0, "Creating AMG - Init Smoothers");

		//for(size_t level=0; level < my_multigrid_levels(); level++) //my_multigrid_levels(); level++)
		for(size_t level=0; level < levels.size(); level++) //my_multigrid_levels(); level++)
		{
			AMGLevelType &L = *levels[level];
			PROGRESS_UPDATE(prog, N0-L.pA->num_rows());
			init_smoother(level);
		}
		return true;
	}UG_CATCH_THROW_FUNC()

	return false;
}


template<typename TAlgebra>
void AMGBase<TAlgebra>:: preprocess_init_positions(matrix_operator_type& mat)
{
	AMG_PROFILE_FUNC();
	update_positions();
	m_amghelper.parentIndex = &m_parentIndex;

	if(needs_positions() && m_dbgDimension == 0)
	{
		SmartPtr<IDebugWriter<algebra_type> > dbgW = debug_writer();

		if(dbgW.valid())
		{
			const char *errStr =
					"debug size and size of matrix not matching.\n"
					"if you are using GridFunctionDebugWriter, make sure you selected the right base level:\n"
					"   dbgW=GridFunctionDebugWriter(approxSpace)\n"
					"   dbgW:set_grid_level(GridLevel(baseLev))\n"
					"   amg:set_debug(dbgW)\n";

			dbgW->update_positions();
			m_dbgDimension = dbgW->get_dim();
			if(m_dbgDimension == 2)
			{
				const std::vector<MathVector<2> >&v = dbgW->template get_positions<2>();
				if(v.size() != mat.num_rows()) UG_THROW(errStr);
				set_positions(&v[0], v.size());
			}
			else if(m_dbgDimension == 3)
			{
				const std::vector<MathVector<3> >&v = dbgW->template get_positions<3>();
				if(v.size() != mat.num_rows()) UG_THROW(errStr);
				set_positions(&v[0], v.size());
			}
			else m_dbgDimension = 0;
		}

		if(m_dbgDimension == 0)
		{
			m_writeMatrices = false;
			UG_ERR_LOG("ERROR: AMG set to write Matrices, but no position data given. Not writing grid.\n");
			UG_LOG("ERROR: AMG set to write Matrices, but no position data given. Not writing grid.\n");
		}
	}

	// init m_amghelper for grid printing
	if(needs_positions())
	{
		UG_COND_THROW(m_dbgDimension == 0, "AMG method needs positions");
		THROW_IF_NOT_EQUAL(m_dbgPositions.size(), mat.num_rows());
		if(m_dbgPositions.size() > 0)
		{
			m_amghelper.positions.resize(1);
			m_amghelper.positions[0] = m_dbgPositions;
		}
		m_amghelper.dimension = m_dbgDimension;
	}
}

template<typename TAlgebra>
bool AMGBase<TAlgebra>::preprocess_check_requirements()
{
	try{
	AMG_PROFILE_FUNC();
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

	return true;
	}UG_CATCH_THROW_FUNC();

	return false;
}

template<typename TAlgebra>
bool AMGBase<TAlgebra>::preprocess_create_level(size_t level, double &createAMGlevelTiming, int oldOutputProcess)
{
	try{

	AMG_PROFILE_FUNC();

	m_baseLevel = -1;

	AMGLevelType &L = *levels[level];
	L.pAgglomeratedA = L.pA;
	calculate_level_information(level, createAMGlevelTiming);

	// Check for direct solver condition
	bool bCreateDirectSolver=false;
	if(m_levelInformation[level].get_nr_of_nodes() < m_maxNodesForBase)
	{
		UG_LOG("AMG: Create Direct solver for " << m_levelInformation[level].get_nr_of_nodes() << " nodes (< m_maxNodesForBase = " << m_maxNodesForBase << "\n");
		bCreateDirectSolver=true;
	}
	else if(m_levelInformation[level].get_fill_in() > m_dMaxFillBeforeBase)
	{
		UG_LOG("AMG: Create Direct solver because fill in = " << m_levelInformation[level].get_fill_in() << " > m_dMaxFillBeforeBase = " << m_dMaxFillBeforeBase << "\n");
		bCreateDirectSolver=true;
	}
	else if(level >= m_maxLevels-1)
	{
		UG_LOG("AMG: Create Direct solver because level = " << level+1 << " >= m_maxLevels = " << m_maxLevels << "\n");
			bCreateDirectSolver=true;
	}
	else if(level > 0 && (m_levelInformation[level].get_nr_of_nodes() == m_levelInformation[level-1].get_nr_of_nodes()) )
	{
		UG_LOG("AMG: Create Direct solver for because no coarsening happened.\n");
		bCreateDirectSolver=true;
	}

	if(bCreateDirectSolver)
	{
		create_direct_solver(level);
		return false;
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

	// some AMGs might need data from previous, now agglomerated levels
	precalc_level(level);



	// set output process to a process which was not agglomerated.
	if(oldOutputProcess == 0 && GetLogAssistant().get_output_process() != -1)
	{
		if(m_agglomerateLevel == level)
		{
			// we got agglomerated. disable output (just set to not-pcl:ProcRank())
			GetLogAssistant().set_output_process(pcl::ProcRank()+1); // just to someone else not me
		}
		else
		{
			// we are still computing. we set the output proc to
			// if we are the 0-th process in the new interfaces,
			// set us as the new output proc
			const pcl::ProcessCommunicator &pc = L.pAgglomeratedA->layouts()->proc_comm();
			if(pc.get_proc_id(0) == pcl::ProcRank())
			{
				GetLogAssistant().set_output_process(pcl::ProcRank());
				UG_LOG("Processor " << pcl::ProcRank() << " is now output process.\n");
			}
		}
	}

	// if agglomerated, end here.
	if(m_agglomerateLevel == level)
	{
		AMG_LOG(0, "Processor " << pcl::ProcRank() << " got agglomerated at level " << level << ".\n");
		return false;
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

	levels.push_back(make_sp(new AMGLevelType(level+1)));
	AMGLevelType &nextL = *levels[level+1];
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
	if(m_bChecks)
	{ ConsistencyCheck(levels[level]->is_fine, L.pAgglomeratedA->layouts()->comm(),
			L.pAgglomeratedA->layouts()->proc_comm(),
			L.pAgglomeratedA->layouts()->master(), L.pAgglomeratedA->layouts()->slave(), "fine marks");
	}
//			PrintLayoutIfBroken(AH.get_process_communicator(), nextL.com, nextL.masterLayout, nextL.slaveLayout);

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
	return true;

	}UG_CATCH_THROW("AMGBase::preprocess_create_level(level = "<<level<<")");
	return false;
}
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////

template<typename TAlgebra>
bool AMGBase<TAlgebra>::preprocess(SmartPtr<MatrixOperator<matrix_type, vector_type> > pOp)
{
	try{
	int oldOutputProcess = GetLogAssistant().get_output_process();
	matrix_operator_type& mat = *pOp;

	AMG_PROFILE_BEGIN(amg_base_preprocess);

	if(m_bOneInit && m_bInited)
		return true;

	m_bEmptyProcess=false;
#ifdef UG_PARALLEL
	if(mat.num_rows() == 0 && mat.layouts()->proc_comm().empty())
	{
		m_bInited=true;
		m_bEmptyProcess=true;
		return true;
	}
#endif


	cleanup();

	preprocess_init_positions(mat);

	if(preprocess_check_requirements()==false) return false;

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
	AMGLevelType *pL = new AMGLevelType(0);
	pL->pA = pOp;
#ifdef UG_PARALLEL
//	pL->com = mat.layouts()->comm();
	pL->spLayouts->comm() = mat.layouts()->comm();
	pL->spLayouts->proc_comm() = mat.layouts()->proc_comm();

	pL->bHasBeenMerged = false;
#endif
	pL->spYConvCheck = spYConvCheck;
	levels.push_back(make_sp(pL));

	PROGRESS_START(prog, pOp->num_rows(), "Creating AMG");
	for(; ; level++)

	{
		AMGLevelType &L = *levels[level];
		PROGRESS_UPDATE(prog, pOp->num_rows()-L.pA->num_rows());

		if(preprocess_create_level(level, createAMGlevelTiming, oldOutputProcess) == false)
			break;
	}

	//GetLogAssistant().set_output_process(oldOutputProcess);

	m_usedLevels = level+1;
#ifdef UG_PARALLEL
	const pcl::ProcessCommunicator &pc = mat.layouts()->proc_comm(); //!
	m_totalUsedLevels = pc.allreduce(m_usedLevels, PCL_RO_MAX);
	AMG_LOG(1, "total used levels: " << m_totalUsedLevels << "\n");
	AMGLevelInformation::parallel_update(pc, m_levelInformation, m_totalUsedLevels);
	for(level=m_usedLevels;level<m_totalUsedLevels; level++)
	{
		AMGLevelInformation &li = m_levelInformation[level];
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
		AMGLevelInformation &li = m_levelInformation[i];
		totalNNZs += li.get_nnz();
		totalNrOfNodes += li.get_nr_of_nodes();
	}

	m_dOperatorComplexity = totalNNZs / m_levelInformation[0].get_nnz();
	m_dGridComplexity = totalNrOfNodes / m_levelInformation[0].get_nr_of_nodes();

	AMG_LOG(1, "Operator Complexity: " << m_dOperatorComplexity << " grid complexity: "
			<< m_dGridComplexity << std::endl << std::endl);

	m_bInited = true;

	m_vLvlConvRates.clear();
	m_vLvlConvRates.resize(m_totalUsedLevels);

	}UG_CATCH_THROW_FUNC()
	catch(...)
	{
		UG_THROW("unknown error");
	}
	return true;
}




template<typename TAlgebra>
void AMGBase<TAlgebra>::set_or_copy(const AMGBase<TAlgebra> *parent)
{
#define SET_OR_COPY(variableName, defaultValue) if(parent!=NULL) variableName=parent->variableName; else variableName=defaultValue;
	m_vec4 = NULL;
	m_bInited = false;
	m_usedLevels = 0;
	m_dbgDimension = 0;
	iteration_glboal=0;
	m_bEmptyProcess=false;

	SET_OR_COPY(m_numPreSmooth, 3);
	SET_OR_COPY(m_numPostSmooth, 3);
	SET_OR_COPY(m_presmoother, make_sp(new SymmetricGaussSeidel<TAlgebra>));
	SET_OR_COPY(m_postsmoother, m_presmoother);
	SET_OR_COPY(m_bFSmoothing, true);
	SET_OR_COPY(m_basesolver, make_sp(new LU<TAlgebra>));

	SET_OR_COPY(m_maxLevels, 100);
	if(parent == NULL)
		set_Y_cycle(4, 1e-30, 0.1);

	// base
	SET_OR_COPY(m_maxNodesForBase, 1000);
	SET_OR_COPY(m_minNodesOnOneProcessor, 100);
	SET_OR_COPY(m_preferredNodesOnOneProcessor, 1000);
	SET_OR_COPY(m_dMaxFillBeforeBase, 0.5);

	// debug
	SET_OR_COPY(m_writeMatrices, false);
	SET_OR_COPY(m_bDebugWriteAllVectors, false);
	SET_OR_COPY(m_verbose, 0);
	SET_OR_COPY(m_pPositionProvider2d, SPNULL);
	SET_OR_COPY(m_pPositionProvider3d, SPNULL);

	SET_OR_COPY(m_bOneInit, false);
	SET_OR_COPY(m_dGalerkinTruncation, 1e-12);

#ifdef UG_DEBUG
	SET_OR_COPY(m_bChecks, true);
#else
	SET_OR_COPY(m_bChecks, false);
#endif

	// check
	SET_OR_COPY(m_checkLevelPostIterations, 8);
	SET_OR_COPY(m_iNrOfPreiterationsCheck, 8);
	SET_OR_COPY(m_dPreiterationsMimumDefect, 1e-8);
#undef SET_OR_COPY
}
/// clone constructor
template<typename TAlgebra>
AMGBase<TAlgebra> :: AMGBase(const AMGBase<TAlgebra> &parent) :
	base_type(parent)
{
	set_or_copy(&parent);
}

//!
//! amg constructor
template<typename TAlgebra>
AMGBase<TAlgebra>::AMGBase()
{
	set_or_copy(NULL);
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
	AMGLevelType &L = *levels[level];
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
void AMGBase<TAlgebra>::calc_level_conv_rates(SmartPtr<IConvergenceCheck<vector_type> > pCheck, size_t level)
{
	if(is_fixed_conv_check()) return;

	AMGLevelConvergenceRates &s = m_vLvlConvRates[level];

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
		levels[level]->corr.set(0.0);
		solve_on_base(levels[level]->corr, d, level);
		c += levels[level]->corr;
		return true;
	}

	AMGLevelType &L = *levels[level];


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
		writevec(Stringify() << "AMG_" << iteration_glboal++ << "aa_d_L", d, level);
		writevec(Stringify() << "AMG_" << iteration_glboal++ << "aa_c_L", c, level);
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
		writevec(Stringify() << "AMG_" << iteration_glboal++ << "pa_d_L", d, level);
		writevec(Stringify() << "AMG_" << iteration_glboal++ << "pa_c_L", c, level);
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
		writevec(Stringify() << "AMG_" << iteration_glboal++ << "_dH_L", dH, level+1);

	// apply lmgc on coarser nodes

	{
		// all V/W/Y-Cycles are handled here
		// V-Cycle: spYConvCheck = FixedConvergenceCheck(1)
		// W-Cycle: spYConvCheck = FixedConvergenceCheck(2)
		// Y-Cycle: spYConvCheck = StdConvergenceCheck(20, 1e-24, 0.01)
		SmartPtr<IConvergenceCheck<vector_type> > pCheck = levels[level]->spYConvCheck;
		pCheck->set_offset(level);
		pCheck->set_name((Stringify() << "Cycle " << level).str());
	  //SmartPtr<IConvergenceCheck<vector_type> > pCheck(new FixedConvergenceCheck<vector_type>(1));
		//FixedConvergenceCheck<vector_type> pCheck(m_iYCycle);
		pCheck->start(dH);

		//UG_LOG("\n" <<repeat(' ', level+1) << "." <<"\n");
		do
		{
			add_correction_and_update_defect(cH, dH, level+1);
			pCheck->update(dH);
		}
		while(!pCheck->iteration_ended());

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
		writevec(Stringify() << "AMG_" << iteration_glboal++ << "cc_d_L", d, level);
		writevec(Stringify() << "AMG_" << iteration_glboal++ << "cc_c_L", c, level);
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
		writevec(Stringify() << "AMG_" << iteration_glboal++ << "pp_d_L", d, level);
		writevec(Stringify() << "AMG_" << iteration_glboal++ << "pp_c_L", c, level);
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
bool AMGBase<TAlgebra>::injection(vector_type &vH, const vector_type &v, size_t level)
{
	size_t nH = levels[level]->P.num_cols(); (void) nH;
	size_t n = levels[level]->P.num_rows(); (void) n;
	THROW_IF_NOT_EQUAL(nH, vH.size());
	UG_COND_THROW(n > v.size(), n << " > " << v.size());

	for(size_t j=0; j < nH; j++)
		vH[j] = v[m_parentIndex[level+1][j]];
	return true;
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
	if(v.size() == levels[level]->pA->num_rows())
		return ANorm(*levels[level]->pA, v);
	else
		return ANorm(*levels[level]->pAgglomeratedA, v);
}



template<typename TAlgebra>
void AMGBase<TAlgebra>::check_vector(SmartPtr<ITestvectorCalculator<TAlgebra> > testvectorCalculator, std::string desc)
{
	try{
	AMG_PROFILE_FUNC();
#ifdef UG_PARALLEL
	if(pcl::NumProcs() > 1)
	{
		UG_LINE_LOG();
		UG_LOG("WARNING: NOT AVAILABLE IN PARALLEL");
		return;
	}
#endif

	stdvector< vector_type > testvectors;
	for(size_t level=0; level < levels.size()-1; level++)
	{
		if(levels.size() <= 1) return;
		testvectors.clear();
		MathVector<3> *pos = m_amghelper.has_positions() ? m_amghelper.get_positions(level) : NULL;

		testvectorCalculator->calc(levels[level]->pAgglomeratedA, pos, testvectors);


		std::string path=std::string("/level") + ToString(level) + "/";
		CreateDirectory(Stringify() << m_writeMatrixPath << path);

		for(size_t iTV=0; iTV<testvectors.size(); iTV++)
		{
			SmartPtr<matrix_type> pA, pAH;
			pA = levels[level]->pAgglomeratedA;
			pAH = levels[level+1]->pAgglomeratedA;

			matrix_type &A = *pA;
			matrix_type &AH = *pAH;
	//		prolongation_matrix_type &R = levels[level]->R;
			prolongation_matrix_type &P = levels[level]->P;
			vector_type tvH, Ptv;

		#ifdef UG_PARALLEL
			SetParallelVectorAsMatrix(tvH, A, PST_CONSISTENT);
			SetParallelVectorAsMatrix(Ptv, AH, PST_CONSISTENT);
		#endif
			tvH.resize(AH.num_rows());
			Ptv.resize(A.num_rows());


			ConnectionViewer::WriteVector(Stringify() << m_writeMatrixPath << path << "AMG_TV_" << desc << "_" << iTV << "_L" << level << ".vec",
					testvectors[iTV], pos, m_amghelper.dimension);

			injection(tvH, testvectors[iTV], level);

			P.apply(Ptv, tvH);

			ConnectionViewer::WriteVector(Stringify() << m_writeMatrixPath << path << "AMG_TV_PRinj_" << desc << "_" << iTV << "_L" << level << ".vec",
														tvH, pos, m_amghelper.dimension);

			VecScaleAdd(Ptv, 1.0, Ptv, -1.0, testvectors[iTV]);

			ConnectionViewer::WriteVector(Stringify() << m_writeMatrixPath << path << "AMG_TV_1-PRinj_" << desc << "_" << iTV << "_L" << level << ".vec",
								Ptv, pos, m_amghelper.dimension);

			ConnectionViewer::WriteVector(Stringify() << m_writeMatrixPath << path << "AMG_TV_Rinj_" << desc << "_" << iTV << "_L" << level << ".vec",
											tvH, pos, m_amghelper.dimension);


		}
	}

	}UG_CATCH_THROW("FAMG::check_vector failed");
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
#include "amg_basesolver.h"

#endif //  __H__UG__LIB_ALGEBRA__AMG_SOLVER__AMG_BASE_IMPL_H__

