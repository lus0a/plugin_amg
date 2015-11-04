#ifndef __H__LIB_ALGEBRA__AMG__FAMG_IMPL_H__
#define __H__LIB_ALGEBRA__AMG__FAMG_IMPL_H__

//#include "sparsematrix_util.h"

//#include "famg_nodeinfo.h"
#include "common/stopwatch.h"
#include "common/assert.h"
#include "famg_testvectors.h"
#include "lib_algebra/operator/preconditioner/jacobi.h"
//#include "maxheap.h"

namespace ug{

/// \addtogroup plugins_amg_famg
/// \{

template<typename TAlgebra>
FAMG<TAlgebra>::FAMG() : AMGBase<TAlgebra>()
{
	set_or_copy(NULL);
}


template<typename TAlgebra>
FAMG<TAlgebra>::FAMG(const FAMG<TAlgebra> &parent) : base_type(parent)
{
	set_or_copy(&parent);
}

template<typename TAlgebra>
void FAMG<TAlgebra>::set_or_copy(const FAMG<TAlgebra> *parent)
{
#define SET_OR_COPY(variableName, defaultValue) if(parent!=NULL) variableName=parent->variableName; else variableName=defaultValue;

	SET_OR_COPY(m_theta, 0.9);
	SET_OR_COPY(m_delta, 0.1);
	SET_OR_COPY(mDeltaAll, 0.1);

	SET_OR_COPY(m_bHAsInverse, false);
	SET_OR_COPY(m_bAggressiveCoarsening, false);
	SET_OR_COPY(m_writeTestvectors, false);
	SET_OR_COPY(m_bWriteFValues, false);

	SET_OR_COPY(m_dProlongationTruncation, 0.0);

	SET_OR_COPY(m_dPrereduceAToStrongParameter, 0.0);

	SET_OR_COPY(m_dHReduceInterpolationNodesParameter, 0.0);
	SET_OR_COPY(m_iHReduceMinParents, 1);

	SET_OR_COPY(m_dStrongConnectionExternal, 0.1);
	SET_OR_COPY(m_iPreventFFConnection, 0);

	SET_OR_COPY(m_iStrongConnectionSubMask, -1);

	SET_OR_COPY(m_testvectorCalculator, make_sp(new TestvectorsLinearComponent<TAlgebra>()));

	SET_OR_COPY(m_bExternalCoarsening, false);
	SET_OR_COPY(m_bUsePrecalculate, true);

	SET_OR_COPY(m_bWriteFValues, false);
	SET_OR_COPY(m_bProjectedEVP, false);
	SET_OR_COPY(m_bInjectionRestriction, false);
	SET_OR_COPY(m_fileCoarsening, NONE);
	SET_OR_COPY(m_iMaxParents, 2);
	SET_OR_COPY(m_threeParentsPenalty, 1.1);

	SET_OR_COPY(m_bUseConstantIfFailed, false);
	SET_OR_COPY(m_dDeficitSolveTolerance, 1e-7);
	SET_OR_COPY(m_bLE_XC_Griebel, false);
	SET_OR_COPY(m_bXC_new, false);

	SET_OR_COPY(m_bXCaggressive, false);
	SET_OR_COPY(m_bIndirectInterpolation, true);

	SET_OR_COPY(m_iTestvectorDamps, 2);
	SET_OR_COPY(m_dDampingForSmootherInInterpolationCalculation, 0.8);
	SET_OR_COPY(m_testvectorsmoother, make_sp(new Jacobi<TAlgebra>()));

	SET_OR_COPY(m_iXCCalcInterpolationRoutine, 3);
	SET_OR_COPY(m_bBeam, false);

	SET_OR_COPY(m_bUseN2Interpolation, false);
	SET_OR_COPY(m_bFullBlockJacobi, true);

#ifdef UG_PARALLEL
	SET_OR_COPY(m_pParallelCoarsening, GetColorCoarsening());
#endif

#undef SET_OR_COPY
}

template<typename TAlgebra>
FAMG<TAlgebra>::~FAMG()
{

}

template<typename TAlgebra>
std::string FAMG<TAlgebra>::config_string() const
{
	std::stringstream ss;
	ss << "FAMG Preconditioner:\n";

	// coarsening

	ss << "  Injection Restriction = " << OnOffString(m_bInjectionRestriction)
			<< ", submask = " << get_strong_connection_submask_desc()
			<< ",  m_bIndirectInterpolation = " << OnOffString(m_bIndirectInterpolation)
			<< "\n";


	ss << " TRUNCATION PARAMETERS:\n";
	ss << "  Prereduce A parameter: " << m_dPrereduceAToStrongParameter << "\n";
	ss << "  Truncation of interpolation = " << m_dProlongationTruncation << "\n";
	ss << "  DeficitSolveTolerance = " << m_dDeficitSolveTolerance << "\n";

	ss << " COARSENING:\n";
	ss << "  Aggressive Coarsening is " << OnOffString(m_bAggressiveCoarsening) << "\n";
	if(m_bExternalCoarsening)
	{
		ss << "  External Coarsening: StrongConnection = " << m_dStrongConnectionExternal << ", StrongConnectionSub = " << get_strong_connection_submask_desc() << "\n";
		ss << "  XC Calc Interpolation Routine = " << m_iXCCalcInterpolationRoutine << "\n";
#ifdef UG_PARALLEL
		ss << "  ParallelCoarsening: " << ConfigShift(m_pParallelCoarsening->config_string()) << "\n";
#endif
	}
	else
	{
		if(m_bUsePrecalculate)	ss << "  precalculate Coarsening\n";
		else					ss << "  std FAMG coarsening.\n";
	}

	if(m_bLE_XC_Griebel)
		ss << "  GRIEBEL Linear Elasticity External Neumann-Coarsening.\n";
	if(m_bXC_new)
		ss << "  New External Neumann-Coarsening.\n";

	ss << " INTERPOLATION CALCULATION:\n";
	if(m_bHAsInverse)
		ss << "  USING H AS INVERSE!\n";
	ss << "  Delta: " << mDeltaAll << " (forces interpolation quality measure F < delta.) " << std::endl;
	ss << "  Theta: " << m_theta << " (with multiple parents pairs, discard pairs with m_theta * F > min F.) " << std::endl;
	ss << "  Damping for Smoother in interpolation calculation: " << m_dDampingForSmootherInInterpolationCalculation << std::endl;
	ss << "  max parents: " << m_iMaxParents << ", ";
	ss << "  H reduce: [ interpolation nodes parameter: " << m_dHReduceInterpolationNodesParameter << ", min parents: " << m_iHReduceMinParents << "]\n";
	ss << "  Use Neighboorhood 2 Interpolation if necessary: " << OnOffString(m_bUseN2Interpolation) << "\n";
	ss << "  Use Full Block Jacobi: " << OnOffString(m_bFullBlockJacobi) << "\n";

	ss << " TESTVECTORS:\n";
	if(m_testvectorsmoother == SPNULL)
		ss << "  No testvectorsmoother set!!!;\n";
	else
		ss << "  Testvector Smoother (" << m_iTestvectorDamps << "x): " << m_testvectorsmoother->config_string() << "\n";

	if(m_testvectorCalculator == SPNULL)
			ss << "  No testvector set!!!;\n";
		else
			ss << "  TestvectorCalculator: " << ConfigShift(m_testvectorCalculator->config_string()) << "\n";
	ss << "  UseConstantIfFailed = " << OnOffString(m_bUseConstantIfFailed) << "\n";
	if(m_writeMatrices && m_writeTestvectors)
		ss << "  Write Testvectors is on.\n";




//	ss << m_testvectors.size() + m_vVectorWriters.size() << " test vectors" << std::endl; ?
	ss << " AMGBase: " << ConfigShift(AMGBase<TAlgebra>::config_string()) << "\n";

	return ss.str();
}


template<typename TAlgebra>
void FAMG<TAlgebra>::precalc_level(size_t level)
{
	try{
	AMG_PROFILE_FUNC();
	//UG_LOG("\n\n\nprecalc!\n\n\n");
	UG_COND_THROW(m_testvectorsmoother == SPNULL, "please provide a testvector smoother.");
	UG_COND_THROW(m_testvectorCalculator == SPNULL, "please provide a testvectors.");
	UG_COND_THROW(m_testvectorCalculator->needs_positions() == true && m_amghelper.has_positions() == false,
			"your selected testvector calculator needs position, but they are not there?.");




#ifdef UG_PARALLEL
	// check this: do FE tv need to be known?
	if(levels[level]->bHasBeenMerged)
	{
		vector_type t;
		if(super::isMergingMaster(level))
			t.resize(levels[level]->pAgglomeratedA->num_rows());
		for(size_t i=0; i<testvectors.size(); i++)
		{
			t.set_storage_type(PST_CONSISTENT);
			testvectors[i].set_storage_type(PST_CONSISTENT);
			super::gather_vertical(testvectors[i], t, level, PST_CONSISTENT);
			testvectors[i] = t;
		}

		if(super::isMergingMaster(level) == false)
			return;
	}
#endif
	if(level == 0)
	{
		AMG_PROFILE_BEGIN(testvectorCalculator);
		try{
			testvectors.clear();
			MathVector<3> *pos = m_amghelper.has_positions() ? &(m_amghelper.positions[0][0]) : NULL;
			m_testvectorCalculator->calc(levels[level]->pAgglomeratedA, pos, testvectors);
		}UG_CATCH_THROW("testvectorCalculator " << m_testvectorCalculator->name() << " failed");
	}
	}UG_CATCH_THROW("FAMG::precalc_level failed, level = " << level);
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// createFAMGLevel:
//-------------------------
/**
 * create FAMG matrix R, P, and AH = R A P
 * \param AH
 * \param R
 * \param A
 * \param P
 * \param level
 */
template<typename TAlgebra>
void FAMG<TAlgebra>::create_AMG_level(matrix_type &AH, prolongation_matrix_type &R, const matrix_type &A,
		prolongation_matrix_type &P, size_t level, AMGLayouts &nextLevelLayouts)
{
	try{
	AMG_PROFILE_FUNC();


	if(mDeltaLevel.find(level) != mDeltaLevel.end())
		m_delta = mDeltaLevel[level];
	else
		m_delta = mDeltaAll;

	UG_COND_THROW(m_testvectorsmoother == SPNULL, "please provide a testvector smoother.");

	if(m_bProjectedEVP && level == 0)
	{

	}
	UG_COND_THROW(testvectors.size() == 0, "testvectors?");
	UG_COND_THROW(testvectors[0].size() != A.num_rows(), "testvectorsize = " << testvectors[0].size() << " != A.num_rows() = " << A.num_rows() << " ?");

	if(m_writeMatrices && m_writeTestvectors)
	{
		for(size_t i=0; i<testvectors.size(); i++)
			ConnectionViewer::WriteVector(Stringify() << m_writeMatrixPath << "testvector" << i << "_L" << level << ".vec",
					testvectors[i], &m_amghelper.positions[level][0], m_dbgDimension);
	}

	//UG_ASSERT(testvectors.size() > 0, "we need at least one testvector.");

	// testvectors will be altered by FAMGLevelCalculator

	// UG_SET_DEBUG_LEVEL(LIB_ALG_AMG, 4);

	if(m_dPrereduceAToStrongParameter > 0.0)
	{
		assert(0); // no thanks
		matrix_type Aeps;
		ReduceToStrongConnections(Aeps, A, m_dPrereduceAToStrongParameter);

		FAMGLevelCalculator<algebra_type>
			dummy(*this, AH, R, Aeps, A, P, level, testvectors, omega, nextLevelLayouts);
		dummy.do_calculation();
	}
	else
	{
		FAMGLevelCalculator<algebra_type> dummy(*this, AH, R, A, A, P, level, testvectors, omega, nextLevelLayouts);
		dummy.do_calculation();
	}

	}UG_CATCH_THROW("FAMG::create_AMG_level failed, level = " << level);

}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//       check_testvector
//--------------------------------
/**
 */
template<typename TAlgebra>
bool FAMG<TAlgebra>::check_testvector(size_t fromlevel, size_t tolevel)
{
	try{
	AMG_PROFILE_FUNC();
#ifdef UG_PARALLEL
	if(pcl::NumProcs() > 1)
	{
		UG_LINE_LOG();
		UG_LOG("WARNING: NOT AVAILABLE IN PARALLEL");
		return false;
	}
#endif
	UG_LOG("\n");
	UG_LOG("            check_testvector\n");
	UG_LOG("==========================================\n")
	UG_LOG("check if testvectors are interpolated exactly\n")
	if (0){
		int level = super::m_usedLevels-1;
		matrix_type &A = *super::levels[level]->pA;
		for(size_t i=0; i<A.num_rows(); i++)
		{
			typename matrix_type::value_type sum=0;
			for(typename matrix_type::row_iterator it = A.begin_row(i); it != A.end_row(i); ++it)
				sum += it.value();
			if(BlockNorm(sum) > 0.1)
			{
				UG_LOG("Row " << i << " has sum " << sum << ": ");
				A.pr(i);
				UG_LOG("\n");
			}
		}
	}

	size_t preSmooth = super::get_num_presmooth();
	size_t postSmooth = super::get_num_postsmooth();
	super::set_num_presmooth(get_testvector_smooths());
	super::set_num_postsmooth(get_testvector_smooths());

	for(size_t level=0; ; level++)
	{
		if(level > tolevel) break;
		// get testvectors, get agglomerated tv
		precalc_level(level);
#ifdef UG_PARALLEL
		if(super::isMergingSlave(level))
		{
			UG_LOG("merged.\n");
			break;
		}
#endif
		if(level >=fromlevel)
		{
			UG_LOG("Check Level " << level << "\n-----------------------------------\n")
		}

		SmartPtr<matrix_type> pA, pAH;
		pA = super::levels[level]->pAgglomeratedA;
		pAH = super::levels[level+1]->pAgglomeratedA;

		matrix_type &A = *pA;
		matrix_type &AH = *pAH;
		prolongation_matrix_type &R = super::levels[level]->R;
		prolongation_matrix_type &P = super::levels[level]->P;
		vector_type c, d, tv;

	#ifdef UG_PARALLEL
		SetParallelVectorAsMatrix(c, A, PST_CONSISTENT);
		SetParallelVectorAsMatrix(d, A, PST_CONSISTENT);
		SetParallelVectorAsMatrix(tv, A, PST_CONSISTENT);
	#endif
		tv.resize(A.num_rows());
		d.resize(A.num_rows());
		c.resize(A.num_rows());
		for(size_t i=0; i<d.size(); i++)
			tv[i] = testvectors[0][i];
		A.apply(d, tv);
		c.set(0.0);
		typename super::checkResult res;

		/*PRINTLAYOUT(A.get_process_communicator(), A.get_communicator(), d.get_master_layout(), d.get_slave_layout());
		d.change_storage_type(PST_CONSISTENT);

		UG_LOG("\n\ntv level " << level << "\n");
		tv.print();
		UG_LOG("\n\nA level " << level << "\n");
		A.print();
		UG_LOG("\n\ntestvectors level " << level << "\n");
		testvectors[0].print();
		UG_LOG("\n\nd level " << level << "\n");
		d.print();

		d.change_storage_type(PST_ADDITIVE);*/

		if(level >=fromlevel)
			super::check_level(c, d, A, level, res, &tv);

		AMGLayouts nextLevelLayouts;
		if(level+1 < super::m_usedLevels-1)
		{
			matrix_type Aeps;
			ReduceToStrongConnections(Aeps, A, m_dPrereduceAToStrongParameter);
			FAMGLevelCalculator<algebra_type>
				dummy(*this, AH, R, Aeps, A, P, level, testvectors, omega, nextLevelLayouts);
			dummy.onlyTV();
		}
		else break;

	}
	super::set_num_presmooth(preSmooth);
	super::set_num_postsmooth(postSmooth);

	}UG_CATCH_THROW("FAMG::check_testvector failed");
	return true;
}




template<typename TAlgebra>
void FAMG<TAlgebra>::init_for_EVP(const vector_type &u, const vector_type &b)
{
	CloneVector(evpU, u);
	CloneVector(evpB, b);
	m_bProjectedEVP=true;
}

// end group plugins_amg_famg
/// \}

} // namespace ug

#include "famg_level_calculator.h"
#endif //  __H__LIB_ALGEBRA__AMG__FAMG_IMPL_H__

