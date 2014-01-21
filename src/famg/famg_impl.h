/**
 * \file plugins/experimental/amg/src/famg/famg_impl.h
 *
 * \author Martin Rupp
 *
 * \date 16.11.2010
 *
 * implementation file for amg
 *
 * Goethe-Center for Scientific Computing 2009-2010.
 */

#ifndef __H__LIB_ALGEBRA__AMG__FAMG_IMPL_H__
#define __H__LIB_ALGEBRA__AMG__FAMG_IMPL_H__

//#include "sparsematrix_util.h"

//#include "famg_nodeinfo.h"
#include "common/stopwatch.h"
#include "common/assert.h"
#include "famg_testvectors.h"
//#include "maxheap.h"

namespace ug{

/// \addtogroup plugins_amg_famg
/// \{

template<typename TAlgebra>
FAMG<TAlgebra>::FAMG() : AMGBase<TAlgebra>()
{
	m_theta = 0.9;
	m_delta = 0.1;

	m_dDampingForSmootherInInterpolationCalculation = 0.8;
	m_bAggressiveCoarsening = false;
	m_writeTestvectors = false;
	m_bWriteFValues = false;

	m_dProlongationTruncation = 0.0;
	
	m_dPrereduceAToStrongParameter = 0.0;

	m_dHReduceInterpolationNodesParameter = 0.0;
	m_iHReduceMinParents = 1;

	m_dStrongConnectionExternal = 0.1;
	m_iPreventFFConnection=0;

	m_iStrongConnectionSubMask=-1;

	m_iTestvectors = TESTVECTORS_DIRICHLET_MATRIX_ROWS;

	m_bExternalCoarsening = false;
	m_bUsePrecalculate = true;

	iDebugLevelOverlapAMG = iDebugLevelOverlapMatrix = iDebugLevelTestvectorCalc = iDebugLevelPhase3 =
			iDebugLevelCalculateParentPairs = iDebugLevelColoring = iDebugLevelRecvCoarsening = iDebugLevelGetRatings =
			iDebugLevelPrecalculateCoarsening =	iDebugLevelAggressiveCoarsening = iDebugLevelSendCoarsening =
			iDebugLevelCommunicateProlongation = iDebugLevelAfterCommunicateProlongation = iDebugLevelExternalCoarsening = 0;

	m_bWriteFValues = false;
	m_bProjectedEVP = false;
	m_bInjectionRestriction = false;
	m_fileCoarsening = NONE;
	m_iMaxParents = 2;
	m_threeParentsPenalty = 1.1;
#ifdef UG_PARALLEL
	m_pParallelCoarsening = NULL;
#endif
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
	ss << " Aggressive Coarsening is " << (m_bAggressiveCoarsening ? "[ON]" : "OFF") << "\n";
	if(m_bExternalCoarsening)
		ss << " External Coarsening: StrongConnection = " << m_dStrongConnectionExternal << ", StrongConnectionSub = " << get_strong_connection_submask_desc() << "\n";
	if(m_bUsePrecalculate)	ss << " precalculate Coarsening\n";
	else					ss << " std FAMG coarsening.\n";

	ss << " Truncation of interpolation = " << m_dProlongationTruncation << "\n";
	ss << " Prereduce A parameter: " << m_dPrereduceAToStrongParameter << "\n";

	ss << " Injection Restriction " << (m_bInjectionRestriction ? "[ON]" : "OFF") << "\n";

	ss << " Delta: " << mDeltaAll << " (forces interpolation quality measure F < delta.) " << std::endl;
	ss << " Theta: " << m_theta << " (with multiple parents pairs, discard pairs with m_theta * F > min F.) " << std::endl;
	ss << " Damping for Smoother in interpolation calculation: " << m_dDampingForSmootherInInterpolationCalculation << std::endl;
	ss << " max parents: " << m_iMaxParents << ", ";
	ss << " H reduce: [ interpolation nodes parameter: " << m_dHReduceInterpolationNodesParameter << ", min parents: " << m_iHReduceMinParents << "]\n";

	if(m_testvectorsmoother.invalid())
		ss << " No testvectorsmoother set!!!;\n";
	else
		ss << " Testvector Smoother (" << m_iTestvectorDamps << "x): " << m_testvectorsmoother->config_string() << "\n";

	ss << " Testvectors: " << testvector_desc() << "\n";
//	ss << m_testvectors.size() + m_vVectorWriters.size() << " test vectors" << std::endl; ?
	if(m_writeMatrices && m_writeTestvectors)
		ss << " Write Testvectors is on.\n";
	ss << " \n";
	ss << " AMGBase:\n" << ConfigShift(AMGBase<TAlgebra>::config_string());

	return ss.str();
}


template<typename TAlgebra>
void FAMG<TAlgebra>::get_testvectors(stdvector<vector_type> &testvectors, stdvector<double> &omega)
{
	AMG_PROFILE_FUNC();

	testvectors.resize(m_vVectorWriters.size() + m_testvectors.size());
	omega.resize(m_vVectorWriters.size() + m_testvectors.size());

	for(size_t i=0; i<m_vVectorWriters.size(); i++)
	{
		size_t index = i;
		m_vVectorWriters[i]->update(testvectors[index]);
		omega[index] = m_omegaVectorWriters[i];

	}
	for(size_t i=0; i<m_testvectors.size(); i++)
	{
		size_t index = i+m_vVectorWriters.size();
		testvectors[index] = m_testvectors[i];
		omega[index] = m_omegaVectors[i];
	}


}

template<typename TAlgebra>
void FAMG<TAlgebra>::get_testvectors_from_matrix_rows
	(const matrix_type &A, stdvector<vector_type> &testvectors, stdvector<double> &omega)
{
	AMG_PROFILE_FUNC();

	testvectors.resize(1);
	omega.resize(1);
	vector_type &v = testvectors[0];
	v.resize(A.num_rows());
	for(size_t i=0; i<A.num_rows(); i++)
	{
		if(A.is_isolated(i))
			v[i] = 0.0;
		else
			v[i] = 1.0;
	}
}

template<typename TAlgebra>
void FAMG<TAlgebra>::precalc_level(size_t level)
{

	AMG_PROFILE_FUNC();
	//UG_LOG("\n\n\nprecalc!\n\n\n");
	UG_ASSERT(m_testvectorsmoother != NULL, "please provide a testvector smoother.");

	if(level == 0)
	{
		testvectors.clear();
		if(m_iTestvectors == TESTVECTORS_DIRICHLET_MATRIX_ROWS)
		{
			UG_LOG("USING TESTVECTORS_DIRICHLET_MATRIX_ROWS\n");
			get_testvectors_from_matrix_rows(*levels[level]->pAgglomeratedA, testvectors, omega);
		}
		else if(m_iTestvectors == TESTVECTORS_LE_2D)
		{
			UG_LOG("USING TESTVECTORS_LE_2D\n");
			UG_ASSERT(m_amghelper.dimension == 2, "dimension needs to be 2");
			GetTestvectorsLE_2d(*levels[level]->pAgglomeratedA, testvectors, omega,
					m_amghelper.positions[0], m_weightC, m_weightZ);
		}
		else if(m_iTestvectors == TESTVECTORS_LE_3D)
		{
			UG_LOG("USING TESTVECTORS_LE_3D\n");
			UG_ASSERT(m_amghelper.dimension == 3, "dimension needs to be 3");
			GetTestvectorsLE_3d(*levels[level]->pAgglomeratedA, testvectors, omega,
					m_amghelper.positions[0], m_weightC, m_weightX, m_weightY, m_weightZ);
		}
		else if(m_iTestvectors == TESTVECTORS_LINEAR)
		{
			UG_LOG("USING TESTVECTORS_LINEAR\n");
			UG_ASSERT(m_amghelper.dimension == 2 || m_amghelper.dimension == 3, "dimension " << m_amghelper.dimension << " not supported");
			if(m_amghelper.dimension == 3)
				GetTestvectorsLinear_3d(*levels[level]->pAgglomeratedA, testvectors, omega,
								m_amghelper.positions[0], m_weightC, m_weightX, m_weightY, m_weightZ);
			else
				GetTestvectorsLinear_2d(*levels[level]->pAgglomeratedA, testvectors, omega,
							m_amghelper.positions[0], m_weightC, m_weightX, m_weightY);


		}
		else
		{
			UG_ASSERT(m_iTestvectors == TESTVECTORS_EXTERNAL, "");
			UG_LOG("USING TESTVECTORS_EXTERNAL\n");
			get_testvectors(testvectors, omega);
		}
	}

#ifdef UG_PARALLEL
	// check this: do FE tv need to be known?
	if((level != 0 || m_iTestvectors == TESTVECTORS_EXTERNAL)
			&& levels[level]->bHasBeenMerged)
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
	}
#endif
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
	AMG_PROFILE_FUNC();


	if(mDeltaLevel.find(level) != mDeltaLevel.end())
		m_delta = mDeltaLevel[level];
	else
		m_delta = mDeltaAll;

	UG_ASSERT(m_testvectorsmoother != NULL, "please provide a testvector smoother.");

	if(m_bProjectedEVP && level == 0)
	{

	}
	UG_ASSERT(testvectors.size() != 0, "testvectors?");
	UG_ASSERT(testvectors[0].size() == A.num_rows(), "testvectorsize = " << testvectors[0].size() << " != A.num_rows() = " << A.num_rows() << " ?");

	if(m_writeMatrices && m_writeTestvectors)
	{
		for(size_t i=0; i<testvectors.size(); i++)
			ConnectionViewer::WriteVector((m_writeMatrixPath + ToString("testvector") + ToString(i) + ToString("_L") + ToString(level) + ".vec").c_str(),
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

}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//       check_testvector
//--------------------------------
/**
 */
template<typename TAlgebra>
bool FAMG<TAlgebra>::check_testvector(size_t fromlevel, size_t tolevel)
{
	AMG_PROFILE_FUNC();

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

