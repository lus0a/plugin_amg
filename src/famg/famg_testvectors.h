/**
 * \file amg_debug.h
 *
 * \author Martin Rupp
 *
 * \date 24.11.10
 *
 * Goethe-Center for Scientific Computing 2010.
 */

#include "common/common.h" // MathVector<3>
#include "famg_level_calculator.h"
#ifdef UG_PARALLEL
#include "lib_algebra/parallelization/parallelization_util.h"
#endif

#include "lib_algebra/operator/eigensolver/additional_math.h"

#ifndef __H__LIB_ALGEBRA__FAMG_SOLVER__FAMG_TESTVECTORS_H__
#define __H__LIB_ALGEBRA__FAMG_SOLVER__FAMG_TESTVECTORS_H__

namespace ug {

/// \addtogroup plugins_amg_famg
/// \{

#if 0
static void get_testvector_constant(stdvector<double> &testvector, stdvector<size_t> &N, double &i_value, size_t i)
{
	for(size_t j=0; j<N.size(); j++)
		testvector[j] = 1.0;
	i_value = 1.0;
}
#endif


/*
template<typename value_type, typename pos_type>
void get_testvector_value_special(value_type &value, const pos_type &myPos, const pos_type &referencePos)
{
	BlockRef(value, TcomponentOut) = myPos[TcomponentIn] - referencePos[TcomponentIn];
}


inline void get_testvector_xx(DenseVector<VariableArray1<double> > &testvector, stdvector<size_t> &N2, int i_index)
{
	return get_testvector_special<0, 0>(testvector, N2, i_index);
}
inline void get_testvector_xy(DenseVector<VariableArray1<double> > &testvector, stdvector<size_t> &N2, int i_index)
{
	return get_testvector_special<0, 1>(testvector, N2, i_index);
}
...
2d :
(x 0), (y 0), (0 x), (0 y)
3d :
(x 0 0) (0 y 0) (0 0 z)  ..?

*/




template<typename algebra_type>
void FAMGLevelCalculator<algebra_type>::smooth_testvectors()
{
	AMG_PROFILE_FUNC();

	// todo: all global?
	Stopwatch SW;
	if(m_bTiming) { UG_LOG("\ncalculating testvector... "); SW.start(); }

	vector_type d; d.resize(A.num_rows());
	vector_type c; c.resize(A.num_rows());
#ifdef UG_PARALLEL
	c.set_layouts(A.layouts());
	d.set_layouts(A.layouts());
#endif

	// output before smoothing of testvector
	if(m_famg.m_writeMatrices && m_famg.m_writeTestvectors)
		for(size_t i=0; i<m_testvectors.size(); i++)
			ConnectionViewer::WriteVectorPar((m_famg.m_writeMatrixPath + ToString("/testvector_pS") + ToString(i) + ToString("_L") + ToString(level) + ".vec").c_str(),
					m_testvectors[i], &m_famg.m_amghelper.positions[level][0], m_famg.m_dbgDimension);

	for(size_t i=0; i<m_testvectors.size(); i++)
	{
		UG_ASSERT(m_testvectors[i].size() == A.num_rows(), m_testvectors[i].size() << " != " << A.num_rows());
#ifdef UG_PARALLEL
		m_testvectors[i].set_storage_type(PST_CONSISTENT);
#endif

		m_famg.m_testvectorsmoother->init(m_famg.levels[level]->pAgglomeratedA);

		for(size_t jj=0; jj < m_famg.get_testvector_smooths(); jj++)
		{
			d.set(0.0);
			c.set(0.0);
			A.apply(d, m_testvectors[i]);
			// don't change dirichlet nodes
			for(size_t i=0; i<A.num_rows(); i++)
				if(A.is_isolated(i)) d[i] = 0.0;
			m_famg.m_testvectorsmoother->apply(c, d);
			for(size_t i=0; i<A.num_rows(); i++)
				if(A.is_isolated(i)) c[i] = 0.0;
			m_testvectors[i] -= c;
		}

		//double e=EnergyProd(m_testvectors[i], A, m_testvectors[i]);;
		//UG_LOG("EnergyProd " << i << " = " << e << "\n");
//		m_testvectors[i] *= 1/e;

		//VecScaleAdd(m_testvectors[i], 1.0, m_testvectors[i], -1.0, d);
	}

	// output after smoothing of testvector
	if(m_famg.m_writeMatrices && m_famg.m_writeTestvectors)
			for(size_t i=0; i<m_testvectors.size(); i++)
				ConnectionViewer::WriteVectorPar((m_famg.m_writeMatrixPath + ToString("/testvector_S") + ToString(i) + ToString("_L") + ToString(level) + ".vec").c_str(),
					m_testvectors[i], &m_famg.m_amghelper.positions[level][0], m_famg.m_dbgDimension);


	if(m_bTiming) {UG_LOG(" took " << SW.ms() << " ms"); }
}
template<typename algebra_type>
void FAMGLevelCalculator<algebra_type>::get_testvectors_on_OL2()
{
#ifdef UG_PARALLEL
	for(size_t i=0; i<m_testvectors.size(); i++)
	{
		SetLayoutValues(&m_testvectors[i], A.layouts()->slave(), 0.0);
		m_testvectors[i].resize(A_OL2.num_rows());
		m_testvectors[i].set_layouts(A_OL2.layouts());
		for(size_t j=A.num_rows(); j<A_OL2.num_rows(); j++)
			m_testvectors[i][j]=0.0;

		m_testvectors[i].set_storage_type(PST_ADDITIVE);
		m_testvectors[i].change_storage_type(PST_CONSISTENT);
	}
#endif
}


template<typename algebra_type>
void FAMGLevelCalculator<algebra_type>::calculate_next_testvectors()
{
	vector_type t;
	t.resize(AH.num_rows());
	// todo: remove dynamic cast, change big_testvector to parallel
	for(size_t i=0; i<m_testvectors.size(); i++)
	{
		m_famg.injection(t, m_testvectors[i], level);
		m_testvectors[i] = t;
	}
}

// end group plugins_amg_famg
/// \}

}
#endif // __H__LIB_ALGEBRA__FAMG_SOLVER__FAMG_TESTVECTORS_H__
