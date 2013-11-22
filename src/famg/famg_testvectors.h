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
void FAMGLevelCalculator<algebra_type>::calculate_testvectors()
{
	UG_SET_DEBUG_LEVEL(LIB_ALG_AMG, m_famg.iDebugLevelTestvectorCalc);

	// todo: all global?

	Stopwatch SW;
	if(m_bTiming) { UG_LOG("\ncalculating testvector... "); SW.start(); }

	vector_type d; d.resize(A.num_rows());
	vector_type c; c.resize(A.num_rows());
#ifdef UG_PARALLEL
	c.set_layouts(A.layouts());
	d.set_layouts(A.layouts());
#endif

	if(m_famg.m_writeMatrices && m_famg.m_writeTestvectors)
		for(size_t i=0; i<m_testvectors.size(); i++)
			ConnectionViewer::WriteVector((m_famg.m_writeMatrixPath + ToString("/testvector_pS") + ToString(i) + ToString("_L") + ToString(level) + ".vec").c_str(),
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
			m_famg.m_testvectorsmoother->apply(c, d);
			m_testvectors[i] -= c;
			if(false && m_famg.m_writeMatrices && m_famg.m_writeTestvectors)
			{
					ConnectionViewer::WriteVector((m_famg.m_writeMatrixPath + ToString("/testvector_St") + ToString(i) + "____" + ToString(jj) + ToString("_L") + ToString(level) + ".vec").c_str(),
						m_testvectors[i], &m_famg.m_amghelper.positions[level][0], m_famg.m_dbgDimension);
					ConnectionViewer::WriteVector((m_famg.m_writeMatrixPath + ToString("/testvector_Sc") + ToString(i) + "____" + ToString(jj) + ToString("_L") + ToString(level) + ".vec").c_str(),
										c, &m_famg.m_amghelper.positions[level][0], m_famg.m_dbgDimension);
					ConnectionViewer::WriteVector((m_famg.m_writeMatrixPath + ToString("/testvector_Sd") + ToString(i) + "____" + ToString(jj) + ToString("_L") + ToString(level) + ".vec").c_str(),
										d, &m_famg.m_amghelper.positions[level][0], m_famg.m_dbgDimension);
			}

		}

		//double e=EnergyProd(m_testvectors[i], A, m_testvectors[i]);;
		//UG_LOG("EnergyProd " << i << " = " << e << "\n");
//		m_testvectors[i] *= 1/e;

		//VecScaleAdd(m_testvectors[i], 1.0, m_testvectors[i], -1.0, d);
	}

	if(m_famg.m_writeMatrices && m_famg.m_writeTestvectors)
			for(size_t i=0; i<m_testvectors.size(); i++)
				ConnectionViewer::WriteVector((m_famg.m_writeMatrixPath + ToString("/testvector_S") + ToString(i) + ToString("_L") + ToString(level) + ".vec").c_str(),
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


/**
 * @return [0, z, -y], the rotation around the x-axes
 */
MathVector<3> GetRotX(MathVector<3> &pos)
{
	return MathVector<3>(0, pos.z(), -pos.y());
}

/**
 * @return [z, 0, -x], the rotation around the x-axes
 */
MathVector<3> GetRotY(MathVector<3> &pos)
{
	return MathVector<3>(pos.z(), 0, -pos.x());
}

/**
 * @return [y, -x, 0], the rotation around the z-axes
 */
MathVector<3> GetRotZ(MathVector<3> &pos)
{
	return MathVector<3>(pos.y(), -pos.x(), 0);
}

template<typename T>
inline void SetAlgVector2d(T &v, MathVector<3> mv)
{
	BlockRef(v, 0) = mv.x();
	BlockRef(v, 1) = mv.y();
}

template<typename T>
inline void SetAlgVector3d(T &v, MathVector<3> mv)
{
	BlockRef(v, 0) = mv.x();
	BlockRef(v, 1) = mv.y();
	BlockRef(v, 2) = mv.z();
}


template<typename TMatrix, typename TVector>
void GetTestvectorsLE_2d	(const TMatrix &A, stdvector<TVector> &testvectors, stdvector<double> &omega, stdvector<MathVector<3> > &pos,
		double weightC, double weightZ)
{
	testvectors.resize(3);
	for(size_t i=0; i<testvectors.size(); i++)
		testvectors[i].resize(A.num_rows());

	omega.resize(3);
	omega[0] = weightC;
	omega[1] = weightC;
	omega[2] = weightZ;

	for(size_t i=0; i<A.num_rows(); i++)
	{
		testvectors[0][i] = 0.0;
		testvectors[1][i] = 0.0;
		testvectors[2][i] = 0.0;

		if(!A.is_isolated(i))
		{
			BlockRef(testvectors[0][i], 0) = 1.0;
			BlockRef(testvectors[1][i], 1) = 1.0;
			SetAlgVector2d(testvectors[2][i], GetRotZ(pos[i]));
		}
	}
}


template<typename TMatrix, typename TVector>
void GetTestvectorsLE_3d	(const TMatrix &A, stdvector<TVector> &testvectors, stdvector<double> &omega, stdvector<MathVector<3> > &pos,
		double weightC, double weightX, double weightY, double weightZ)

{
	testvectors.resize(6);
	for(size_t i=0; i<testvectors.size(); i++)
		testvectors[i].resize(A.num_rows());

	omega.resize(6);
	omega[0] = weightC;
	omega[1] = weightC;
	omega[2] = weightC;
	omega[3] = weightX;
	omega[4] = weightY;
	omega[5] = weightZ;

	for(size_t i=0; i<A.num_rows(); i++)
	{
		for(size_t j=0; j<6; j++)
			testvectors[j][i] = 0.0;

		if(!A.is_isolated(i))
		{
			BlockRef(testvectors[0][i], 0) = 1.0;
			BlockRef(testvectors[1][i], 1) = 1.0;
			BlockRef(testvectors[2][i], 2) = 1.0;
			SetAlgVector3d(testvectors[3][i], GetRotX(pos[i]));
			SetAlgVector3d(testvectors[4][i], GetRotY(pos[i]));
			SetAlgVector3d(testvectors[5][i], GetRotZ(pos[i]));
		}
	}
}

template<typename TMatrix, typename TVector>
void GetTestvectorsLinear_3d(const TMatrix &A, stdvector<TVector> &testvectors, stdvector<double> &omega, stdvector<MathVector<3> > &pos,
		double weightC, double weightX, double weightY, double weightZ)

{
	testvectors.resize(4);
	for(size_t i=0; i<testvectors.size(); i++)
		testvectors[i].resize(A.num_rows());

	omega.resize(4);
	omega[0] = weightC;
	omega[1] = weightX;
	omega[2] = weightY;
	omega[3] = weightZ;

	for(size_t i=0; i<A.num_rows(); i++)
	{
		testvectors[0][i] = 0.0;
		testvectors[1][i] = 0.0;
		testvectors[2][i] = 0.0;
		testvectors[3][i] = 0.0;
		if(!A.is_isolated(i))
		{
			testvectors[0][i] = 1.0;
			testvectors[1][i] = pos[i].x();
			testvectors[2][i] = pos[i].y();
			testvectors[3][i] = pos[i].z();
		}
	}
}


template<typename TMatrix, typename TVector>
void GetTestvectorsLinear_2d(const TMatrix &A, stdvector<TVector> &testvectors, stdvector<double> &omega, stdvector<MathVector<3> > &pos,
		double weightC, double weightX, double weightY)

{
	testvectors.resize(3);
	for(size_t i=0; i<testvectors.size(); i++)
		testvectors[i].resize(A.num_rows());

	omega.resize(3);
	omega[0] = weightC;
	omega[1] = weightX;
	omega[2] = weightY;

	for(size_t i=0; i<A.num_rows(); i++)
	{
		testvectors[0][i] = 0.0;
		testvectors[1][i] = 0.0;
		testvectors[2][i] = 0.0;
		if(!A.is_isolated(i))
		{
			testvectors[0][i] = 1.0;
			testvectors[1][i] = pos[i].x();
			testvectors[2][i] = pos[i].y();
		}
	}
}

// end group plugins_amg_famg
/// \}

}
#endif // __H__LIB_ALGEBRA__FAMG_SOLVER__FAMG_TESTVECTORS_H__
