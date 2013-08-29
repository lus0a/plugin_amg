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



template<typename matrix_type, typename vector_type>
double ScalProd(const vector_type &v1, const matrix_type &M, const vector_type &v2)
{
	double sum=0.0;
	UG_ASSERT(v1.size() == M.num_rows() && M.num_cols() == v2.size(), "size mismatch");
	for(size_t i=0; i < v1.size(); ++i)
	{
		for(typename matrix_type::const_row_iterator it = M.begin_row(i); it != M.end_row(i); ++it)
			sum += v1[i] * it.value() * v2[i];
	}
	return sum;
}

template<typename matrix_type, typename vector_type>
double EnergyProd(const vector_type &v1, const matrix_type &M)
{
	double sum=0.0;
	UG_ASSERT(v1.size() == M.num_rows() && M.num_cols() == v1.size(), "size mismatch");
	for(size_t i=0; i < v1.size(); ++i)
	{
		double s=0.0;
		for(typename matrix_type::const_row_iterator it = M.begin_row(i); it != M.end_row(i); ++it)
			s += v1[it.index()] * it.value();
		sum += v1[i]*s;
	}
	return sum;
}



template<typename matrix_type, typename vector_type>
void CalculateTestvector(const matrix_type &A_OL2, vector_type &big_testvector,
		size_t iTestvectorDamps)
{
	AMG_PROFILE_FUNC();
	vector_type d; d.resize(A_OL2.num_rows());
	for(size_t jj=0; jj < iTestvectorDamps; jj++)
	{
		MatMult(d, 1.0, A_OL2, big_testvector);
		for(size_t i=0; i<A_OL2.num_rows(); i++)
			big_testvector[i] = big_testvector[i] - 0.6*d[i]/A_OL2(i,i);
	}

	//big_testvector.set_storage_type(PST_CONSISTENT);
	//big_testvector.change_storage_type(PST_UNIQUE);
	//UG_LOG("Norm of testvector " << big_testvector.norm() << "\n");
	//big_testvector.change_storage_type(PST_CONSISTENT);
	//VecScaleAssign(big_testvector, 1/sqrt(EnergyProd(big_testvector, A_OL2)), big_testvector);
}


template<typename algebra_type>
void FAMGLevelCalculator<algebra_type>::calculate_testvectors()
{
	UG_SET_DEBUG_LEVEL(LIB_ALG_AMG, m_famg.iDebugLevelTestvectorCalc);

	// todo: all global?
	UG_DLOG(LIB_ALG_AMG, 1, "\ncalculating testvector... ");
	Stopwatch SW;
	if(bTiming) SW.start();

	vector_type d; d.resize(A.num_rows());
	vector_type c; c.resize(A.num_rows());
#ifdef UG_PARALLEL
	c.set_layouts(A.layouts());
	d.set_layouts(A.layouts());
#endif

	for(size_t i=0; i<m_testvectors.size(); i++)
	{

		UG_ASSERT(m_testvectors[i].size() == A.num_rows(), m_testvectors[i].size() << " != " << A.num_rows());

#ifdef UG_PARALLEL
		m_testvectors[i].set_storage_type(PST_CONSISTENT);
#endif

		m_famg.m_testvectorsmoother->init(m_famg.levels[level]->pAgglomeratedA);

		if(m_famg.m_/Matrices && m_famg.m_writeTestvectors)
					for(size_t i=0; i<m_testvectors.size(); i++)
						ConnectionViewer::WriteVector((m_famg.m_writeMatrixPath + ToString("/testvector_pS") + ToString(i) + ToString("_L") + ToString(level) + ".vec").c_str(),
							m_testvectors[i], &m_famg.m_amghelper.positions[level][0], m_famg.m_dbgDimension);

		for(size_t jj=0; jj < m_famg.get_testvector_smooths(); jj++)
		{
			d.set(0.0);
			c.set(0.0);
			A.apply(d, m_testvectors[i]);
			m_famg.m_testvectorsmoother->apply(c, d);
			m_testvectors[i] -= c;
		}

		if(m_famg.m_writeMatrices && m_famg.m_writeTestvectors)
			for(size_t i=0; i<m_testvectors.size(); i++)
				ConnectionViewer::WriteVector((m_famg.m_writeMatrixPath + ToString("/testvector_S") + ToString(i) + ToString("_L") + ToString(level) + ".vec").c_str(),
					m_testvectors[i], &m_famg.m_amghelper.positions[level][0], m_famg.m_dbgDimension);
		//VecScaleAdd(m_testvectors[i], 1.0, m_testvectors[i], -1.0, d);
	}


	if(bTiming) UG_DLOG(LIB_ALG_AMG, 1, "took " << SW.ms() << " ms");
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
	testvectors.resize(2);
	for(size_t i=0; i<testvectors.size(); i++)
		testvectors[i].resize(A.num_rows());

	omega.resize(2);
	omega[0] = weightC;
	omega[1] = weightZ;

	for(size_t i=0; i<A.num_rows(); i++)
	{
		if(A.is_isolated(i))
		{
			testvectors[0][i] = 0.0;
			testvectors[1][i] = 0.0;
		}
		else
		{
			testvectors[0][i] = 1.0;
			SetAlgVector2d(testvectors[1][i], GetRotZ(pos[i]));
		}
	}
}


template<typename TMatrix, typename TVector>
void GetTestvectorsLE_3d	(const TMatrix &A, stdvector<TVector> &testvectors, stdvector<double> &omega, stdvector<MathVector<3> > &pos,
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
			SetAlgVector3d(testvectors[1][i], GetRotX(pos[i]));
			SetAlgVector3d(testvectors[2][i], GetRotY(pos[i]));
			SetAlgVector3d(testvectors[3][i], GetRotZ(pos[i]));
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
