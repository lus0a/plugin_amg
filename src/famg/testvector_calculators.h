/*
 * testvector_calculators.h
 *
 *  Created on: 07.02.2014
 *      Author: mrupp
 */

#ifndef TESTVECTOR_CALCULATORS_H_
#define TESTVECTOR_CALCULATORS_H_
namespace ug{

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

template<typename TAlgebra>
class ITestvectorCalculator
{
	typedef typename TAlgebra::matrix_type matrix_type;
	typedef typename TAlgebra::vector_type vector_type;
	typedef MatrixOperator<matrix_type, vector_type> matrix_operator_type;
public:
	virtual void calc(SmartPtr<matrix_operator_type> A, MathVector<3> *pos, std::vector<vector_type> &testvectors) = 0;
	virtual std::string name() const  =0 ;
	virtual std::string config_string() const { return name(); }
	virtual bool needs_positions() const = 0;
	virtual ~ITestvectorCalculator() {}
};


template<typename TAlgebra>
class TestvectorsConstant : public ITestvectorCalculator<TAlgebra>
{
	typedef typename TAlgebra::matrix_type matrix_type;
	typedef typename TAlgebra::vector_type vector_type;
	typedef MatrixOperator<matrix_type, vector_type> matrix_operator_type;
public:
	virtual void calc(SmartPtr<matrix_operator_type> A, MathVector<3> *pos, std::vector<vector_type> &testvectors)
	{
		testvectors.resize(1);
		vector_type &v = testvectors[0];
		v.resize(A->num_rows());
		for(size_t i=0; i<A->num_rows(); i++)
		{
//			if(A->is_isolated(i))
//				v[i] = 0.0;
//			else
				v[i] = 1.0;
		}
#ifdef UG_PARALLEL
		v.set_storage_type(PST_CONSISTENT);
#endif
	}
	virtual bool needs_positions() const { return false; }
	virtual std::string name() const { return "TestvectorsConstant"; }
};

template<typename TAlgebra>
class TestvectorsConstantWobbled : public ITestvectorCalculator<TAlgebra>
{
	typedef typename TAlgebra::matrix_type matrix_type;
	typedef typename TAlgebra::vector_type vector_type;
	typedef MatrixOperator<matrix_type, vector_type> matrix_operator_type;
public:
	double m_delta;
	TestvectorsConstantWobbled(double delta)
	{
		m_delta = delta;
	}

	virtual void calc(SmartPtr<matrix_operator_type> A, MathVector<3> *pos, std::vector<vector_type> &testvectors)
	{
		testvectors.resize(1);
		vector_type &v = testvectors[0];
		v.resize(A->num_rows());
		for(size_t i=0; i<A->num_rows(); i++)
		{
			/*if(A->is_isolated(i))
				v[i] = 0.0;
			else*/
				v[i] = 1.0+urand(-m_delta, m_delta);
		}
#ifdef UG_PARALLEL
		v.set_storage_type(PST_ADDITIVE);
		v.change_storage_type(PST_CONSISTENT);
#endif
	}
	virtual bool needs_positions() const { return false; }
	virtual std::string name() const { return "TestvectorsConstant"; }
};

template<typename TAlgebra>
class TestvectorsLinearElasticity2d : public ITestvectorCalculator<TAlgebra>
{
	typedef typename TAlgebra::matrix_type matrix_type;
	typedef typename TAlgebra::vector_type vector_type;
	typedef MatrixOperator<matrix_type, vector_type> matrix_operator_type;
public:
	virtual void calc(SmartPtr<matrix_operator_type> A, MathVector<3> *pos, std::vector<vector_type> &testvectors)
	{
		UG_COND_THROW(pos == NULL, "TestvectorsLinearElasticity2d needs positions!");
		size_t N = A->num_rows();
		testvectors.resize(2);
		for(size_t i=0; i<testvectors.size(); i++)
			testvectors[i].resize(N);

		for(size_t i=0; i<N; i++)
		{
			testvectors[0][i] = 0.0;
			testvectors[1][i] = 0.0;

//			if(!A->is_isolated(i))
			{
				testvectors[0][i] = 1.0;
				SetAlgVector2d(testvectors[1][i], GetRotZ(pos[i]));
			}
		}

#ifdef UG_PARALLEL
		for(size_t i=0; i<testvectors.size(); i++)
			testvectors[i].set_storage_type(PST_CONSISTENT);
#endif
	}
	virtual bool needs_positions() const { return true; }
	virtual std::string name() const { return "TestvectorsLinearElasticity2d"; }
	virtual std::string config_string() const { return "TestvectorsLinearElasticity2d: Rigid Body Modes of Linear Elasticity 2d: (1, 0), (0, 1), (-y, x)"; }
};

template<typename TAlgebra>
class TestvectorsOnlyRotation2d : public ITestvectorCalculator<TAlgebra>
{
	typedef typename TAlgebra::matrix_type matrix_type;
	typedef typename TAlgebra::vector_type vector_type;
	typedef MatrixOperator<matrix_type, vector_type> matrix_operator_type;
public:
	virtual void calc(SmartPtr<matrix_operator_type> A, MathVector<3> *pos, std::vector<vector_type> &testvectors)
	{
		UG_COND_THROW(pos == NULL, "TestvectorsOnlyRotation2d needs positions!");
		size_t N = A->num_rows();
		testvectors.resize(1);
		for(size_t i=0; i<testvectors.size(); i++)
			testvectors[i].resize(N);

		for(size_t i=0; i<N; i++)
		{
			testvectors[1][0] = 0.0;
//			if(!A->is_isolated(i))
				SetAlgVector2d(testvectors[0][i], GetRotZ(pos[i]));
		}
#ifdef UG_PARALLEL
		for(size_t i=0; i<testvectors.size(); i++)
			testvectors[i].set_storage_type(PST_CONSISTENT);
#endif
	}

	virtual bool needs_positions() const { return true; }
	virtual std::string name() const { return "TestvectorsOnlyRotation2d"; }
	virtual std::string config_string() const { return "TestvectorsOnlyRotation2d: (-y, x)"; }
};

template<typename TAlgebra>
class TestvectorsLinearElasticity3d : public ITestvectorCalculator<TAlgebra>
{
	typedef typename TAlgebra::matrix_type matrix_type;
	typedef typename TAlgebra::vector_type vector_type;
	typedef MatrixOperator<matrix_type, vector_type> matrix_operator_type;
public:
	virtual void calc(SmartPtr<matrix_operator_type> A, MathVector<3> *pos, std::vector<vector_type> &testvectors)
	{
		UG_COND_THROW(pos == NULL, "TestvectorsLinearElasticity3d needs positions!");
		size_t N = A->num_rows();
		testvectors.resize(4);
		for(size_t i=0; i<testvectors.size(); i++)
			testvectors[i].resize(N);

		for(size_t i=0; i<N; i++)
		{
			for(size_t j=0; j<4; j++)
				testvectors[j][i] = 0.0;

//			if(!A->is_isolated(i))
			{
				testvectors[0][i] = 1.0;
				SetAlgVector3d(testvectors[1][i], GetRotX(pos[i]));
				SetAlgVector3d(testvectors[2][i], GetRotY(pos[i]));
				SetAlgVector3d(testvectors[3][i], GetRotZ(pos[i]));
//				BlockRef(testvectors[0][i], 0) = 1.0;
//				BlockRef(testvectors[1][i], 1) = 1.0;
//				BlockRef(testvectors[2][i], 2) = 1.0;
//				SetAlgVector3d(testvectors[3][i], GetRotX(pos[i]));
//				SetAlgVector3d(testvectors[4][i], GetRotY(pos[i]));
//				SetAlgVector3d(testvectors[5][i], GetRotZ(pos[i]));
			}
		}
#ifdef UG_PARALLEL
		for(size_t i=0; i<testvectors.size(); i++)
			testvectors[i].set_storage_type(PST_CONSISTENT);
#endif
	}
	virtual bool needs_positions() const { return true; }
	virtual std::string name() const { return "TestvectorsLinearElasticity3d"; }
	virtual std::string config_string() const { return "TestvectorsLinearElasticity2d: Rigid Body Modes of Linear Elasticity 3d: (1, 0, 0), (0, 1, 0), (0, 0, 1), (0, -z, y), (-z, 0, x), (-y, x, 0)]"; }
};





template<typename TAlgebra>
class TestvectorsLinear3d : public ITestvectorCalculator<TAlgebra>
{
	typedef typename TAlgebra::matrix_type matrix_type;
	typedef typename TAlgebra::vector_type vector_type;
	typedef MatrixOperator<matrix_type, vector_type> matrix_operator_type;
public:
	virtual void calc(SmartPtr<matrix_operator_type> A, MathVector<3> *pos, std::vector<vector_type> &testvectors)
	{
		UG_COND_THROW(pos == NULL, "TestvectorsLinear3d needs positions!");
		size_t N = A->num_rows();
		testvectors.resize(4);
		for(size_t i=0; i<testvectors.size(); i++)
			testvectors[i].resize(N);
		for(size_t i=0; i<N; i++)
		{
			testvectors[0][i] = 0.0;
			testvectors[1][i] = 0.0;
			testvectors[2][i] = 0.0;
			testvectors[3][i] = 0.0;
//			if(!A->is_isolated(i))
			{
				testvectors[0][i] = 1.0;
				testvectors[1][i] = pos[i].x();
				testvectors[2][i] = pos[i].y();
				testvectors[3][i] = pos[i].z();
			}
		}
#ifdef UG_PARALLEL
		for(size_t i=0; i<testvectors.size(); i++)
			testvectors[i].set_storage_type(PST_CONSISTENT);
#endif
	}
	virtual bool needs_positions() const { return true; }
	virtual std::string name() const { return "TestvectorsLinear3d"; }
	virtual std::string config_string() const { return "TestvectorsLinear3d: constants and linear functions: 1, x, y, z in each component"; }
};


template<typename TAlgebra>
class TestvectorsLinear2d : public ITestvectorCalculator<TAlgebra>
{
	typedef typename TAlgebra::matrix_type matrix_type;
	typedef typename TAlgebra::vector_type vector_type;
	typedef MatrixOperator<matrix_type, vector_type> matrix_operator_type;
public:
	virtual void calc(SmartPtr<matrix_operator_type> A, MathVector<3> *pos, std::vector<vector_type> &testvectors)
	{
		UG_COND_THROW(pos == NULL, "TestvectorsLinear2d needs positions!");
		size_t N = A->num_rows();
		testvectors.resize(3);
		for(size_t i=0; i<testvectors.size(); i++)
			testvectors[i].resize(N);

		for(size_t i=0; i<N; i++)
		{
			testvectors[0][i] = 0.0;
			testvectors[1][i] = 0.0;
			testvectors[2][i] = 0.0;
//			if(!A->is_isolated(i))
			{
				testvectors[0][i] = 1.0;
				testvectors[1][i] = pos[i].x();
				testvectors[2][i] = pos[i].y();
			}
		}
#ifdef UG_PARALLEL
		for(size_t i=0; i<testvectors.size(); i++)
			testvectors[i].set_storage_type(PST_CONSISTENT);
#endif
	}
	virtual bool needs_positions() const { return true; }
	virtual std::string name() const { return "TestvectorsLinear2d"; }
	virtual std::string config_string() const { return "TestvectorsLinear2d: constants and linear functions: 1, x, y in each component"; }
};


template<typename TAlgebra>
class TestvectorsLinearComponent : public ITestvectorCalculator<TAlgebra>
{
	typedef typename TAlgebra::matrix_type matrix_type;
	typedef typename TAlgebra::vector_type vector_type;
	typedef MatrixOperator<matrix_type, vector_type> matrix_operator_type;
public:
	virtual void calc(SmartPtr<matrix_operator_type> A, MathVector<3> *pos, std::vector<vector_type> &testvectors)
	{
		UG_COND_THROW(pos == NULL, "TestvectorsLinearComponent needs positions!");
		size_t N = A->num_rows();
		testvectors.resize(2);
		for(size_t i=0; i<testvectors.size(); i++)
			testvectors[i].resize(N);

		bool b3 = GetSize(testvectors[0][0]) == 3;
		for(size_t i=0; i<N; i++)
		{
			testvectors[0][i] = 0.0;
//			if(!A->is_isolated(i))
			{
				testvectors[0][i] = 1.0;
				BlockRef(testvectors[1][i], 0) = pos[i].x();
				BlockRef(testvectors[1][i], 1) = pos[i].y();
				if(b3)
					BlockRef(testvectors[1][i], 2) = pos[i].z();
			}
		}
#ifdef UG_PARALLEL
		for(size_t i=0; i<testvectors.size(); i++)
			testvectors[i].set_storage_type(PST_CONSISTENT);
#endif
	}
	virtual bool needs_positions() const { return true; }
	virtual std::string name() const { return "TestvectorsLinearComponent"; }
	virtual std::string config_string() const { return "TestvectorsLinear2d: constants and linear functions: 1 and x_i in each component i"; }
};

template<typename TAlgebra>
class TestvectorsAuto : public ITestvectorCalculator<TAlgebra>
{
	typedef typename TAlgebra::matrix_type matrix_type;
	typedef typename TAlgebra::vector_type vector_type;
	typedef MatrixOperator<matrix_type, vector_type> matrix_operator_type;

	SmartPtr<ILinearIterator<vector_type> > m_smoother;
	size_t m_numSmooths;
	size_t m_numTestvectors;
public:
	TestvectorsAuto(size_t numTestvectors, SmartPtr<ILinearIterator<vector_type> > smoother, size_t numSmooths)
	{
		m_numTestvectors = numTestvectors;
		m_numSmooths = numSmooths;
		m_smoother = smoother;
	}

	virtual void calc(SmartPtr<matrix_operator_type> A, MathVector<3> *pos, std::vector<vector_type> &testvectors)
	{
		size_t N = A->num_rows();
		testvectors.resize(m_numTestvectors);

		// since we only need the smoother in this function, we want the smoother
		// to clear himself up after this function
		SmartPtr<ILinearIterator<vector_type> > smoother = m_smoother->clone();
		smoother->init(A);

		vector_type x; x.resize(N);
		vector_type d; d.resize(N);
		vector_type c; c.resize(N);
	#ifdef UG_PARALLEL
		x.set_layouts(A->layouts());
		c.set_layouts(A->layouts());
		d.set_layouts(A->layouts());
	#endif

		for(size_t i=0; i<testvectors.size(); i++)
		{
	#ifdef UG_PARALLEL
			testvectors[i].set_storage_type(PST_CONSISTENT);
			x.set_storage_type(PST_ADDITIVE);
	#endif
			x.set_random(-1.0, 1.0);
			A->apply(d, x);

			for(size_t jj=0; jj < m_numSmooths; jj++)
			{
				c.set(0.0);
				smoother->apply(c, d);
				x -= c;
				A->apply_sub(d, c);
			}
			testvectors[i] = x;
		}

#ifdef UG_PARALLEL
		for(size_t i=0; i<testvectors.size(); i++)
			testvectors[i].change_storage_type(PST_CONSISTENT);
#endif
	}
	virtual bool needs_positions() const { return false; }
	virtual std::string name() const { return "TestvectorsAuto"; }
	virtual std::string config_string() const
	{
		std::stringstream ss;
		ss << "TestvectorsAuto: numTestvectors = " << m_numTestvectors << ", " << m_numSmooths << "x ";
		ss << ConfigShift(m_smoother->config_string());
		return ss.str();
	}
};

/*
 * for external:
 * template<typename TAlgebra>
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


	stdvector<double> m_omegaVectors; // testvectorWeights
	stdvector< IVectorWriter<vector_type> * > m_vVectorWriters;
	stdvector<double> m_omegaVectorWriters; // testvectorWeights


}
 *
 * 	void add_testvector(vector_type& c, double weight)
	{
		m_iTestvectors = TESTVECTORS_EXTERNAL;
		UG_ASSERT(m_testvectors.size() == m_omegaVectors.size(), "???");
		m_testvectors.push_back(c);
		m_omegaVectors.push_back(weight);
	}

	void add_testvector(IVectorWriter<vector_type> *vw, double weight)
	{
		m_iTestvectors = TESTVECTORS_EXTERNAL;
		m_vVectorWriters.push_back(vw);
		m_omegaVectorWriters.push_back(weight);
	}*/

}
 #endif /* TESTVECTOR_CALCULATORS_H_ */
