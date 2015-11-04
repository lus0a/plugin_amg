#ifndef __H__UG__LIB_ALGEBRA__AMG_SOLVER__SEMIAMG_H__
#define __H__UG__LIB_ALGEBRA__AMG_SOLVER__SEMIAMG_H__

#include <vector>
#include <iostream>

#include "../amg_base.h"
#include "../rsamg/rsamg_nodeinfo.h"

//#include "semi_amg_strategy.h"

#include "tensor_positions.h"

#ifdef UG_PARALLEL
//#include "../rsamg/rsamg_parallel_coarsening.h"
#endif


namespace ug{

/**
 * \defgroup plugins_amg_rsamg SemiAMG
 * \brief Algebraic Multigrid Functions.
 * \ingroup plugins_amg
 * \{
 */

// AMG
//---------------------------------
/**
 * algebraic multigrid (AMG) class.
 *
 * scientific references:
 * [AMGKS99]: Algebraic Multigrid (AMG): An Introduction with Applications, K Stueben. GMD Report 70, November 1999
 *
 */
template <typename TAlgebra>
class SemiAMG : public AMGBase<TAlgebra>
{
public:
	typedef AMGBase<TAlgebra> super;
	using super::m_amghelper;
	using super::m_parentIndex;
	using super::m_writeMatrices;
	using super::m_writeMatrixPath;
	using super::levels;
	using super::m_dGalerkinTruncation;
	using super::m_verbose;


public:
	typedef TAlgebra algebra_type;                      ///< Algebra type
	typedef typename TAlgebra::vector_type vector_type; ///< Vector type
 	typedef typename TAlgebra::matrix_type matrix_type; ///< Matrix type
	typedef typename TAlgebra::matrix_type prolongation_matrix_type;
	typedef typename matrix_type::value_type value_type;


//  functions
	SemiAMG();
	virtual SmartPtr<ILinearIterator<vector_type> > clone()
	{
		SmartPtr<SemiAMG<algebra_type> > clone(new SemiAMG<algebra_type>);
		return clone;
	}
	virtual ~SemiAMG() { }

	virtual const char* name() const {return "SemiAMGPreconditioner";}
	virtual void precalc_level(size_t level) { }

	void tostring() const;

	virtual bool needs_positions() const { return true; }

	/// set tensor axis
	void set_tensor_axis(unsigned int axis, const stdvector<double> &coord0, const stdvector<int> &level)
	{ m_tli.set_axis_level(axis, level, coord0);}

	static const bool bTiming =true;
	void set_virtual_nodes(bool flag)
	{ m_bCreateVirtualNodes = flag;}
protected:
//  functions
	virtual void create_AMG_level(matrix_type &AH, prolongation_matrix_type &R, const matrix_type &A,
					prolongation_matrix_type &P, size_t level, AMGLayouts &nextLevelLayouts);

	size_t user_defined_coarsening(const cgraph &Sl, size_t level, AMGNodes &rating);
	void user_defined_prolongation(const cgraph &Sl, const matrix_type &Al, prolongation_matrix_type &P, size_t level, AMGNodes &rating);



#ifdef UG_PARALLEL

#endif

// data
	TensorLevelInformation m_tli;
	double m_dTheta;
	//semi_amg_user_strategy m_strategy;
	typedef std::map<MathVector<3>, size_t> TVirtualCoarseMap;
	typedef std::pair<MathVector<3>, size_t> TVirtualCoarseMapPair;

	TVirtualCoarseMap vcoarse[20];
	size_t vcount[20];
	bool m_bCreateVirtualNodes;

};
	

// end group plugins_amg_rsamg
/// \}

} // namespace ug


#include "semi_amg_impl.h"

#endif // __H__UG__LIB_ALGEBRA__AMG_SOLVER__SEMIAMG_H__
