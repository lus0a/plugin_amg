/**
 * \file famg.h
 *
 * \author Martin Rupp
 *
 * \date 16.11.2010
 *
 * class declaration for famg
 *
 * Goethe-Center for Scientific Computing 2010.
 */


#ifndef __H__UG__LIB_DISC__AMG_SOLVER__FAMG_H__
#define __H__UG__LIB_DISC__AMG_SOLVER__FAMG_H__

#include <vector>
#include <iostream>

#include "../amg_base.h"
#include "../amg_profiling.h"
#ifdef UG_PARALLEL
#include "../rsamg/rsamg_parallel_coarsening.h"
#endif

#include "testvector_calculators.h"
#include "debug_ids.h"

namespace ug{

/**
 * \defgroup plugins_amg_famg FAMG
 * \brief Filtering Algebraic Multigrid Functions.
 * \ingroup plugins_amg
 * \{
 */

template<typename TAlgebra>
class FAMGLevelCalculator;

template<typename TAlgebra>
class FAMGInterpolationCalculator;

// AMG
//---------------------------------
//! algebraic multigrid class.
//!

struct levelnode
{
	levelnode(size_t l, size_t n) {level =l; node =n; }
	size_t level;
	size_t node;
};


template <typename TAlgebra>
class FAMG:
	public AMGBase< TAlgebra >
{
public:
	typedef AMGBase<TAlgebra> super;
	using super::m_amghelper;
	using super::m_parentIndex;
	using super::m_writeMatrices;
	using super::m_writeMatrixPath;
	using super::levels;
	using super::m_dbgDimension;
	using super::m_dGalerkinTruncation;

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

	typedef AMGBase<TAlgebra> base_type;

//  functions
	FAMG() ;
	/// clone constructor
	FAMG(const FAMG<TAlgebra> &parent);
	virtual SmartPtr<ILinearIterator<vector_type> > clone()
	{
		SmartPtr<FAMG<algebra_type> > clone(new FAMG<algebra_type>(*this));
		return clone;
	}
	//	Name of preconditioner
	virtual ~FAMG();
	void cleanup();

	virtual const char* name() const {return "FAMGPreconditioner";}
	virtual std::string config_string() const;

	void set_aggressive_coarsening(bool bAggressiveCoarsening) { m_bAggressiveCoarsening = bAggressiveCoarsening; }
	bool get_aggressive_coarsening() const { return m_bAggressiveCoarsening; }

	void set_damping_for_smoother_in_interpolation_calculation(double dDampingForSmootherInInterpolationCalculation)
	{
		m_dDampingForSmootherInInterpolationCalculation = dDampingForSmootherInInterpolationCalculation;
	}

	void set_use_constant_if_failed(bool b)
	{
		m_bUseConstantIfFailed = b;
	}

	void set_external_coarsening(bool bUse) { m_bExternalCoarsening = bUse; }
	void set_use_precalculate(bool bUse) { m_bUsePrecalculate = bUse; }
private:
	void set_or_copy(const FAMG<TAlgebra> *parent);

public:
#ifdef UG_PARALLEL
	void set_parallel_coarsening(IParallelCoarsening *pPC)
	{
		m_pParallelCoarsening = pPC;
	}
#endif

	double get_damping_for_smoother_in_interpolation_calculation() const { return m_dDampingForSmootherInInterpolationCalculation; }

	std::map<size_t, double> mDeltaLevel;
	double mDeltaAll;
	bool m_bHAsInverse;

	void set_H_as_inverse(bool b) { m_bHAsInverse = b; }
	//! "Interpolation quality" F may not be worse than this (F < m_delta)
	void set_delta(double delta) { mDeltaLevel.clear(); mDeltaAll = delta; }
	void set_delta(double delta, size_t level) { mDeltaLevel[level] = delta; }

	//! with multiple parents, discard pairs with m_theta * F > min F.
	void set_theta(double theta) { m_theta = theta; }

	double get_delta() const { return m_delta;}
	double get_theta() const { return m_theta;}


	enum{ WRITE_COARSENING, READ_COARSENING, NONE} m_fileCoarsening;

	void set_write_coarsening() { m_fileCoarsening=WRITE_COARSENING;}
	void set_read_coarsening() { m_fileCoarsening=READ_COARSENING;}

	void set_testvector_smooths(size_t testvectordamps) { m_iTestvectorDamps = testvectordamps; }
	size_t get_testvector_smooths() const { return m_iTestvectorDamps; }

	// void add_std_testvector
	// with callback calculate(vector_type &c, double &weight, stdvector<positions>, stdvector<bool> isinner) o�.

	void set_write_testvectors(bool wt)
	{
		m_writeTestvectors = wt;
	}

	void init_for_EVP(const vector_type &u, const vector_type &b);


	void set_testvector_smoother(SmartPtr<ILinearIterator<vector_type> > testvectorsmoother) { m_testvectorsmoother = testvectorsmoother; }

	//!		sets epsilon_trunction, used in truncation of the interpolation [AMGKS99] 7.2.4
	//!		prolongation truncation: set all P(i,j) = 0 if they are too small compared to greatest P(i,.)
	void 	set_prolongation_truncation(double prolongationTr) 	{ m_dProlongationTruncation = prolongationTr; }
	double 	get_prolongation_truncation() const				{ return m_dProlongationTruncation; }

	void 	set_prereduce_A_parameter(double d) 	{	m_dPrereduceAToStrongParameter = d; }
	double 	get_prereduce_A_parameter() const		{	return m_dPrereduceAToStrongParameter; }

	void 	set_H_reduce_interpolation_nodes_parameter(double d)	{	m_dHReduceInterpolationNodesParameter = d; 	}
	double	get_H_reduce_interpolation_nodes_parameter() const		{	return m_dHReduceInterpolationNodesParameter; }

	void 	set_H_reduce_min_parents(int min)	{ m_iHReduceMinParents = min; }
	int	get_H_reduce_min_parents() const	{ return m_iHReduceMinParents; }

	void	set_strong_connection_external(double d) 	{ m_dStrongConnectionExternal = d; }
	double	get_strong_connection_external() 			{ return m_dStrongConnectionExternal; }

	void	set_strong_connection_submask(int m) 	{ m_iStrongConnectionSubMask = m; }
	int		get_strong_connection_submask() 		{ return m_iStrongConnectionSubMask; }
	const char *get_strong_connection_submask_desc() const
	{
		switch (m_iStrongConnectionSubMask) {
		case 0: 	return "submatrix(0)";
		case 1: 	return "submatrix(1)";
		case 2: 	return "submatrix(2)";
		case -2: 	return "full primary matrix (norm based)";
		default: 	return "standard trace";
		}

	}

	void	set_prevent_ff_connection_external(int i) 	{ m_iPreventFFConnection = i; }
	double	get_prevent_ff_connection_external() 		{ return m_iPreventFFConnection; }

	void	set_write_f_values(bool b)	{ m_bWriteFValues = b; }
	void	set_injection_restriction(bool b) { m_bInjectionRestriction = b; }

	void 	set_max_parents(size_t n) 	{ m_iMaxParents = n;}

	void	set_XC_LE_griebel(bool b)
	{
		m_bLE_XC_Griebel = b;
	}

	void	set_indirect_interpolation(bool b)
	{
		m_bIndirectInterpolation = b;
	}

	void	set_XC_new(bool b)
	{
		m_bXC_new = b;
	}

	void check_vector(SmartPtr<ITestvectorCalculator<TAlgebra> > testvectorCalculator, std::string desc);

	bool check_testvector(size_t fromlevel, size_t tolevel);
	bool check_testvector(size_t fromlevel)
	{
		return check_testvector(fromlevel, AMGBase<TAlgebra>::m_usedLevels);
	}
	bool check_testvector()
	{
		return check_testvector(0);
	}

	virtual bool needs_positions() const
	{
		if(m_bWriteFValues || m_writeTestvectors) return true;
		UG_COND_THROW(m_testvectorCalculator == SPNULL, "please provide testvectors.");
		return m_testvectorCalculator->needs_positions();
	}




	void set_three_parents_penalty(double d)
	{
		m_threeParentsPenalty = d;
	}

public:
	void set_testvectors(SmartPtr<ITestvectorCalculator<TAlgebra> > testvectorCalculator)
	{
		m_testvectorCalculator = testvectorCalculator;
	}

	SmartPtr<ITestvectorCalculator<TAlgebra> > get_testvector_calculator() const { return m_testvectorCalculator; }


	void add_watchnode(size_t level, size_t node)
	{
		debugWatchNodes.push_back(levelnode(level, node));
	}


	void set_deficit_solve_tolerance(double d)
	{
		m_dDeficitSolveTolerance = d;
	}

	void adaptive_filtering(SmartPtr<ILinearOperator<vector_type, vector_type> > linOp, vector_type &x, SmartPtr<IConvergenceCheck<vector_type> > spConvCheck,
				double expectedRate)
	{
		SmartPtr<TestvectorsConstant<TAlgebra> > spTV = make_sp(new TestvectorsConstant<TAlgebra>());
		adaptive_filtering(linOp, x, spConvCheck, expectedRate, spTV);
	}

	void adaptive_filtering(SmartPtr<ILinearOperator<vector_type, vector_type> > linOp, vector_type &x, SmartPtr<IConvergenceCheck<vector_type> > spConvCheck,
			double expectedRate, SmartPtr<ITestvectorCalculator<TAlgebra> > spTV);

	std::vector<levelnode> debugWatchNodes;

private:
//  functions
	virtual void create_AMG_level(matrix_type &AH, prolongation_matrix_type &R, const matrix_type &A,
				prolongation_matrix_type &P, size_t level, AMGLayouts &nextLevelLayouts);

	virtual void precalc_level(size_t level);
private:
// data
	SmartPtr<ITestvectorCalculator<TAlgebra> > m_testvectorCalculator;
	double m_dDampingForSmootherInInterpolationCalculation;
	double m_delta;				///< "Interpolation quality" F may not be worse than this (F < m_delta)
	double m_theta;				///< with multiple parents, discard pairs with m_theta * F > min F.

	double m_dProlongationTruncation;		///< parameter used for truncation of interpolation

	double m_dHReduceInterpolationNodesParameter; ///< parameter used to reduce the number of interpolation nodes by looking at the H(i,j) value.
	int m_iHReduceMinParents; ///< parameter used to reduce the number of interpolation nodes by looking at the H(i,j) value.

	double m_dPrereduceAToStrongParameter;	///< parameter used to reduce the matrix A to


	double m_dStrongConnectionExternal;
	int m_iStrongConnectionSubMask;

	int m_iPreventFFConnection;  ///<


	size_t m_iTestvectorDamps;
	bool m_writeTestvectors;

	size_t m_iAutoTestvectorsSize;
	double m_threeParentsPenalty;
		double m_dDeficitSolveTolerance;

	bool m_bAggressiveCoarsening;
	bool m_bExternalCoarsening;
	bool m_bUsePrecalculate;
	bool m_bWriteFValues;
	bool m_bInjectionRestriction;

	bool m_bXC_new;				///< new strategy, \sa rs_amg_coarsening_new
	bool m_bLE_XC_Griebel;		///< Griebel LE strategy for neumann-coarsening \sa rs_amg_coarsening_new

	bool m_bUseConstantIfFailed;

	SmartPtr<ILinearIterator<vector_type> > m_testvectorsmoother;

	friend class FAMGLevelCalculator<algebra_type>;
	friend class FAMGInterpolationCalculator<algebra_type>;

	stdvector< vector_type > testvectors;
	stdvector<double> omega;

	bool m_bXCaggressive;
	bool m_bIndirectInterpolation;

	vector_type evpU, evpB;
	bool m_bProjectedEVP;

#ifdef UG_PARALLEL
	IParallelCoarsening *m_pParallelCoarsening;
#endif
public:
	size_t m_iMaxParents;
};

// end group plugins_amg_famg
/// \}

} // namespace ug


#include "adaptive_filtering.h"
#include "famg_impl.h"

#endif // __H__UG__LIB_DISC__AMG_SOLVER__FAMG_H__
