/*
 * Copyright (c) 2010-2014:  G-CSC, Goethe University Frankfurt
 * Author: Arne Nägel
 * 
 * This file is part of UG4.
 * 
 * UG4 is free software: you can redistribute it and/or modify it under the
 * terms of the GNU Lesser General Public License version 3 (as published by the
 * Free Software Foundation) with the following additional attribution
 * requirements (according to LGPL/GPL v3 §7):
 * 
 * (1) The following notice must be displayed in the Appropriate Legal Notices
 * of covered and combined works: "Based on UG4 (www.ug4.org/license)".
 * 
 * (2) The following notice must be displayed at a prominent place in the
 * terminal output of covered works: "Based on UG4 (www.ug4.org/license)".
 * 
 * (3) The following bibliography is recommended for citation and must be
 * preserved in all covered files:
 * "Reiter, S., Vogel, A., Heppner, I., Rupp, M., and Wittum, G. A massively
 *   parallel geometric multigrid solver on hierarchically distributed grids.
 *   Computing and visualization in science 16, 4 (2013), 151-164"
 * "Vogel, A., Reiter, S., Rupp, M., Nägel, A., and Wittum, G. UG4 -- a novel
 *   flexible software system for simulating pde based models on high performance
 *   computers. Computing and visualization in science 16, 4 (2013), 165-179"
 * 
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 * GNU Lesser General Public License for more details.
 */


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
