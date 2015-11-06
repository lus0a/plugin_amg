/*
 * Copyright (c) 2012-2015:  G-CSC, Goethe University Frankfurt
 * Author: Martin Rupp
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


#include <iostream>
#include <sstream>
#include <string>

//#define UG_USE_AMG // temporary switch until AMG for systems works again

// include bridge
#include "bridge/bridge.h"
#include "bridge/util.h"
#include "bridge/util_algebra_dependent.h"

#include "lib_algebra/lib_algebra.h"

#include "../src/rsamg/rsamg.h"

using namespace std;

namespace ug{
namespace bridge{
namespace registerRSAMG{

/**
 * Class exporting the functionality. All functionality that is to
 * be used in scripts or visualization must be registered here.
 */
struct Functionality
{

/**
 * Function called for the registration of Algebra dependent parts.
 * All Functions and Classes depending on Algebra
 * are to be placed here when registering. The method is called for all
 * available Algebra types, based on the current build options.
 *
 * @param reg				registry
 * @param parentGroup		group for sorting of functionality
 */
template <typename TAlgebra>
static void Algebra(Registry& reg, string grp)
{
	string suffix = GetAlgebraSuffix<TAlgebra>();
	string tag = GetAlgebraTag<TAlgebra>();

//	typedefs for this algebra
	typedef TAlgebra algebra_type;
	typedef typename algebra_type::vector_type vector_type;
	typedef typename algebra_type::matrix_type matrix_type;

	typedef RSAMG<algebra_type> T;

	reg.add_class_<T, AMGBase<algebra_type> > (string("RSAMGPreconditioner").append(suffix), grp, "Ruge-Stueben Algebraic Multigrid Preconditioner")
		.add_constructor()
		.add_method("set_epsilon_strong", &T::set_epsilon_strong, "", "epsilon_str", "sets epsilon_strong, a measure for strong connectivity")
		.add_method("get_epsilon_strong", &T::get_epsilon_strong, "epsilon_strong", "")
		.add_method("set_epsilon_strong_positive", &T::set_epsilon_strong, "", "epsilon_str", "sets epsilon_strong, a measure for strong connectivity")
		.add_method("get_epsilon_strong_positive", &T::get_epsilon_strong, "epsilon_strong", "")
		.add_method("set_prolongation_truncation", &T::set_prolongation_truncation, "", "prolongation_tr", "sets the truncation of interpolation")
		.add_method("get_prolongation_truncation", &T::get_prolongation_truncation, "prolongation_tr")

		.add_method("enable_aggressive_coarsening_A", &T::enable_aggressive_coarsening_A, "", "nrOfPaths", "enables aggressive coarsening (A1 or A2) on the first level.")
		.add_method("disable_aggressive_coarsening", &T::disable_aggressive_coarsening, "", "", "disables aggressive coarsening")
		.add_method("is_aggressive_coarsening", &T::is_aggressive_coarsening)
		.add_method("is_aggressive_coarsening_A", &T::is_aggressive_coarsening_A)
		.add_method("check_vector", &T::check_vector)

#ifdef UG_PARALLEL
		.add_method("set_parallel_coarsening", &T::set_parallel_coarsening)
#endif
		.set_construct_as_smart_pointer(true)
		;

	reg.add_class_to_group(string("RSAMGPreconditioner").append(suffix), "RSAMGPreconditioner", tag);

}
}; // end Functionality
}// end AMG

void RegisterBridge_RSAMG(Registry& reg, string grp)
{
	
	grp.append("/Algebra/Preconditioner");

#ifdef UG_PARALLEL
	reg.add_class_<	IParallelCoarsening > ("IParallelCoarsening", grp, "Parallel Coarsening Interface (RSAMG)");
	reg.add_function("GetFullSubdomainBlockingCoarsening", GetFullSubdomainBlockingCoarsening, grp);
	reg.add_function("GetColorCoarsening", GetColorCoarsening, grp);
	reg.add_function("GetRS3Coarsening", GetRS3Coarsening, grp);
	reg.add_function("GetCLJPCoarsening", GetCLJPCoarsening, grp);
	reg.add_function("GetFalgoutCoarsening", GetFalgoutCoarsening, grp);
	reg.add_function("GetMinimumSubdomainBlockingCoarsening", GetMinimumSubdomainBlockingCoarsening, grp);
	reg.add_function("GetCoarseGridClassificationCoarsening", GetCoarseGridClassificationCoarsening, grp);
	reg.add_function("GetSimpleParallelCoarsening", GetSimpleParallelCoarsening, grp);
#endif

	typedef registerRSAMG::Functionality Functionality;
#ifdef UG_CPU_1
	typedef boost::mpl::list<CPUAlgebra> AlgList;
	try{
		RegisterAlgebraDependent<Functionality, AlgList>(reg,grp);
	}
	UG_REGISTRY_CATCH_THROW(grp);
#endif
/*#ifdef UG_CRS_1
    typedef boost::mpl::list<CRSAlgebra> AlgList2;
	try{
		RegisterAlgebraDependent<Functionality, AlgList2>(reg,grp);
	}
	UG_REGISTRY_CATCH_THROW(grp);
#endif*/

}

} // namespace bridge
} // namespace ug
