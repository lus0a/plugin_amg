/*
 * rsamg_bridge.cpp
 *
 *  Created on: 23.03.2011
 *      Author: mrupp
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

#include "../src/semiamg/semi_amg.h"

using namespace std;

namespace ug{
namespace bridge{
namespace registerSemiAMG{

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
	{
		typedef SemiAMG<algebra_type> T;
		typedef AMGBase<algebra_type> TBase;
		reg.add_class_<T, TBase > (string("SemiAMGPreconditioner").append(suffix), grp, "Ruge-Stueben Algebraic Multigrid Preconditioner")
				.add_constructor()
				.add_method("set_virtual_nodes",&T::set_virtual_nodes, "", "flag", "if true, then we will construct virtual")
				//.add_method("set_tensor_axis", &T::set_tensor_axis, "", "axis|coord|level", "sets epsilon_strong, a measure for strong connectivity")
#ifdef UG_PARALLEL
		//.add_method("set_parallel_coarsening", &SemiAMG<algebra_type>::set_parallel_coarsening)
#endif
		.set_construct_as_smart_pointer(true)
		;

	reg.add_class_to_group(string("SemiAMGPreconditioner").append(suffix), "SemiAMGPreconditioner", tag);
	}
}
}; // end Functionality
}// end AMG

void RegisterBridge_SemiAMG(Registry& reg, string grp)
{
	
	grp.append("/Algebra/Preconditioner");

#ifdef UG_PARALLEL
	//reg.add_class_<	IParallelCoarsening > ("IParallelCoarsening", grp, "Parallel Coarsening Interface (RSAMG)");
	/*reg.add_function("GetFullSubdomainBlockingCoarsening", GetFullSubdomainBlockingCoarsening, grp);
	reg.add_function("GetColorCoarsening", GetColorCoarsening, grp);
	reg.add_function("GetRS3Coarsening", GetRS3Coarsening, grp);
	reg.add_function("GetCLJPCoarsening", GetCLJPCoarsening, grp);
	reg.add_function("GetFalgoutCoarsening", GetFalgoutCoarsening, grp);
	reg.add_function("GetMinimumSubdomainBlockingCoarsening", GetMinimumSubdomainBlockingCoarsening, grp);
	reg.add_function("GetCoarseGridClassificationCoarsening", GetCoarseGridClassificationCoarsening, grp);
	reg.add_function("GetSimpleParallelCoarsening", GetSimpleParallelCoarsening, grp);*/
#endif

	typedef registerSemiAMG::Functionality Functionality;
	//typedef boost::mpl::list<CPUAlgebra> AlgList;
	try{
		RegisterAlgebraDependent<Functionality>(reg,grp);
	}
	UG_REGISTRY_CATCH_THROW(grp);

}

} // namespace bridge
} // namespace ug
