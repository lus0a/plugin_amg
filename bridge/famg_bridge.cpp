/*
 * amg_bridge.cpp
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
#include "bridge/util_overloaded.h"

#include "lib_algebra/lib_algebra.h"
#include "named_infos.h"

#include "../src/famg/famg.h"

using namespace std;

namespace ug{
namespace bridge{
namespace registerFAMG{

template<typename TAlgebra, typename TClass>
ExportedClass<TClass> &addClassGroup(Registry& reg, std::string name, std::string grp)
{
	string suffix = GetAlgebraSuffix<TAlgebra>();
	string tag = GetAlgebraTag<TAlgebra>();
	ExportedClass<TClass>  &c = reg.add_class_<	TClass > (name+suffix, grp);
	reg.add_class_to_group(name+suffix, name, tag);
	return c;
}


template<typename TAlgebra, typename TClass, typename TBaseClass>
ExportedClass<TClass> &addClassGroup(Registry& reg, std::string name, std::string grp)
{
	string suffix = GetAlgebraSuffix<TAlgebra>();
	string tag = GetAlgebraTag<TAlgebra>();
	ExportedClass<TClass>  &c = reg.add_class_<	TClass, TBaseClass > (name+suffix, grp);
		reg.add_class_to_group(name+suffix, name, tag);
	return c;
}

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
		typedef FAMG<algebra_type> T;
		typedef AMGBase<algebra_type> TBase;

		reg.add_class_<T, TBase> (string("FAMGPreconditioner").append(suffix), grp, "Filtering Algebraic Multigrid")
			.add_constructor()
			.add_method("set_aggressive_coarsening", &T::set_aggressive_coarsening)
			.add_method("set_delta", OVERLOADED_METHOD_PTR(void, T, set_delta, (double delta)),
					"", "delta", "\"Interpolation quality\" F may not be worse than this (F < m_delta). e.g. 0.5")
			.add_method("set_delta", OVERLOADED_METHOD_PTR(void, T, set_delta, (double delta, size_t level)),
					"", "delta#level", "\"Interpolation quality\" F may not be worse than this (F < m_delta). e.g. 0.5")
			.add_method("get_delta", &T::get_delta, "delta")
			.add_method("set_theta", &T::set_theta, "" , "theta", "with multiple parents pairs, discard pairs with m_theta * F > min F. e.g. 0.9")
			.add_method("get_theta", &T::get_theta, "theta")


			.add_method("set_damping_for_smoother_in_interpolation_calculation",&T::set_damping_for_smoother_in_interpolation_calculation)

//			.add_method("add_testvector", OVERLOADED_METHOD_PTR(void, T, add_testvector, (vector_type& c, double weight)),"", "testVector#weight","adds a testvector with weight")
//			.add_method("add_testvector", OVERLOADED_METHOD_PTR(void, T, add_testvector, (IVectorWriter<vector_type> *vw, double weight)),"", "testVector#weight", "adds a testvector with weight by using the IVectorWriter interface")
			.add_method("set_write_testvectors", &T::set_write_testvectors, "", "bWrite", "if true, write testvectors to path specified in set_matrix_write_path")

			.add_(Method("set_testvector_smoother", &T::set_testvector_smoother)
								 .parameter("smoother") .tooltip("sets the smoother to smooth testvectors"))

			.add_(Method("set_testvector_smooths", &T::set_testvector_smooths)
					 .parameter("n") .tooltip("number of smoothing steps to smooth testvectors"))

			.add_method("set_testvectors", &T::set_testvectors, "", "tvc")
			.add_method("get_testvector_calculator", &T::get_testvector_calculator, "tvc")
			.add_method("add_watchnode", &T::add_watchnode)

//			.add_method("reset_testvectors", &T::reset_testvectors, "", "removes all added testvectors")

			.add_method("set_prolongation_truncation", &T::set_prolongation_truncation, "", "prolongation_tr", "sets prolongation_truncation, a parameter used for truncation of interpolation. use with care! (like 1e-5)")
			.add_method("get_prolongation_truncation", &T::get_prolongation_truncation, "prolongation_tr")

			.add_method("set_galerkin_truncation", &T::set_galerkin_truncation, "", "galerkin_tr", "sets galerkin truncation, a parameter used to truncate the galerkin product. use with care! (like 1e-9)")
			.add_method("get_galerkin_truncation", &T::get_galerkin_truncation, "galerkin_tr")

			.add_method("set_H_reduce_interpolation_nodes_parameter", &T::set_H_reduce_interpolation_nodes_parameter, "", "HreduceParameter", "we can restrict the number of parent nodes by looking at the entries of H(i,j) to prevent high fill in rates (e.g. 1e-3)")
			.add_method("get_H_reduce_interpolation_nodes_parameter", &T::get_H_reduce_interpolation_nodes_parameter, "HreduceParameter")

			.add_method("set_prereduce_A_parameter", &T::set_prereduce_A_parameter, "", "prereduceA", "by setting this != 0.0, we reduce the matrix A before using it to its strong connections. (e.g. 1e-3)")
			.add_method("get_prereduce_A_parameter", &T::get_prereduce_A_parameter, "prereduceA")

			.add_method("set_external_coarsening", &T::set_external_coarsening, "", "bExternalCoarsening", "You need to set_parallel_coarsening in parallel.")
			.add_method("set_strong_connection_external", &T::set_strong_connection_external, "", "epsilon", "set strong_connection value for coarsening (like set_epsilon_strong in RSAMG)")
			.add_method("get_strong_connection_external", &T::get_strong_connection_external)
			.add_method("set_strong_connection_submask", &T::set_strong_connection_submask, "", "mask", "May specify sub-matrix (< 0 : std. matrix)")
			.add_method("get_strong_connection_submask", &T::get_strong_connection_submask)

			.add_method("set_write_coarsening", &T::set_write_coarsening)
			.add_method("set_read_coarsening", &T::set_read_coarsening)

			.add_method("set_prevent_ff_connection_external", &T::set_prevent_ff_connection_external)
			.add_method("get_prevent_ff_connection_external", &T::get_prevent_ff_connection_external)
			.add_method("set_injection_restriction", &T::set_injection_restriction)
			.add_method("set_max_parents", &T::set_max_parents)
			.add_method("set_H_as_inverse", &T::set_H_as_inverse)
			.add_method("set_three_parents_penalty", &T::set_three_parents_penalty)
			.add_method("set_deficit_solve_tolerance", &T::set_deficit_solve_tolerance)
			.add_method("set_XC_new", &T::set_XC_new)
			.add_method("set_XC_LE_griebel", &T::set_XC_LE_griebel)
			.add_method("set_indirect_interpolation", &T::set_indirect_interpolation)
			.add_method("set_use_constant_if_failed", &T::set_use_constant_if_failed)
			.add_method("set_xc_calc_interpolation_routine", &T::set_xc_calc_interpolation_routine)
			.add_method("set_beam", &T::set_beam)
			.add_method("adaptive_filtering",
					OVERLOADED_METHOD_PTR(void, T, adaptive_filtering,
							(SmartPtr<ILinearOperator<vector_type, vector_type> > linOp, vector_type &x, SmartPtr<IConvergenceCheck<vector_type> > spConvCheck, double expectedRate) ) )
			.add_method("adaptive_filtering",
					OVERLOADED_METHOD_PTR(void, T, adaptive_filtering,
							(SmartPtr<ILinearOperator<vector_type, vector_type> > linOp, vector_type &x, SmartPtr<IConvergenceCheck<vector_type> > spConvCheck, double expectedRate,
									SmartPtr<ITestvectorCalculator<TAlgebra> > spTV) ) )




	#ifdef UG_PARALLEL
			.add_method("set_parallel_coarsening", &T::set_parallel_coarsening, "", "parallelCoarsening", "e.g. GetColorCoarsening()")
	#endif
			.add_method("set_use_precalculate", &T::set_use_precalculate, "", "bUsePrecalculate", "experimental way of coarsening. beta.")

			.add_method("set_write_f_values", &T::set_write_f_values)

			.add_method("check_testvector",	OVERLOADED_METHOD_PTR(bool, T, check_testvector, ()) )
			.add_method("check_testvector",	OVERLOADED_METHOD_PTR(bool, T, check_testvector, (size_t fromlevel)) , "", "fromlevel")
			.add_method("check_testvector",	OVERLOADED_METHOD_PTR(bool, T, check_testvector, (size_t, size_t)) , "", "fromlevel#tolevel")

			.add_method("set_N2_interpolation", &T::set_N2_interpolation)
			.add_method("set_full_block_jacobi", &T::set_full_block_jacobi)

			.set_construct_as_smart_pointer(true)

			;
	}
	reg.add_class_to_group(string("FAMGPreconditioner").append(suffix), "FAMGPreconditioner", tag);


	addClassGroup<TAlgebra, ITestvectorCalculator<TAlgebra> >(reg, "ITestvectorCalculator", grp)
			.add_method("config_string", &ITestvectorCalculator<TAlgebra>::config_string);

	addClassGroup<TAlgebra, TestvectorsLinearElasticity2d<TAlgebra>, ITestvectorCalculator<TAlgebra> >(reg, "TestvectorsLinearElasticity2d", grp)
			.add_constructor() .set_construct_as_smart_pointer(true);
	addClassGroup<TAlgebra, TestvectorsLinearElasticity3d<TAlgebra>, ITestvectorCalculator<TAlgebra> >(reg, "TestvectorsLinearElasticity3d", grp)
				.add_constructor() .set_construct_as_smart_pointer(true);
	addClassGroup<TAlgebra, TestvectorsConstant<TAlgebra>, ITestvectorCalculator<TAlgebra> >(reg, "TestvectorsConstant", grp)
				.add_constructor() .set_construct_as_smart_pointer(true);
	addClassGroup<TAlgebra, TestvectorsLinear3d<TAlgebra>, ITestvectorCalculator<TAlgebra> >(reg, "TestvectorsLinear3d", grp)
				.add_constructor() .set_construct_as_smart_pointer(true);
	addClassGroup<TAlgebra, TestvectorsLinear2d<TAlgebra>, ITestvectorCalculator<TAlgebra> >(reg, "TestvectorsLinear2d", grp)
					.add_constructor() .set_construct_as_smart_pointer(true);
	addClassGroup<TAlgebra, TestvectorsLinearComponent<TAlgebra>, ITestvectorCalculator<TAlgebra> >(reg, "TestvectorsLinearComponent", grp)
						.add_constructor() .set_construct_as_smart_pointer(true);
	addClassGroup<TAlgebra, TestvectorsOnlyRotation2d<TAlgebra>, ITestvectorCalculator<TAlgebra> >(reg, "TestvectorsOnlyRotation2d", grp)
						.add_constructor() .set_construct_as_smart_pointer(true);
	addClassGroup<TAlgebra, TestvectorsAuto<TAlgebra>, ITestvectorCalculator<TAlgebra> >(reg, "TestvectorsAuto", grp)
				.ADD_CONSTRUCTOR( (size_t numTestvectors, size_t numKeep, SmartPtr<ILinearIterator<vector_type> > smoother, size_t numSmooths) )
					("numTestvectors#keep#smoother#numSmooths")
				.set_construct_as_smart_pointer(true);

	addClassGroup<TAlgebra, TestvectorsConstantWobbled<TAlgebra>, ITestvectorCalculator<TAlgebra> >(reg, "TestvectorsConstantWobbled", grp)
					.ADD_CONSTRUCTOR( (double delta) )
						("delta")
					.set_construct_as_smart_pointer(true);

	{
		typedef TestvectorsExternal<TAlgebra> T;
		addClassGroup<TAlgebra, TestvectorsExternal<TAlgebra>, ITestvectorCalculator<TAlgebra> >(reg, "TestvectorsExternal", grp)
						.add_constructor()
						.ADD_CONSTRUCTOR( (vector_type &t1) ) ()
						.ADD_CONSTRUCTOR( (vector_type &t1, vector_type &t2) ) ()
						.ADD_CONSTRUCTOR( (vector_type &t1, vector_type &t2, vector_type &t3) ) ()
						.ADD_CONSTRUCTOR( (std::vector<vector_type *> t) ) ()
						.add_method("add", &T::add)
						.set_construct_as_smart_pointer(true);
	}

}


}; // end Functionality
}// end AMG

void RegisterBridge_FAMG(Registry& reg, string grp)
{
	grp.append("/Algebra/Preconditioner");
	typedef registerFAMG::Functionality Functionality;
	typedef boost::mpl::list<CPUAlgebra> AlgList;

	try{
		RegisterAlgebraDependent<Functionality>(reg,grp);
	}
	UG_REGISTRY_CATCH_THROW(grp);
}

} // namespace bridge
} // namespace ug
