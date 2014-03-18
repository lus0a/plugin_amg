// new
/**
ý * \file famg_interpolation_calculator.h
 *
 * \author Martin Rupp
 *
 * \date 24.11.10
 *
 * Goethe-Center for Scientific Computing 2010.
 */

#ifndef __H__LIB_ALGEBRA__FAMG_SOLVER__FAMG_INTERPOLATION_CALCULATOR_H__
#define __H__LIB_ALGEBRA__FAMG_SOLVER__FAMG_INTERPOLATION_CALCULATOR_H__

#include "famg.h"
#include "famg_nodeinfo.h"
#include "../rsamg/rsamg.h"
#include "lib_algebra/algebra_common/local_helper.h"
#include "lib_algebra/small_algebra/lapack/eigenvalue2.h"
#include "lib_algebra/small_algebra/additional_math.h"
#include "neighborhood.h"

#include "common/util/ostream_util.h"
#include "common/util/sort_util.h"
#include "lib_algebra/small_algebra/solve_deficit.h"
#include "lib_algebra/small_algebra/small_matrix/print.h"

#define FAMG_SCALING
#define FAMG_FULL_BLOCK
#define FAMG_USE_DIRCHLET_AS_PARENTS 1

#define FAMG_FULL_BLOCK

namespace ug
{

/// \addtogroup plugins_amg_famg
/// \{

#define SYS_AS_SYS 0
//#define SYS_AS_SYS 1
#ifdef UG_DEBUG
#define IFI(s) if(is_watch_node(i)) { s; }
#else
#define IFI(s)
#endif

#define UG_LOGI(s) IFI(UG_LOG(s));


/**
 * FAMGInterpolationCalculator
 *
 * This class calculates the interpolation weights. It should perhaps be merged with FAMGLevelCalculator
 * its only public member functions are
 * void get_possible_parent_pairs(size_t i, stdvector<neighborstruct2> &possible_neighbors, famg_nodes &rating)
 * void get_all_neighbors_interpolation(size_t i, TPMatrix &P, famg_nodes &rating)
 */
template<typename TAlgebra>
class FAMGInterpolationCalculator
{
private:
	typedef typename TAlgebra::vector_type vector_type;
	typedef typename TAlgebra::matrix_type matrix_type;
	typedef typename matrix_type::value_type value_type;
	typedef typename TAlgebra::matrix_type prolongation_matrix_type;

	typedef typename matrix_type::value_type smallmat_type;
	typedef typename vector_type::value_type smallvec_type;
	typedef typename prolongation_matrix_type::value_type smallmat_ptype;
	typedef prolongation_matrix_type TPMatrix;
	typedef matrix_type TMatrix;
	typedef vector_type TVector;

	FAMGInterpolationCalculator(
			const FAMGInterpolationCalculator<TAlgebra> &other);

	class two_indices : public FixedArray1<size_t, 2>
	{
		using FixedArray1<size_t, 2>::at;
	public:
		two_indices(size_t i0, size_t i1)
		{
			at(0) = i0;
			at(1) = i1;
		}
	};
	class three_indices : public FixedArray1<size_t, 3>
	{
		using FixedArray1<size_t, 3>::at;
	public:
		three_indices(size_t i0, size_t i1, size_t i2)
		{
			at(0) = i0;
			at(1) = i1;
			at(2) = i2;
		}
	};

	class index_range
	{
	public:
		index_range(size_t s) :
				m_size(s)
		{
		}
		inline size_t size() const
		{
			return m_size;
		}
		inline size_t operator[](size_t i) const
		{
			return i;
		}
	private:
		const size_t m_size;
	};

public:

	std::string m_matrixWritePath;


	FAMGInterpolationCalculator(const matrix_type &A_,
			const matrix_type &A_OL2_, const FAMG<TAlgebra> &famg,
			stdvector<vector_type> &testvectors, stdvector<double> &omega,
			stdvector<double> &_fvalues) :
			A(A_), A_OL2(A_OL2_), m_testvectors(testvectors), m_omega(omega), fvalues(_fvalues),
			m_amghelper(famg.m_amghelper),
			debugWatchNodes(famg.debugWatchNodes)
	{
		m_matrixWritePath = famg.m_writeMatrixPath;
		m_delta = famg.get_delta();
		m_theta = famg.get_theta();
		m_bHAsInverse = famg.m_bHAsInverse;
		m_damping =
				famg.get_damping_for_smoother_in_interpolation_calculation();
		m_dProlongationTruncation = famg.get_prolongation_truncation();
		m_iMaxParents = famg.m_iMaxParents;
		m_threeParentsPenalty = famg.m_threeParentsPenalty;

		m_dHReduceInterpolationNodesParameter = famg.get_H_reduce_interpolation_nodes_parameter();
		m_iHReduceMinParents = famg.get_H_reduce_min_parents();
		testvectorsExtern = (m_testvectors.size() > 0);
		m_bUseConstantIfFailed = famg.m_bUseConstantIfFailed;
		m_dDeficitSolveTolerance = famg.m_dDeficitSolveTolerance;
	}

	void init()
	{
		bvisited.resize(A_OL2.num_rows(), false);
	}

	void check_weights(TPMatrix &P, size_t i)
	{
		IF_DEBUG(DBG_FAMG_INTERPOLATION, 2)
		{
			for (typename TPMatrix::row_iterator it =
					P.begin_row(i); it != P.end_row(i); ++it)
			{
				for(size_t k=0; k<GetRows(it.value()); k++)
				{
					double d = BlockRef(it.value(), k, k);
					if (d < 0.01 || d > 1)
						UG_LOG(	"P(" << i << ", " << it.index() << ")[" << k << "," << k <<"] = " << d << "\n");
				}
			}
		}
	}

	/**
	 * it could be that this processor can not set node i coarse because
	 * a) it is master on another processor
	 * b) it is already fine
	 * @param rating
	 * @param i
	 * @return true if i could be parent
	 */
	bool asParent(FAMGNodes<TPMatrix> &rating, size_t i)
	{
		return rating.i_can_set_coarse(i)
#ifdef FAMG_USE_DIRCHLET_AS_PARENTS
				|| rating[i].is_dirichlet()
#endif
		;
	}

#include "local_eigenvalue_test_computation.h"
#if 1
	/**
	 * get possible interpolations for i with 2 interpolation parents.
	 * @param i
	 * @param rating used to decide if node can be parent
	 * @return
	 */
	double adjust(double &doo)
	{
		if(doo < 0.6 && doo > 0.4) doo = 0.5;
	}
	double get_possible_parent_pairs(size_t i, FAMGNodes<TPMatrix> &rating)
	{
		AMG_PROFILE_NEXT(FAMG_getPPP_minCalc);
		int i_min = -1;
		double f_min = 1e12;


		stdvector<neighborstruct_var<smallmat_ptype> > np;
		neighborstruct_var<smallmat_ptype> s;
		s.parents.resize(2);

		set_sizes(2);
		//size_t blockSize = GetRows(H(0,0));
		double total_f_min=1e12;

		IFI(PrintVector(onlyN1, "N1"))

//		IFI(UG_LOG(CPPString(H, "H")));
//		UG_LOG("\n\n----------------------- " << i << " -------------------------\n");

		size_t blockSize = GetRows(H(0,0));
		UG_LOGI("-----------------" << i << "---------------------------\n");
		for (size_t p0 = 0; p0 < onlyN1.size(); p0++)
		{
			size_t gp0 = onlyN1[p0];
			if (!asParent(rating, gp0))
				continue;
			for (size_t p1 = p0 + 1; p1 < onlyN1.size(); p1++)
			{
				size_t gp1 = onlyN1[p1];
				if (!asParent(rating, gp1))
					continue;

				UG_LOGI("parents " << gp0 << ", " << gp1 << "\n");
				if (solve_KKT_copy(i, two_indices(p0, p1)) != 0)
					continue;
				UG_LOGI("F = " << F << "\n");

				/*



				//check_one_interpolation(i, p0, p1);
				UG_DLOG(DBG_FAMG_INTERPOLATION, 2, "F: " << F << "\n");*/

				s.F = F;
				s.Fb = Fb;

				s.parents[0].value = q[0];
				s.parents[1].value = q[1];

				s.parents[0].from = gp0;
				s.parents[1].from = gp1;

				IF_DEBUG(DBG_FAMG_INTERPOLATION, 5)
				{
					UG_LOG("interpolating " << gp0 << " with " << q[0] << "\n");
					UG_LOG("interpolating " << gp1 << " with " << q[1] << "\n");
				}
				if(blockSize == 1)
				{
					if ((F < 10 || F < 10*m_delta) && F < total_f_min)
							total_f_min = F;
					if (s.F < f_min)
					{
						f_min = s.F;
						i_min = i_neighborpairs.size();
					}
					if (F > m_delta)
						continue;
					i_neighborpairs.push_back(s);
				}
				else
				{
					if(s.Fb[0] < m_delta || s.Fb[1] < m_delta || (blockSize == 2 || s.Fb[2] < m_delta))
						np.push_back(s);
				}
				//std::cout << s.parents[0].value << "   " << s.parents[1].value << " : " << s.F << "\n";



				// F = q^T H q
			}
		}

		if(blockSize == 2)
		{
			for(size_t k1=0; k1<np.size(); k1++)
			{
				for(size_t k2=0; k2<np.size(); k2++)
				{
	//				s.F = std::max(np[k1].Fb[0], np[k2].Fb[1]);
					//s.F = sqrt(np[k1].Fb[0]*np[k1].Fb[0], np[k2].Fb[1]*np[k2].Fb[1]);
					s.F = (np[k1].Fb[0] + np[k2].Fb[1])/2;

					if(np[k1].Fb[0] > m_delta || np[k2].Fb[1] > m_delta) continue;
					if(np[k1].Fb[0] == 0.0 || np[k2].Fb[1] == 0.0) continue;
					std::map<size_t, smallmat_type> m;
					s_interpolation<smallmat_ptype, size_t> &p0k1 = np[k1].parents[0];
					s_interpolation<smallmat_ptype, size_t> & p1k1 = np[k1].parents[1];
					s_interpolation<smallmat_ptype, size_t> &p0k2 = np[k2].parents[0];
					s_interpolation<smallmat_ptype, size_t> &p1k2 = np[k2].parents[1];
					m[p0k1.from] = 0.0;
					m[p1k1.from] = 0.0;
					m[p0k2.from] = 0.0;
					m[p1k2.from] = 0.0;
					BlockRef(m[p0k1.from], 0, 0) = BlockRef(p0k1.value, 0, 0);
					BlockRef(m[p1k1.from], 0, 0) = BlockRef(p1k1.value, 0, 0);
					BlockRef(m[p0k2.from], 1, 1) = BlockRef(p0k2.value, 1, 1);
					BlockRef(m[p1k2.from], 1, 1) = BlockRef(p1k2.value, 1, 1);

					s.parents.resize( m.size() );

					s.Fb[0] = np[k1].Fb[0];
					s.Fb[1] = np[k2].Fb[1];

					size_t k=0;
					for(typename std::map<size_t, smallmat_type>::iterator it = m.begin(); it != m.end(); ++it)
					{
//						if(BlockNorm(it->second) < 1e-3) continue; ///// !!!!!!
						s.parents[k].from = it->first;
						s.parents[k].value = it->second;
						k++;
					}
					s.parents.resize( k );

					if(m.size() != 2)
						s.F *= m_threeParentsPenalty;

					if (s.F < f_min)
					{
						f_min = s.F;
						i_min = i_neighborpairs.size();
					}
	//				s.print(rating);
					i_neighborpairs.push_back(s);
				}
			}
		}
		else if(blockSize == 3)
		{
			for(size_t k1=0; k1<np.size(); k1++)
			{
				for(size_t k2=0; k2<np.size(); k2++)
				{
					for(size_t k3=0; k3<np.size(); k3++)
					{
						s.F = (np[k1].Fb[0] + np[k2].Fb[1] + np[k3].Fb[2])/3;

						if(np[k1].Fb[0] > m_delta || np[k2].Fb[1] > m_delta || np[k3].Fb[2] > m_delta) continue;
						if(np[k1].Fb[0] == 0.0 || np[k2].Fb[1] == 0.0 || np[k3].Fb[1] == 0.0) continue;
						std::map<size_t, smallmat_type> m;
						s_interpolation<smallmat_ptype, size_t> &p0k1 = np[k1].parents[0];
						s_interpolation<smallmat_ptype, size_t> & p1k1 = np[k1].parents[1];
						s_interpolation<smallmat_ptype, size_t> &p0k2 = np[k2].parents[0];
						s_interpolation<smallmat_ptype, size_t> &p1k2 = np[k2].parents[1];
						s_interpolation<smallmat_ptype, size_t> &p0k3 = np[k3].parents[0];
						s_interpolation<smallmat_ptype, size_t> &p1k3 = np[k3].parents[1];
						m[p0k1.from] = 0.0;
						m[p1k1.from] = 0.0;
						m[p0k2.from] = 0.0;
						m[p1k2.from] = 0.0;
						m[p0k3.from] = 0.0;
						m[p1k3.from] = 0.0;
						BlockRef(m[p0k1.from], 0, 0) = BlockRef(p0k1.value, 0, 0);
						BlockRef(m[p1k1.from], 0, 0) = BlockRef(p1k1.value, 0, 0);
						BlockRef(m[p0k2.from], 1, 1) = BlockRef(p0k2.value, 1, 1);
						BlockRef(m[p1k2.from], 1, 1) = BlockRef(p1k2.value, 1, 1);
						BlockRef(m[p0k3.from], 2, 2) = BlockRef(p0k2.value, 2, 2);
						BlockRef(m[p1k3.from], 2, 2) = BlockRef(p1k2.value, 2, 2);

						s.parents.resize( m.size() );

						s.Fb[0] = np[k1].Fb[0];
						s.Fb[1] = np[k2].Fb[1];
						s.Fb[2] = np[k3].Fb[2];

						size_t k=0;
						for(typename std::map<size_t, smallmat_type>::iterator it = m.begin(); it != m.end(); ++it)
						{
							s.parents[k].from = it->first;
							s.parents[k].value = it->second;
							k++;
						}

						if(m.size() != 2)
							s.F *= m_threeParentsPenalty;

						if (s.F < f_min)
						{
							f_min = s.F;
							i_min = i_neighborpairs.size();
						}
		//				s.print(rating);
						i_neighborpairs.push_back(s);
					}
				}
			}
		}
		else if(blockSize == 1)
		{

		}
		else UG_THROW("blockSize = " << blockSize << " not supported");
		if(i_min != -1 && i_neighborpairs.size() > 0 && i_neighborpairs[0].F > i_neighborpairs[i_min].F)
			std::swap(i_neighborpairs[0],  i_neighborpairs[i_min]);
		IFI(UG_LOG("best 2-F is " << total_f_min << "\n"));
		//return total_f_min != 1e12 ? total_f_min : 0.0;
		return 0.0;
	}

	/**
	 * get possible interpolations for i with 3 interpolation parents.
	 * @param i
	 * @param rating used to decide if node can be parent
	 * @return
	 */
	double get_possible_parent_threes(size_t i, FAMGNodes<TPMatrix> &rating)
	{
		AMG_PROFILE_FUNC();
		AMG_PROFILE_NEXT(FAMG_getPPP_minCalc);
		int i_min = -1;
		double f_min = 1e12;

		neighborstruct_var<smallmat_ptype> s;
		s.parents.resize(3);
		set_sizes(3);
		//size_t blockSize = GetRows(H(0,0));
		double total_f_min=1e12;
		size_t gp[3];
		for (size_t p0 = 0; p0 < onlyN1.size(); p0++)
		{
			gp[0] = onlyN1[p0];
			if (!asParent(rating, gp[0]))
				continue;
			for (size_t p1 = p0 + 1; p1 < onlyN1.size(); p1++)
			{
				gp[1] = onlyN1[p1];
				if (!asParent(rating, gp[1]))
					continue;
				for (size_t p2 = p1 + 1; p2 < onlyN1.size(); p2++)
				{
					gp[2] = onlyN1[p2];

					UG_LOGI("parents " << gp[0] << ", " << gp[1] << ", " << gp[2] << "\n");
					if (solve_KKT_copy(i, three_indices(p0, p1, p2), true) != 0)
						continue;

					if ((F < 10 || F < 10*m_delta) && F < total_f_min)
						total_f_min = F;

					if (Fb[0] > m_delta || Fb[1] > m_delta || F > m_delta)
						continue;

					F *= m_threeParentsPenalty;
					//check_one_interpolation(i, p0, p1);
					UG_DLOG(DBG_FAMG_INTERPOLATION, 2, "F: " << F << "\n");

					s.F = F;
					s.Fb = Fb;


					size_t numP=0;
					s.parents.resize(3);
					for(size_t k=0; k<3; k++)
					{
						//if(BlockNorm(q[k]) < 1e-3) continue; ///// !!!!!!

						s.parents[numP].value = q[k];
						s.parents[numP].from = gp[k];
						numP++;
					}
					s.parents.resize(numP);
					/*IF_DEBUG(DBG_FAMG_INTERPOLATION, 5)

					IF_DEBUG(DBG_FAMG_INTERPOLATION, 5)
					{
						UG_LOG("interpolating " << gp0 << " with " << q[0] << "\n");
						UG_LOG("interpolating " << gp1 << " with " << q[1] << "\n");
					}
					if (s.F < f_min)
					{
						f_min = s.F;
						i_min = i_neighborpairs.size();
					}*/

					//std::cout << s.parents[0].value << "   " << s.parents[1].value << " : " << s.F << "\n";

					i_neighborpairs.push_back(s);
				}
				// F = q^T H q
			}
		}
		/*if(i_min != -1 && i_neighborpairs.size() > 0 && i_neighborpairs[0].F > i_neighborpairs[i_min].F)
			std::swap(i_neighborpairs[0],  i_neighborpairs[i_min]);
		IFI(UG_LOG("best 3-F is " << total_f_min << "\n"));
		return f_min;*/
		return 0.0;
	}
#endif

	/**
	 * prepare interpolation calculation
	 * that is: Get Neighborhood of i, calculate H-matrix.
	 * @param i
	 * @param rating
	 * @return
	 */
	bool prepare_interpolation_calculation(size_t i, FAMGNodes<TPMatrix> rating)
	{
		UG_DLOG(DBG_FAMG_INTERPOLATION, 2,
					"node " << i << " [" << rating.get_original_index(i) << "]\n");

		UG_ASSERT(rating[i].is_valid_rating(),
				"node " << i << " has no valid rating " << rating[i].get_rating());

		if (testvectorsExtern && BlockNorm2(m_testvectors[0][i]) < 1e-12)
		{
			UG_DLOG(DBG_FAMG_INTERPOLATION, 2,
					"testvector in " << i << " nearly 0, set as fine\n");
			rating.set_fine(i);
			return false;
		}

		// see local_eigenvalue_test_computation.h
		calc_special(i);

		GetNeighborhood(A_OL2, i, onlyN1, onlyN2, bvisited);
		if (onlyN1.size() == 0)
		{
			UG_DLOG(DBG_FAMG_INTERPOLATION, 2, "node has no connections, set as fine\n");
			rating.set_fine(i);
			return false;
		}

		// get H on onlyN1,i_index,onlyN2.
		get_H(i, rating);


		IF_DEBUG(DBG_FAMG_INTERPOLATION, 5)
		{
			for (size_t n = 0; n < onlyN1.size(); n++)
				UG_DLOG(DBG_FAMG_INTERPOLATION, 4,
						onlyN1[n] << " = " << rating.info(onlyN1[n]) << " ");
			UG_DLOG(DBG_FAMG_INTERPOLATION, 4, "\n");
		}
		return true;
	}

	size_t get_level()
	{
		return m_amghelper.positions.size()-1;
	}

	/** get_possible_parents: calculates possible interpolations for index i
	 * \param i
	 * \param possible_neighbors	here goes list of possible interpolations (2 or 3 parents)
	 * \param rating				fine/coarse infos of the nodes
	 */
	double get_possible_parents(size_t i,
			stdvector<neighborstruct_var<smallmat_ptype> > &possible_neighbors, FAMGNodes<TPMatrix> &rating)
	{

		UG_LOGI("----------  LEVEL " << get_level() << " node " << i << "---------------\n");
		if(false && (get_level() == 0 && (i == 1057 || i == 1056)) )
		{
			rating.set_uninterpolateable(i);
			return 0.0;
		}
		//if(i==0)
		{
//			UG_LOG("\n\n####################################\n\n");

//			UG_LOG("\n\n####################################\n\n");
		}

		AMG_PROFILE_FUNC();
//		UG_DLOG(DBG_FAMG_INTERPOLATION, 2, "\n\n\n\n============================\n\n\n");
		if(prepare_interpolation_calculation(i, rating) == false) return 0;
//		IFI(UG_LOG("----" << i << "\n"));
		//UG_LOG("\n\n");
		i_neighborpairs.clear();

		if(is_watch_node(i))
		{
		double bestF = get_all_neighbors_interpolation(i);
		}
//		IFI(UG_LOG("Best F is " << get_all_neighbors_interpolation(i) << "\n"));

		get_possible_parent_pairs(i, rating);
		double bestf = 1e3*m_delta;
		if(i_neighborpairs.size() > 0) bestf = i_neighborpairs[0].F;
		if(m_iMaxParents == 3 || (m_iMaxParents == 4 && bestf > 0.5*m_delta))
			get_possible_parent_threes(i, rating);

		/*if(i_neighborpairs.size() > 0)
		{
			bestf = i_neighborpairs[0].F;

			for(size_t i=0; i<i_neighborpairs.size(); i++)
			{
				UG_LOG("2-Fb " << i << ": " << i_neighborpairs[i]  << "\n");
			}
		}*/
		/*if(m_iMaxParents == 3 || (m_iMaxParents == 4 && bestf > 0.5*m_delta))
			get_possible_parent_threes(i, rating);*/
//		UG_LOG("\n\n");

		if(i_neighborpairs.size() == 0)
		{
			UG_DLOG(DBG_FAMG_INTERPOLATION, 2,
			std::endl << rating.get_original_index(i) << ": UNINTERPOLATEABLE" << std::endl);
			//UG_ASSERT(0, "node has no parents :'-(");
			UG_ASSERT(rating.i_must_assign(i), i);
			rating.set_uninterpolateable(i);
			return 0.0;
		}
		else
			return extract_best_parents(i, possible_neighbors);

	}

	/**
	 * from the list of possible parents (in i_neighborpairs), extract the best (with F=f_min), and
	 * the ones with F * m_theta < f_min. (e.g. m_theta = 0.9)
	 * @param i
	 * @param possible_neighbors
	 * @return
	 */
	double extract_best_parents(size_t i, stdvector<neighborstruct_var<smallmat_ptype> > &possible_neighbors)
	{
		//AMG_PROFILE_NEXT(FAMG_getPPP_cleanup);
		possible_neighbors.clear();

		// minimal element is always first in possible_neighbors[i] list.
		int imin = 0;
		double f_min = i_neighborpairs[0].F;


		for (size_t j = 0; j < i_neighborpairs.size(); j++)
			if (i_neighborpairs[j].F < f_min)
				f_min = i_neighborpairs[j].F;

		for (size_t j = 0; j < i_neighborpairs.size(); j++)
			if (m_theta * i_neighborpairs[j].F <= f_min)
			{
				possible_neighbors.push_back(i_neighborpairs[j]);
				imin = j;
			}
		if(imin != 0)
			std::swap(i_neighborpairs[imin], i_neighborpairs[0]);

		//#ifdef FAMG_PRINT_POSSIBLE_PARENTS
/*		IF_DEBUG(DBG_FAMG_INTERPOLATION, 2)
		{
			UG_LOG(std::endl << i << ": onlyN1.size = " << onlyN1.size() << " f_min = " << f_min << std::endl);
			UG_LOG(possible_neighbors.size() << " accepted pairs:" << std::endl);
			for (size_t j = 0; j < i_neighborpairs.size(); j++)
				if (m_theta * i_neighborpairs[j].F <= f_min)
					i_neighborpairs[j].print(rating);
			UG_LOG("not accepted pairs: " << std::endl);
			for (size_t j = 0; j < i_neighborpairs.size(); j++)
				if (m_theta * i_neighborpairs[j].F > f_min)
					i_neighborpairs[j].print(rating);
			UG_LOG(std::endl);
		}*/

		i_neighborpairs.clear();
		return f_min;
	}


	/**
	 * resizes needed matrices
	 * @param N
	 */
	void set_sizes(size_t N, int nrTV=-1)
	{
		q.resize(N);
		// assert blocksize is FIXED!!!
#if SYS_AS_SYS
		size_t blockSize = GetRows(H(0, 0));
		lq.resize(blockSize * N + 1);
		rhs.resize(blockSize * N + 1);
		KKT.resize(blockSize * N + 1, blockSize * N + 1);
#else
		size_t M = localTestvectors.size();
		if(nrTV == -1)
			nrTV = M;

		m_nrTestvectors = nrTV;
		lq.resize(N+m_nrTestvectors);
		rhs.resize(N+m_nrTestvectors);
		KKT.resize(N+m_nrTestvectors, N+m_nrTestvectors);
#endif
	}


	/**
	 * reduceKKTSystem : force values to zero for non-strong couplings because they would become very small
	 * the size of the KKT system stays the same.
	 * @sa markAllStrong, H_strong_measure
	 * @param indices the indices of the KKT system
	 * @param ntv number of testvectors
	 * @param b the block to be used in H_strong_measure
	 */
	template<typename TIndices>
	  void reduceKKTSystem(TIndices indices, size_t b)
	{
		const size_t N = indices.size();
		std::vector<bool> strong(N, false);
		markAllStrong(indices, strong, m_dHReduceInterpolationNodesParameter, m_iHReduceMinParents, b);

		for(size_t r=0; r<N; r++)
		{
			// fix entries to zero for non-strong
			if (strong[r]) continue;
			for(size_t c=0; c < KKT.num_rows(); c++)
				KKT(r, c)= 0.0;
			KKT(r, r)= 1.0;
			rhs[r] = 0.0;

		}
	}


	/**
	 * if Hii*aii is really small, the smoother is sufficient,
	 * so we don' t need to interpolate this unknown
	 * @param b	block
	 * @param i_index index of i in H
	 * @return true if smoother sufficient, false otherwise
	 */
	template<typename TIndices>
	bool check_smoother_sufficient(TIndices &indices, size_t b, size_t blockSize, size_t i_index)
	{
		Hii = BlockRef(H(i_index, i_index), b, b);
		const double Hiiscale = BlockRef(aii, b, b)*Hii;

		// quick check: do we really need to solve? (i.e., dirichlet component?)
		UG_DLOG(DBG_FAMG_INTERPOLATION, 5, "Hiiscale="<< Hiiscale << "("<< b<<")"<< std::endl);

		// if (Hii<-1e-10){ F+=1e+99; continue; }
		UG_ASSERT(Hii>-1e-10, "ERROR: negative values Hii="<< Hii <<" are impossible.");

		// skip, if entry can be treated by smoothing
		if (fabs(Hiiscale) < 0.01*m_delta/blockSize)
		{
			UG_LOG("Skipping comp="<<b<< " (Hiiscale="<< Hiiscale <<")"<< std::endl);
			//  no interpolation for this component
			for(size_t k=0; k<indices.size(); k++)
				BlockRef(q[k], b, b) = 0.0;
			F += Hiiscale;
			return true;
		}
		else
			return false;
	}

	/**
	 * assemble the KKT system from H and solve for a set of interpolation parents
	 * note that H has local indices 0 to onlyN1.size()-1
	 * their global indices are k -> onlyN1[k]
	 * indices is a subset of {0, ..., onlyN1.size()-1} of the local indices,
	 * so indices[k] has global index onlyN1[indices[k]].
	 * @param i
	 * @param indices set of parents.
	 * @param bCalcNegatives
	 * @return 0: ok, -1 : negative interpolation weights. -2 : could not solve KKT system
	 */
	template<typename TIndices>
	int solve_KKT(size_t i, const TIndices &indices, bool bCalcNegatives=false)
	{
		bCalcNegatives = false;
//		AMG_PROFILE_FUNC();

		const size_t blockSize=GetRows(H(0,0)); // assert blocksize is FIXED!!!
		const size_t N = indices.size();
		const size_t M = localTestvectors.size();
		size_t i_index = onlyN1.size();
		for(size_t k=0; k<indices.size(); k++)
			q[k] = 0.0;
		UG_DLOG(DBG_FAMG_INTERPOLATION, 3, "\n");
		F=0;

		Fb.resize(blockSize);
		bool bNegative=false;
		for(size_t b=0; b<blockSize; b++)
		{
			UG_LOGI("BLOCK [ " << b << " ]\n");

			if(check_smoother_sufficient(indices, b, blockSize, i_index) == true)
			{
				UG_LOGI("SMOOTHER SUFFICIENT\n");
				continue;

			}
			const double aii_bb  = BlockRef(aii, b, b);

			// fill KKT matrix
			KKT = 0.0;
			for(size_t r=0; r<N; r++)
				for(size_t c=0; c<N; c++)
					KKT(r, c)=aii_bb*BlockRef(H(indices[r], indices[c]), b, b);

			KKT(N, N) = 0;

			for(size_t k=0; k< m_nrTestvectors; k++)
				for(size_t j=0; j < N; j++)
					KKT(N+k, j) = KKT(j, N+k)= BlockRef(localTestvectors[k][indices[j]], b);


//			UG_LOG(JuliaString(KKT, "KKT"));
			for(size_t j=0; j<N; j++)
				rhs[j] = aii_bb*BlockRef(H(i_index, indices[j]), b, b);

			for(size_t k=0; k < m_nrTestvectors; k++)
				rhs[N+k] = BlockRef(localTestvectors[k][i_index],b);


//			IFI(rhs.maple_print("rhs"));
//			IF_DEBUG(DBG_FAMG_INTERPOLATION, 5) UG_LOG(JuliaString(KKT, "KKT1"));
//			IF_DEBUG(DBG_FAMG_INTERPOLATION, 5) rhs.maple_print("rhs1");
//			IFI(UG_LOG(JuliaString(KKT, "KKT1")));
//			IFI(rhs.maple_print("rhs1"));


			// reduce to essential entries
			// todo: do less work here (KKT system is too large)

//			if(indices.size() > 4) reduceKKTSystem(indices, b);
//			IF_DEBUG(DBG_FAMG_INTERPOLATION, 5) UG_LOG(JuliaString(KKT, "KKT"))
//			IF_DEBUG(DBG_FAMG_INTERPOLATION, 5) rhs.maple_print("rhs");

//			IFI(UG_LOG(JuliaString(KKT, "rKKT")));
//			IFI(rhs.maple_print("red_rhs"));

			// solve system
			DenseMatrix<VariableArray2<double> > KKTbefore = KKT;
			DenseVector<VariableArray1<double> > RHSbefore = rhs;
			if(SolveDeficit(KKT, lq, rhs, m_dDeficitSolveTolerance)==true)
			{
//				IFI(KKT.maple_print("KKT_After"));
//				IFI(rhs.maple_print("rhs_after"));
//				IFI(lq.maple_print("lq"));

				IF_DEBUG(DBG_FAMG_INTERPOLATION, 5) lq.maple_print("lq");

				bool bHasNegative=false;

				for(size_t k=0; k<indices.size(); k++)
				{
//					if(lq[k] < 1e-8 && lq[k] > -1e-8)
					if(lq[k] < 1e-5 && lq[k] > -1e-5)
						lq[k] = 0.0;
					if(lq[k] < 0)
					{
						IFI(UG_LOG("negative interpolation weights. aborting.\n"))
						IF_DEBUG(DBG_FAMG_INTERPOLATION, 2)
						{
							UG_LOG("WARNING: Negative entries at node " << i << ": ");
							for(size_t l=0; l<indices.size(); l++)
							{	UG_LOG("p" << l << " = " << onlyN1[ indices[l] ] << " "); }
							UG_LOG(" solve_KKT: lq[" << k << "] = " << lq[k] << "\n");
						}

						if(bCalcNegatives==false)
						{
							F+=1e12;
							Fb[b] = 1e12;
							//continue; //return -1;
						}
						bNegative=true;
						bHasNegative = true;
					}
					BlockRef(q[k], b, b) = lq[k];
				}
				if(bHasNegative) continue;
//				IF_DEBUG(DBG_FAMG_INTERPOLATION, 5) rhs.maple_print("rhs");
				//IF_DEBUG(DBG_FAMG_INTERPOLATION, 5)

				Hii = BlockRef(H(i_index, i_index), b, b);
				double lF=Hii;

				for(size_t k=0; k<N; k++)
					lF -= 2*lq[k]*BlockRef(H(i_index, indices[k]), b, b);
				for(size_t k1=0; k1<N; k1++)
					for(size_t k2=0; k2<N; k2++)
						lF += lq[k1]*BlockRef(H(indices[k1], indices[k2]), b, b)*lq[k2];
				lF *= BlockRef(aii, b, b);
				Fb[b] = lF;
				IFI(UG_LOGN("Fb[" << b << "] = " << lF << ")");
						for(size_t k1=0; k1<N; k1++) { UG_LOG("lq[0] = " << lq[0] << ", lq[1] = " << lq[1]); }
						UG_LOG("\n");
						)
				F += lF;
//				UG_ASSERT(lF>-1e-10, "ERROR: negative values lF="<<lF<<" (b="<<b <<") are impossible!");
			}
			else
			{
				Fb[b] = 1e12;
				F += 1e12;
				/*if(indices.size() >= 3)
				{
					UG_LOG("indices = " << indices.size() << " test = " << localTestvectors.size()
							<< " KKTNu = " << KKT.num_rows() << "\n");
					UG_LOG(JuliaString(KKTbefore, "KKT"));
					rhs.maple_print("rhs");
					UG_LOG(JuliaString(KKT, "KKT_after"));
					rhs.maple_print("rhs_after");
					UG_LOG(" solve_KKT: could not solve KKT (decomposition failed) at node " << i << "\n--------------------------------------------------\n\n");
				}*/
				IF_DEBUG(DBG_FAMG_INTERPOLATION, 2)
				{
				  UG_LOG(" solve_KKT: could not solve KKT (decomposition failed) at node " << i << "\n");
				}

//				return -2;
			}
		}
		IFI(UG_LOG("SUCCESS, F = " << F << "\n"));
		// TODO: Sinvoll?
		F /= blockSize;

		return (bCalcNegatives==false && bNegative == true) ? -1 : 0;
		//return 0;
	}

	// non-reference version
	template<typename TIndices>
	inline int solve_KKT_copy(size_t i, const TIndices indices, bool bCalcNegatives=false)
	{
		return solve_KKT(i, indices, bCalcNegatives);
	}

	/**
	 * check if interpolation from one parent is possible
	 * @param i node to be interpolated
	 * @param j parent
	 * @param newF interpolation measure F
	 * @param newQ interpolation (= w_{i,j} )
	 * @return
	 */
	double calc_one_interpolation(size_t i, size_t j, double &newF, smallmat_type &newQ)
	{
		newF = 0;
		size_t i_index = onlyN1.size();
		DenseVector<stdvector<smallvec_type> > &t0 = localTestvectors[0];
		newQ = 0.0;
		const size_t blockSize=GetRows(H(0,0)); // assert blocksize is FIXED!!!
		for(size_t b=0; b<blockSize; b++)
		{
			double w = BlockRef(t0[i_index], b) / BlockRef(t0[j], b);
			BlockRef(newQ, b, b) = w;
//			UG_LOG(BlockRef(H(i_index, i_index), b, b)  << " w " << w << " ij " <<  BlockRef(H(i_index, j), b, b) << "\n");
			newF += BlockRef(H(i_index, i_index), b, b)  - 2 * w *  BlockRef(H(i_index, j), b, b) + BlockRef(H(j, j), b, b) *w*w;
		}
		return newF;
	}

	void check_one_interpolation(size_t i, size_t p0, size_t p1)
	{
		double newF=0;
		smallmat_type newQ;
		const double limit=1;
		const double limit2=2;
		if(BlockRef(q[0], 0, 0) > limit*BlockRef(q[1], 0, 0))
		{
			calc_one_interpolation(i, p0, newF, newQ);
			if(newF < m_delta && newF < F*limit2)
			{
				q[0] = newQ;
				q[1] = 0.0;
				F = newF;
			}
		}
		else if(limit*BlockRef(q[0], 0, 0) < BlockRef(q[1], 0, 0))
		{
			calc_one_interpolation(i, p1, newF, newQ);
			if(newF < m_delta && newF < F*limit2)
			{
				q[1] = newQ;
				q[0] = 0.0;
				F = newF;
			}
		}

	}


	void printH(size_t i)
	{
		std::vector<size_t> N;
		for (size_t j = 0; j < onlyN1.size(); j++)
			N.push_back(onlyN1[j]);
		N.push_back(i);

		for (size_t r = 0; r < H.num_rows(); r++)
		{
			UG_LOG(N[r] <<": ");
			for (size_t c = 0; c < H.num_cols(); c++)
				UG_LOG("(" << N[c] << "): " << H(r, c) << " ");
			UG_LOG("\n");
		}
	}

	bool get_N2(size_t i, std::vector<size_t> &myN1)
	{
		AMG_PROFILE_FUNC();
		swap(onlyN1, myN1);

		for (size_t j = 0; j < bvisited.size(); j++)
			bvisited[j] = false;
		bvisited[i] = true;
		for (size_t j = 0; j < onlyN1.size(); j++)
			bvisited[onlyN1[j]] = true;

		onlyN2.clear();
		for (size_t j = 0; j < onlyN1.size(); j++)
		{
			UG_LOG("Neighbhors of " << onlyN1[j] << ": ");
			size_t k = onlyN2.size();
			GetNeighborhoodRec(A_OL2, onlyN1[j], onlyN2, bvisited);
			for (; k < onlyN2.size(); k++)
				UG_LOG(onlyN2[k] << " ");
		}

		CleanupFlagArray(bvisited, onlyN1);
		CleanupFlagArray(bvisited, onlyN2);
		bvisited[i] = false;
		return true;
	}

	template <typename TIndices>
	void maple_print(TIndices v, const char *title){
		std::cout << title << " = [";
		for(size_t i = 0; i<v.size(); ++i)
			std::cout << v[i] << ", ";
		std::cout << "]"<< std::endl;
	}

	/**
	 * the bigger, the stronger is the connection
	 *
	 * if j was the only node where we would interpolate from, then F would be
	 * F = H_{ii} - 2 H_{ij} w_j + H_{jj} w_j^2
	 * this is minimized for w_j = H_{ij}/H_{jj}
	 *
	 * @param j index in onlyN1.
	 * @param b block component
	 * @return H_{ij} / H_{jj}        (  (b,b) - components  )
	 */
	double H_strong_measure(int j, int b)
	{
		size_t i_index = onlyN1.size();
		return BlockRef(H(i_index, j), b,b) / BlockRef(H(j, j), b,b);
	}


	// select only coarse neighbors w/ sufficiently large couplings
	template<typename TIndices>
	void markAllStrong(const TIndices &coarseNeighbors,
			   std::vector<bool> &strongMarks, double eps, size_t nnodes, int b=0)
	{

		UG_DLOG(DBG_FAMG_INTERPOLATION, 3, "markAllStrong with sigma = "<<eps <<" for comp=" << b << std::endl);
		const size_t ncoarse= coarseNeighbors.size();
		strongMarks.resize(ncoarse);

		if(nnodes >= ncoarse)
		{
			for (size_t j = 0; j < ncoarse; j++) strongMarks[j] = true;
			return;
		}

		// determine maximum of offdiag/diag in i-th row
		double maxVal = 0.0;

		if (nnodes == 1)
		{
			for (size_t j = 0; j < ncoarse; j++)
			{
				const double val = H_strong_measure(coarseNeighbors[j], b);
				maxVal = std::max(maxVal, val);
				UG_DLOG(DBG_FAMG_INTERPOLATION, 3, val<< " ")
				strongMarks[j]=false;
			}
		}
		else
		{
			// determine maximum of offdiag/diag in i-th row (quick select)
			std::vector<double> sortVec(ncoarse);
			for (size_t j = 0; j <ncoarse; j++)
			{
				double val= H_strong_measure(coarseNeighbors[j], b);
				sortVec[j] =val;
				UG_DLOG(DBG_FAMG_INTERPOLATION, 3, val<< " ")
				strongMarks[j]=false;
			}
			std::sort(sortVec.begin(), sortVec.end());
			maxVal=sortVec[ncoarse-nnodes];
		}
		// eliminate coarse connections with small contributions
		maxVal *=eps;
		for (size_t j = 0; j < coarseNeighbors.size(); j++)
			if (H_strong_measure(coarseNeighbors[j], b) >= maxVal)
				strongMarks[j]=true;

		IF_DEBUG(DBG_FAMG_INTERPOLATION, 3) {
			maple_print(coarseNeighbors, "coarseNghbs");
			maple_print(strongMarks,     "strongMarks");
		}

	}


	// select only coarse neighbors w/ sufficiently large couplings
	template<typename TIndices>
	void markAllStrongAbs(const TIndices &coarseNeighbors,
				std::vector<bool> &strongMarks, double eps, int b=0)
	{

		UG_DLOG(DBG_FAMG_INTERPOLATION, 3, "markAllStrongAbs with sigma = "<<eps <<" for comp=" << b << std::endl);
		strongMarks.resize(coarseNeighbors.size());

		// determine maximum of offdiag/diag in i-th row
		double maxVal = 0.0;
		for (size_t j = 0; j < coarseNeighbors.size(); j++)
		{
			const double val = fabs( H_strong_measure(coarseNeighbors[j], b));
			maxVal = std::max(maxVal, val);
			UG_DLOG(DBG_FAMG_INTERPOLATION, 3, val<< " ")
			strongMarks[j]=false;
		}

		// eliminate coarse connections with small contributions
		maxVal *=eps;
		for (size_t j = 0; j < coarseNeighbors.size(); j++)
			if (fabs( H_strong_measure(coarseNeighbors[j], b)) >= maxVal)
				strongMarks[j]=true;

		IF_DEBUG(DBG_FAMG_INTERPOLATION, 3) {
			maple_print(coarseNeighbors, "coarseNghbs");
			maple_print(strongMarks,     "strongMarks");
		}

	}






	template<typename TIndices>
	std::vector<size_t> sortByStrength(const TIndices &coarseNeighbors, int b=-1)
	{
		std::vector<double> strongMeasure(coarseNeighbors.size());

		// determine maximum of offdiag/diag in i-th row
//		double maxVal = 0.0;
		size_t blockSize=GetRows(H(0,0));
		// invert with - because bigger is better -> sort is small to big.
		for (size_t j = 0; j < coarseNeighbors.size(); j++)
		{
			strongMeasure[j] = 0;
			if(b == -1)
				for(size_t bi=0; bi<blockSize; bi++)
					strongMeasure[j] -= fabs(H_strong_measure(coarseNeighbors[j], bi));
			else
				strongMeasure[j] = -fabs(H_strong_measure(coarseNeighbors[j], b));
		}

		/*UG_LOG("Strong Measures:\n");
		for(size_t i=0; i<strongMeasure.size(); i++)
			UG_LOG(onlyN1[coarseNeighbors[i]] << " = " << strongMeasure[i] << "\n")*/
		return GetSortedIndices(strongMeasure);
	}

	std::vector<size_t> sortByStrength(int b=0)
	{
		index_range r(onlyN1.size());
		return sortByStrength(r, b);
	}

	void getAllStrong(std::vector<size_t> &res, double tol, int nnodes=1)
	{
		std::vector<size_t> indices;
		for(size_t i=0; i<onlyN1.size(); i++) indices.push_back(i);
		std::vector<bool> strongMarksAll(indices.size(), false);
		std::vector<bool> strongMarks(indices.size(), false);
		size_t blockSize=GetRows(H(0,0));
		for(size_t b=0; b<blockSize; b++)
		{
		  markAllStrong(indices, strongMarks, tol, nnodes, b);
		  for(size_t i=0; i<strongMarks.size(); i++)
		    if(strongMarksAll[i]==false) strongMarksAll[i] = strongMarks[i];
		}
		res.clear();
		for(size_t i=0; i<strongMarks.size(); i++)
			if(strongMarks[i]) res.push_back(i);

	}


	// select only coarse neighbors w/ sufficiently large couplings
	void reduceToStrong(const std::vector<size_t> &coarseNeighbors,
			std::vector<size_t> &strongNeighbors, double eps, int b=0)
	{

		UG_DLOG(DBG_FAMG_INTERPOLATION, 3, "reduceToStrong with sigma = "<<eps <<" for comp=" << b << std::endl);
		IF_DEBUG(DBG_FAMG_INTERPOLATION, 3) {
			maple_print(coarseNeighbors, "coarse");
		}

		// determine maximum of offdiag/diag in i-th row
		double maxVal = 0.0;
		for (size_t j = 0; j < coarseNeighbors.size(); j++)
		{
			const double val = H_strong_measure(coarseNeighbors[j], b);
			maxVal = std::max(maxVal, val);
			UG_DLOG(DBG_FAMG_INTERPOLATION, 3, val<< " ")
		}

		// eliminate coarse connections with small contributions
		maxVal *=eps;
		for (size_t j = 0; j < coarseNeighbors.size(); j++)
			if (fabs(H_strong_measure(coarseNeighbors[j], b)) >= maxVal )
				strongNeighbors.push_back(coarseNeighbors[j]);


		IF_DEBUG(DBG_FAMG_INTERPOLATION, 3) {
			maple_print(strongNeighbors, "strong");
		}

	}

	// select only coarse neighbors w/ sufficiently large couplings
		void reduceToStrongNew(const std::vector<size_t> &coarseNeighbors,
				std::vector<size_t> &strongNeighbors, double eps)
		{
			// determine maximum among coarse neighbors
			double maxVal = BlockNorm(H(coarseNeighbors[0],coarseNeighbors[0]));
			for (size_t j = 1; j < coarseNeighbors.size(); j++)
				maxVal = std::max(BlockNorm(H(coarseNeighbors[j], coarseNeighbors[j])), maxVal);

			// select strongly coupled coarse neighbors
			for (size_t j = 0; j < coarseNeighbors.size(); j++)
				if (BlockNorm(H(coarseNeighbors[j], coarseNeighbors[j])) >= eps * maxVal)
					strongNeighbors.push_back(coarseNeighbors[j]);
		}

	size_t csolve_KKT_non_neg;
	int solve_KKT_non_neg(size_t i, std::vector<size_t> &indices, int tvSize, bool bCalcNegatives=false)
	{
		csolve_KKT_non_neg++;
		if(csolve_KKT_non_neg > 20) return -1;
		set_sizes(indices.size(), tvSize);
		IFI(UG_LOG("tvSize = " << tvSize << " (");		for(size_t j=0; j<indices.size(); j++)		{ UG_LOG(onlyN1[indices[j]] << " "); }		UG_LOG(")\n");
		)
		int ret = solve_KKT(i, indices);
		UG_LOGI(" ret = " << ret << ": " << Fb[0] << " " << Fb[1] << "\n");
		if(ret != -1)
			return ret;
		if(indices.size() <= 3) return -1;
		for(size_t j=0; j<indices.size(); j++)
		{
			std::vector<size_t> mI = indices;
			mI.erase(mI.begin()+j);
			ret = solve_KKT_non_neg(i, mI, tvSize);
			if(ret == 0 && Fb[0] < m_delta && Fb[1] < m_delta)
			{
				indices = mI;
				return ret;
			}
		}
		return -1;
	}

	bool direct_interpolation(size_t i, TPMatrix &P,
			FAMGNodes<TPMatrix> &rating, std::vector<size_t> &coarseNeighbors)
	{
		UG_LOGI("direct interpolation " << i << "\n");
		AMG_PROFILE_FUNC();

		get_H(i, rating); // todo: not really necessary if coarseNeighbors.size() == 1
/*
		if(coarseNeighbors.size() == 1)
		{
			// single parent
			DenseVector<stdvector<smallvec_type> > &t0 = localTestvectors[0];
			get_H(i, rating);
			P(i, onlyN1[coarseNeighbors[0]]) = 0.0;
			const size_t blockSize=GetSize(t0[onlyN1.size()]);
			for(size_t b=0; b<blockSize; b++)
			{

				BlockRef(P(i, onlyN1[coarseNeighbors[0]]), b, b) =
						BlockRef(t0[onlyN1.size()], b) / BlockRef(t0[coarseNeighbors[0]], b);
			}
			rating.set_fine(i);
			return true;
		}
*/
		std::vector<size_t> myNeighbors = coarseNeighbors;
		csolve_KKT_non_neg=0;
		if((solve_KKT_non_neg(i, myNeighbors, -1) == 0  && F < m_delta ) ||
				(m_bUseConstantIfFailed && solve_KKT_non_neg(i, myNeighbors, 1) == 0))
		{
//			if(csolve_KKT_non_neg != 1)UG_LOG("csolve_KKT_non_neg = " << csolve_KKT_non_neg << "\n");
			UG_LOGI("csolve_KKT_non_neg = " << csolve_KKT_non_neg << "\n");
			UG_COND_THROW(F<-1e-10, "ERROR: negative F values are impossible (i.e., find bug and fix it!)");

			if(F > m_delta)
			{

				UG_DLOG(DBG_FAMG_INTERPOLATION, 3, "direct_interpolation, had to set node " << i << " coarse: " << F << ">" << m_delta);
				rating.external_set_coarse(i);
			}
			else
			{
				UG_DLOG(DBG_FAMG_INTERPOLATION, 3, "coarse neighbors, Interpolating from ");
				size_t N = myNeighbors.size();
				if(m_dProlongationTruncation > 0.0)
				{
					UG_THROW("disabled due to systems");
				}
				if(m_dProlongationTruncation > 0.0)
				 {
					UG_THROW("no thanks");
					double dmax = BlockNorm(q[0]);
					for(size_t j=1; j<N; j++)
						dmax = std::max(BlockNorm(q[j]), dmax);
					for(size_t j=0; j<N; j++)
					{
						if(BlockNorm(q[j]) < dmax*m_dProlongationTruncation) continue;
						size_t node = onlyN1[myNeighbors[j]];
						UG_DLOG(DBG_FAMG_INTERPOLATION, 3, node << ": " << q[j] << ", ");
						P(i, node) = q[j];
					}
				 }
				else
				{
					for(size_t j=0; j<N; j++)
					{
						size_t nodej = onlyN1[myNeighbors[j]];
						if(!rating[i].is_coarse()) rating.set_coarse(nodej);
						P(i, nodej) = q[j];
						UG_DLOG(DBG_FAMG_INTERPOLATION, 3, nodej << ": " << q[j] << ", ");
					}
				}
				check_weights(P, i);
				if(fvalues.size() > i) fvalues[i] = F;
//				UG_LOG("setting " << i << " fine.");
				rating.set_fine(i);
			}
		}
		else
		{
//			if(csolve_KKT_non_neg != 1)UG_LOG("csolve_KKT_non_neg = " << csolve_KKT_non_neg << "\n");
			UG_LOGI("csolve_KKT_non_neg = " << csolve_KKT_non_neg << "\n");
			//rating.set_coarse(i);
			rating.external_set_coarse(i);
			UG_DLOG(DBG_FAMG_INTERPOLATION, 5, "direct_interpolation of node " << i << ": could not invert KKT system (coarse neighbors).\n");
		}
		return true;
	}

	bool indirect_interpolation(size_t i, TPMatrix &P,
			FAMGNodes<TPMatrix> &rating,
			std::vector<size_t> &interpolateNeighbors)
	{
		AMG_PROFILE_FUNC();
		get_H(i, rating);

		F = 0;
		set_sizes(interpolateNeighbors.size());
		if (solve_KKT(i, interpolateNeighbors) == 0 && F < m_delta)
		{
			//IF_DEBUG(DBG_FAMG_INTERPOLATION, 5)
				//q.maple_print("q");
			typedef typename std::map<size_t, smallmat_type> smallmat_map;
			smallmat_map localP;
			for (size_t j = 0; j < interpolateNeighbors.size(); j++)
			{
				size_t node = onlyN1[interpolateNeighbors[j]];
				for (typename TPMatrix::row_iterator it = P.begin_row(node); it != P.end_row(node); ++it)
					localP[it.index()] += it.value() * q[j];
			}

			if(m_dProlongationTruncation > 0.0)
			{
				double dmax = 0;
				for (typename smallmat_map::iterator it = localP.begin(); it != localP.end(); ++it)
					dmax = std::max(BlockNorm((*it).second), dmax);
				for (typename smallmat_map::iterator it = localP.begin(); it != localP.end(); ++it)
					if (BlockNorm((*it).second) > dmax * m_dProlongationTruncation)
						P(i, (*it).first) = (*it).second;
			}
			else
			{
				for (typename smallmat_map::iterator it = localP.begin(); it != localP.end(); ++it)
					P(i, (*it).first) = (*it).second;
			}

			check_weights(P, i);
			if (fvalues.size() > i)
				fvalues[i] = F;
			rating.set_aggressive_fine(i);
		}
		else
		{
			rating.external_set_coarse(i);
			UG_DLOG(DBG_FAMG_INTERPOLATION, 3, "had to set node " << i << " coarse: ");
			if (F == 0)
			{
				UG_DLOG(DBG_FAMG_INTERPOLATION, 3, "could not invert KKT system.\n");
			}
			else
			{
				UG_DLOG(DBG_FAMG_INTERPOLATION, 3,
						"F = " << F << " < delta = " << m_delta << "\n");
			}
		}
		return true;
	}

	void get_strong_connected_neighborhood(const TMatrix &A,
			const FAMGNodes<TPMatrix> &rating, size_t i,
			stdvector<size_t> &onlyN1, stdvector<size_t> &onlyN2,
			stdvector<size_t> &coarseNeighbors, stdvector<bool> &bvisited)
	{
		double minConnValue, maxConnValue, diag;
		GetNeighborValues(A, i, minConnValue, maxConnValue, diag);
		onlyN1.clear();
		onlyN2.clear();

		bvisited[i] = true;
		const double theta = 0.3;
		double barrier = theta * std::max(maxConnValue, dabs(minConnValue));
		for (typename TMatrix::const_row_iterator conn = A.begin_row(i);
				conn != A.end_row(i); ++conn)
		{
			if (dabs(amg_offdiag_value(conn.value())) < barrier) // only strong
				continue;
			if (conn.index() == i)
				continue; // skip diagonal
			onlyN1.push_back(conn.index());
			bvisited[conn.index()] = true;
			if (!rating[conn.index()].is_coarse())
				continue; // only coarse

			coarseNeighbors.push_back(onlyN1.size() - 1);
		}

		for (size_t i = 0; i < onlyN1.size(); i++)
		{
			GetNeighborValues(A, onlyN1[i], minConnValue, maxConnValue, diag);
			barrier = theta * std::max(maxConnValue, dabs(minConnValue));
			for (typename TMatrix::const_row_iterator conn = A.begin_row(
					onlyN1[i]); conn != A.end_row(onlyN1[i]); ++conn)
			{
				if (dabs(amg_offdiag_value(conn.value())) < barrier) // only strong
					continue;
				if (bvisited[conn.index()] == false)
				{
					onlyN2.push_back(conn.index());
					bvisited[conn.index()] = true;
				}
			}
			// GetNeighborhoodRec(A, onlyN1[i], onlyN2, bvisited);
		}

		CleanupFlagArray(bvisited, onlyN1);
		CleanupFlagArray(bvisited, onlyN2);
	}

	void get_strong_connected_neighborhood2(const TMatrix &A,
			const FAMGNodes<TPMatrix> &rating, size_t i, stdvector<size_t> &onlyN1,
			stdvector<size_t> &onlyN2, stdvector<size_t> &coarseNeighbors,
			stdvector<bool> &bvisited)
	{

//		bvisited.resize(0);
//		bvisited.resize(A.num_rows(), false);
		GetNeighborhood(A_OL2, i, onlyN1, onlyN2, bvisited);

		for (size_t j = 0; j < onlyN1.size(); j++)
		{
			if (!rating[onlyN1[j]].is_fine())
				coarseNeighbors.push_back(j);
		}
	}

	double get_strong_interpolation(double bestF, size_t i)
	{
		std::vector<size_t> strong = sortByStrength();
		UG_LOG("sorted strong: ");
		for(size_t j=0; j<strong.size(); j++)
		{UG_LOG(onlyN1[strong[j]] << ", ");}
			UG_LOG("\n");


	}

	double get_all_neighbors_interpolation(size_t i)
	{
		std::vector<size_t> strong = sortByStrength();
//		UG_LOG("sorted strong: ");
//		for(size_t j=0; j<strong.size(); j++)
//		{UG_LOG(onlyN1[strong[j]] << ", ");}
//		UG_LOG("\n");

		index_range ir(onlyN1.size());
		set_sizes(onlyN1.size());
		int ret = solve_KKT(i, ir, false);

		neighborstruct_var<smallmat_type> s;
		s.F = F;
		s.Fb = Fb;

		s.parents.resize(q.size());
		for(size_t i=0; i<q.size(); i++)
		{
			s.parents[i].value = q[i];
			s.parents[i].from = onlyN1[i];
			UG_LOG(onlyN1[i] << " = " << q[i] << "\n");
		}

		//UG_LOG(i << " ret = " << ret << " s = " << s << "\n")
		if(ret == 0 || ret == -1)
			return F;
		else if(ret == -2)
			return -1;
		else
			return -1;

		/*UG_LOG("get_all_neighbors_interpolation(" << i << ")\n");
		// H already calculated
		std::vector<size_t> strong = sortByStrength();
		UG_LOG("sorted strong: ");
		for(size_t j=0; j<strong.size(); j++)
		{UG_LOG(onlyN1[strong[j]] << ", ");}
									UG_LOG("\n");
		while(strong.size() > 1)
		{
			set_sizes(strong.size());
			UG_LOG("trying ");
			for(size_t j=0; j<strong.size(); j++)
					{UG_LOG(onlyN1[strong[j]] << ", ");}
												UG_LOG("\n");
			int ret = solve_KKT(i, strong, true);
			if(ret == 0 || ret == -1)
				return F;
			else if(ret == -2)
				return 1000;
			strong.resize(strong.size()-1);
		}
		return 1000;*/
	}

	// get_all_neighbors_interpolation: RSAMG strong
	//---------------------------------------
	/** calculates an interpolation of the node i, when i is interpolating by all his neighbors
	 * if neighbors are fine, their interpolation form coarse nodes is used (indirect)
	 * \param i					node index for which to calculate possible parent pairs
	 * \param P					matrix for the interpolation
	 * \param rating			fine/coarse infos of the nodes
	 * \param indirectInterpolationLevel	if true, only use coarse neighbors
	 */
	bool get_all_neighbors_interpolation(size_t i, TPMatrix &P,
			FAMGNodes<TPMatrix> &rating, int indirectInterpolationLevel = -1)
	{
		//UG_LOG("get_all_neighbors_interpolation " << i << "\n");
		AMG_PROFILE_FUNC();
//		UG_DLOG(DBG_FAMG_INTERPOLATION, 3,
		UG_LOGI("get_all_neighbors_interpolation on node " << rating.get_original_index(i) << "\n")

		stdvector<size_t> coarseNeighbors;


		get_strong_connected_neighborhood2(A_OL2, rating, i, onlyN1, onlyN2,
				coarseNeighbors, bvisited);
		IFI(PrintVector(onlyN1, "N1"))
		IFI(PrintVector(coarseNeighbors, "cN1"))
//		UG_LOG("coarseNeighbors: " << coarseNeighbors.size() << "\n");
		if (onlyN1.size() == 0 || coarseNeighbors.size() == 0)
		{
			UG_DLOG(DBG_FAMG_INTERPOLATION, 3, "no neighbors -> F");
//			UG_LOG("setting1417 " << i << " coarse.");
			rating.set_dirichlet(i);
			return true;
		}

//		UG_LOG("di\n");
		//const double &aii = A_OL2(i,i);

		/*if (coarseNeighbors.size() < 2)
			return indirect_interpolation(i, P, rating, coarseNeighbors);
		else if (indirectInterpolationLevel >= 0)
		{
			UG_DLOG(DBG_FAMG_INTERPOLATION, 3,
					"indirectInterpolationLevel = " << indirectInterpolationLevel << ", but no coarse neighbors of this level");
			return false;
		}
		else*/
		{
			/*std::vector<size_t> innerNodes;
			for (size_t j = 0; j < onlyN1.size(); j++)
				if (rating.is_inner_node(onlyN1[j])	|| rating.is_master(onlyN1[j]))
					innerNodes.push_back(j);*/
			return direct_interpolation(i, P, rating, coarseNeighbors);
		}
	}

private:
	/**
	 * get_H
	 * \param i 		central node i
	 * \param rating	only used for debug output (get_orinal_index)
	 *
	 * this function needs onlyN1 and onlyN2 set
	 */
	bool get_H(size_t i, FAMGNodes<TPMatrix> &rating)
	{
		AMG_PROFILE_FUNC();
		// replace this with
		// stdvector<stdvector<size_t> > neighbors(3);
		/// stdvector<size_t> &onlyN1 = neighbors[1];
		/// stdvector<size_t> &onlyN2 = neighbors[2];
		// GetNeighborhoodHierachy(A, i, neighbors, bvisited);

		// 1. Get Neighborhood N1 and N2 of i.
		//AMG_PROFILE_BEGIN(AMG_H_GetNeighborhood);

		PRINT_VECTOR_AMGDEBUG(5, onlyN1, "\nn1-neighbors");
		PRINT_VECTOR_AMGDEBUG(5, onlyN2, "\nn2-neighbors");

		//AMG_PROFILE_NEXT(AMG_H_CreateN2);
		// create N2: N2 consists of all 2-neighbors of N in this order:
		// onlyN1, i, onlyN2
		N2 = onlyN1;
		N2.push_back(i);
		N2.insert(N2.end(), onlyN2.begin(), onlyN2.end());

		PRINT_VECTOR_AMGDEBUG(2, onlyN1, "\nonlyN1");
		PRINT_VECTOR_AMGDEBUG(2, onlyN2, "\nonlyN2");
		PRINT_VECTOR_AMGDEBUG(2, N2, "\nN2");

		// 2. get submatrix in A_OL2 on N2
		//AMG_PROFILE_NEXT(AMG_H_GetLocalMatrix);

		S.resize(N2.size(), N2.size());
		S = 0.0;

		GetLocalMatrix(A_OL2, S, &N2[0], &N2[0]);

		size_t i_index = onlyN1.size();
		GetDiagSqrt(aiiSqrt, S(i_index, i_index));
		GetDiag(aii, S(i_index, i_index));

		// make S symmetric
		/*for(size_t r=0; r<S.num_rows(); r++)
		 for(size_t c=0; c<S.num_cols(); c++)
		 {
		 S(r, c) = (S(r,c)+S(c,r))/2
		 }*/


#ifdef FAMG_SCALING1
		// do scaling stuff // mr 01.08.2012
		for(size_t r=0; r<S.num_rows(); r++)
		{
			smallmat_type m = S(r,r);
			for(size_t c=0; c<S.num_cols(); c++)
				S(r, c) /= m;
		}
#endif

		IF_DEBUG(DBG_FAMG_INTERPOLATION, 5)
			S.maple_print("\nsubA");

		//AMG_PROFILE_END();
		// 3. calculate H from submatrix A
		calculate_H_from_local_A(i);

		// calculate testvector
		calculate_testvectors(i);




#if SYS_AS_SYS
		DenseVector<stdvector<smallvec_type> > &t0 = localTestvectors[0];
		// b = sqrt(Aii) _ {i}
		// Hii = <b, Hii b>
		Hii = EnergyProd(aiiSqrt, H(i_index, i_index), aiiSqrt);
		t_i = VecProd(t0[i_index], aiiSqrt);
#endif

		return true;
	}

	void calculate_Dinv_from_local_A(size_t i)
	{
		const size_t N = S.num_rows(); // = N2
		
		if (Dinv.size() != N)
			Dinv.resize(N);
		IF_DEBUG(DBG_FAMG_INTERPOLATION, 1)
		{
			for (size_t j = 0; j < N; j++)
				if (BlockNorm(S(j, j)) < 1e-12)
					UG_LOG(i << " j=" << j << " S(j,j)=" << S(j,j) << "\n");
		}

		for (size_t j = 0; j < N; j++)
		{
#ifdef FAMG_FULL_BLOCK

			Dinv[j] = S(j,j);
			Invert(Dinv[j]);
#else
			Dinv[j] = 0.0;
			for (size_t k = 0; k < GetRows(S(j, j)); k++)
				BlockRef(Dinv[j], k, k) = 1 / BlockRef(S(j, j), k, k);

#endif
		}
		PRINT_VECTOR_AMGDEBUG(5, Dinv, "Dinv");
	}

	void calculate_SF_from_Dinv_and_local_A(size_t i)
	{
		const size_t N = S.num_rows(); // = N2
		const size_t i_index = onlyN1.size();

		if (SF.size() != N)
			SF.resize(N);

		// get SF = 1-wDF^{-1} A  (F-smoothing)
		//AMG_PROFILE_BEGIN(AMG_HA_calculate_SF);

		// bei f-smoothing nie und nimmer damping (arne 3.juni)
		if (BlockNorm(S(i_index, i_index)) < 1e-12)
			UG_LOG("i=" << i << "S(i,i) = " << S(i,i) << "\n");

		// todo: only go to onlyN1.size(), since S(i_index, j) = 0 forall j > onlyN1.size()
		// SF[i] = (1 -wD^{-1} A )_( i, N2[i])
		for (size_t j = 0; j < N; j++)
			SF[j] = -Dinv[i_index] * S(i_index, j);
		SF[i_index] += 1.0;

		PRINT_VECTOR_AMGDEBUG(5, SF, "SF");
	}

	void calculate_Sjac_from_Dinv_and_local_A()
	{
		const size_t N = S.num_rows(); // = N2
		// get S = 1-w D^{-1} A
		for (size_t r = 0; r < N; r++)
		{
			// todo: only go to onlyN1.size(), since S(i_index, j) = 0 forall j > onlyN1.size()
			for (size_t c = 0; c < N; c++)
			{
				S(r, c) = -m_damping * Dinv[r] * S(r, c);
				if (r == c)
					S(r, c) += 1.0;
			}
		}

		IF_DEBUG(DBG_FAMG_INTERPOLATION, 5)
			S.maple_print("Sjac");
	}

	void calculate_S_from_SF_and_Sjac()
	{
		const size_t N = S.num_rows(); // = N2
		const size_t i_index = onlyN1.size();

		// get S' = SF S
		// (possible without temporary since SF is mostly Id,
		//  and then SF S is addition of rows of S)

		// S'(i_index, i) = sum_t SF(i_index,t) * S(t, i)
		for (size_t j = 0; j < N; j++)
		{
			//S(i_index, j) = SF[i_index]*S(i_index,j);
			S(i_index, j) *= SF[i_index];
			for (size_t t = 0;  t<i_index; t++){
				S(i_index, j) += SF[t] * S(t, j);
			}

		}

		/* ?? for(size_t r=0; r<N; r++)
		 for(size_t c=0; c < N; c++)
		 {
		 double s = 0;
		 if(c != i_index)
		 s += S(r, c); // SF[c,c] = 1.0
		 //s += S(r, j) * SF(j, c);
		 s += S(r, i_index) * SF[c];
		 S(r, c) = s;
		 }*/

		IF_DEBUG(DBG_FAMG_INTERPOLATION, 5)
		{ S.maple_print("S_SF"); }
	}

	template<typename TMAT>
	void NPrint(TMAT &m, const char *cd)
	{
		UG_LOG("---------------\n")
		for (size_t r = 0; r < m.num_rows(); r++)
			for (size_t c = 0; c < m.num_cols() + 1; c++)
				UG_LOG(cd << "(" << N2[r] << ", " << N2[c] << ") = " << H(r, c) << "\n");
	}

	void calculate_H_from_S_and_Dinv(size_t i)
	{
		const size_t N = S.num_rows(); // = N2
		if (H.num_rows() != onlyN1.size() + 1)
			H.resize(onlyN1.size() + 1, onlyN1.size() + 1);

		// get H = S' Y S'^T
		// todo: use symmetric H.

		for (size_t r = 0; r < onlyN1.size() + 1; r++)
		{
			for (size_t c = r; c < onlyN1.size() + 1; c++)
			{
				smallmat_type s = 0.0;
				for (size_t j = 0; j < N; j++)
					s += S(r, j) * Dinv[j] * MatrixTranspose(S(c, j));
				H(r, c) = s;
				H(c, r) = MatrixTranspose(s);
			}
		}

		//IFI(H.maple_print("H"));
//		IFI(NPrint(H, "H"));

		IF_DEBUG(DBG_FAMG_INTERPOLATION, 6)	{ H.maple_print("H"); }
	}

	template<typename T>
	void BlockToDouble(DenseMatrix<VariableArray2<double> > &Ad, T &A)
	{
		size_t blockSize = GetRows(A(0,0));
		Ad.resize(blockSize*A.num_rows(), blockSize*A.num_cols());

		for(size_t i=0; i<A.num_rows(); i++)
		{
			for(size_t j=0; j<A.num_cols(); j++)
				for(size_t r=0; r<blockSize; r++)
					for(size_t c=0; c<blockSize; c++)
						Ad(i*blockSize+r, j*blockSize+c) = BlockRef(A(i, j), r, c);
		}
	}

	template<typename T>
	void DoubleToBlock(T &A, DenseMatrix<VariableArray2<double> > &Ad, size_t blockSize)
	{
		A.resize(Ad.num_rows()/blockSize, Ad.num_cols()/blockSize);

		for(size_t i=0; i<Ad.num_rows(); i++)
			for(size_t j=0; j<Ad.num_cols(); j++)
				BlockRef(A(i/blockSize, j/blockSize), i%blockSize, j%blockSize) = Ad(i, j);
	}

	/**
	 * calculate_H_from_local_A
	 * calculates locally
	 * - F-Smoothing
	 * 		SF = 1 -D^{-1} A 						on N2 = {i} x {onlyN1, i, onlyN2}
	 * 	so  SF[i] = (1 -D^{-1} A )_( i, N2[i])
	 * - S = SF (1-w D^{-1} A)			 (on N2 x N2)
	 * - H = S Y S^T				 	(on {only1} x {only1})
	 * - Hi[.] = H(., i)				 (on {i} x {only1})
	 *
	 * asserts:
	 * S is precalculated to be A locally on N2 x N2, that is
	 * 		S[i, j] = A( N2[i], N2[j] )
	 *
	 * 	note that S is first A and then overwritten to be S.
	 */
	void calculate_H_from_local_A(size_t i)
	{
		// before this function, S = local A (on only1, i, only2)
		AMG_PROFILE_FUNC();

		UG_ASSERT(S.num_cols() == S.num_rows(), "");
		UG_ASSERT(S.num_cols() == onlyN1.size()+1+onlyN2.size(), "");

		if(m_bHAsInverse == false)
		{
			calculate_Dinv_from_local_A(i);
			calculate_SF_from_Dinv_and_local_A(i);
			calculate_Sjac_from_Dinv_and_local_A();
//			IFI(NPrint(S, "Sjac"));
			calculate_S_from_SF_and_Sjac();
//			IFI(NPrint(S, "S_SF"));
			calculate_H_from_S_and_Dinv(i);
		}
		else
		{
			DenseMatrix<VariableArray2<double> > a, b;
//			S.maple_print("S");
			BlockToDouble(a, S);
			Invert(a);
//			a.maple_print("a");
			DoubleToBlock(H, a, GetRows(S(0,0)));
//			H.maple_print("H");
		}
	}

	// global_to_local_testvectors
	/** calculates the local testvectors
	 * \param node the fine node to be interpolated
	 *	note: das mit den testvektoren ist noch nicht so sicher:
	 *	sie mï¿½ssen ja theoretisch 2 mal geglï¿½ttet werden. 1x normal jacobi,
	 *	und dann ein 2. mal nur auf den feinen knoten. und danach muss noch
	 *	$\frac 1 {\abs{t}_A}$ geteilt werden. Leider weiï¿½ man aber ja vorher noch nicht,
	 *	welche Knoten fein sind. Dh. im Moment ist das so: Man macht einmal Jacobi,
	 *	berechnet dann mal $\frac 1 {\abs{t}_A}$, und macht dann lokal nur noch so eine
	 *	Nachglï¿½ttung, wenn die globalen Testvektoren in lokale Testvektoren umgerechnet werden.
	 */
	void global_to_local_testvectors(size_t node)
	{
		//AMG_PROFILE_FUNC();
		localTestvectors.resize(m_testvectors.size());
		//size_t level = m_amghelper.positions.size()-1;
		//UG_LOG("Getting testvectors for node " << node << ":\n");
		for (size_t k = 0; k < m_testvectors.size(); k++)
		{
			//UG_LOG(" testvector " << k << ":\n");
			localTestvectors[k].resize(onlyN1.size() + 1);

			localTestvectors[k][onlyN1.size()] = m_testvectors[k][node];

			for (size_t j = 0; j < onlyN1.size(); j++)
			{
				localTestvectors[k][j] = m_testvectors[k][onlyN1[j]];
				//UG_LOG(" onlyN1[" << j << "] = " << onlyN1[j] << m_testvectors[k][onlyN1[j]] <<
					//	" position is " << m_amghelper.positions[level][onlyN1[j]] << ":\n");
			}

			// "normalize"
			/*if(k!=0)
				for (size_t j = 0; j < onlyN1.size(); j++)
				{
					localTestvectors[k][j] -= localTestvectors[k][onlyN1.size()];
					//UG_LOG(" onlyN1[" << j << "] = " << onlyN1[j] << m_testvectors[k][onlyN1[j]] <<
						//	" position is " << m_amghelper.positions[level][onlyN1[j]] << ":\n");
				}*/
		}



		// smooth localTestvectors:
		/*for (size_t k = 0; k < m_testvectors.size(); k++)
		{
			smallvec_type s = 0.0;
			for (size_t j = 0; j < onlyN1.size() + 1; j++)
				s += SF[j] * localTestvectors[k][j];

			localTestvectors[k][onlyN1.size()] = s;
			PRINT_VECTOR_AMGDEBUG(5, localTestvectors[k], "\nlocalTestvector " << k);
		}*/
	}

	void global_to_local_testvectors2(size_t node)
	{
		AMG_PROFILE_FUNC();
		localTestvectors.resize(2);
		localTestvectors[0].resize(onlyN1.size() + 1);
		localTestvectors[1].resize(onlyN1.size() + 1);

		size_t level = m_amghelper.positions.size()-1;
		double x=m_amghelper.positions[level][node][0];
		double y=m_amghelper.positions[level][node][1];

		for (size_t j = 0; j < onlyN1.size(); j++)
		{
			localTestvectors[0][j] = 1.0;
			BlockRef(localTestvectors[1][j], 0) = m_amghelper.positions[level][onlyN1[j]][1]-y;
			BlockRef(localTestvectors[1][j], 1) = m_amghelper.positions[level][onlyN1[j]][0]-x;
		}
		localTestvectors[0][onlyN1.size()] = 1.0;
		localTestvectors[1][onlyN1.size()] = 0.0;
		//BlockRef(localTestvectors[1][onlyN1.size()], 0) = m_amghelper.positions[level][node][1];
		//BlockRef(localTestvectors[1][onlyN1.size()], 1) = m_amghelper.positions[level][node][0];

	}

	void calculate_EV_testvectors(size_t node)
	{
		UG_ASSERT(0, "currently disabled");
#if 0
		std::vector<std::vector<size_t> > neighbors(4);
		std::vector<size_t> &onlyN1 = neighbors[1];
		std::vector<size_t> &onlyN2 = neighbors[2];
		std::vector<size_t> &onlyN3 = neighbors[3];
		GetNeighborhoodHierachy(A_OL2, node, 3, neighbors);

		PRINT_VECTOR_AMGDEBUG(3, onlyN1, "\nn1-neighbors");
		PRINT_VECTOR_AMGDEBUG(3, onlyN2, "\nn2-neighbors");
		PRINT_VECTOR_AMGDEBUG(3, onlyN3, "\nn3-neighbors");

		stdvector<size_t> myNeigh;
		myNeigh = onlyN1;
		myNeigh.push_back(node);
		myNeigh.insert(myNeigh.end(), onlyN2.begin(), onlyN2.end());
		myNeigh.insert(myNeigh.end(), onlyN3.begin(), onlyN3.end());

		DenseMatrix<VariableArray2<double> > localA;
		localA.resize(myNeigh.size(), myNeigh.size());

		GetLocalMatrix(A_OL2, localA, &myNeigh[0], &myNeigh[0]);

		localA.maple_print("A");

		size_t N1 = onlyN1.size();
		size_t N2 = onlyN2.size();
		DenseMatrix<VariableArray2<double> > localS;
		localS.resize(N1+1+N2, N1+1+N2);

		for(size_t r=0; r<N1+1+N2; r++)
		{
			for(size_t c=0; c < N1+1+N2; c++)
			localS(r, c) = localA(r, c);
			for(size_t c=N1+1+N2; c < localA.num_cols(); c++)
			localS(r, r) += localA(r, c);
		}

		localS.maple_print("S1");

		// get S = 1-D^{-1} A
		for(size_t r=0; r < N1+1+N2; r++)
		{
			double diagInv = 1/localS(r,r);
			for(size_t c=0; c < N1+1+N2; c++)
			{
				localS(r, c) = -diagInv*localS(r,c);
				if(r==c) localS(r, c) += 1.0;
				else if(localS(r, c) < 1e-12)
				localS(r, c) = 0.0;
			}
		}

		localS.maple_print("S2");

		size_t N = N1+1+N2;
		DenseMatrix<VariableArray2<double> > X;
		X.resize(N, N);
		DenseMatrix<VariableArray2<double> > B;
		B.resize(N, N);
		B = 1.0;
		DenseVector<VariableArray1<double> > lambda;
		lambda.resize(N);

		UG_ASSERT(0, "generalizedEigenvalueProblem has to be fixed.");
		/*int res =*/
		//GeneralizedEigenvalueProblem(localS, X, lambda, B, true);

		// X.maple_print("X");
		// lambda.maple_print("lambda");
		localTestvectors.resize(1);
		localTestvectors[0].resize(onlyN1.size()+1);
		for(size_t i=0; i<onlyN1.size()+1; i++)
		localTestvectors[0][i] = X(i, N-1);
#endif
	}

	// calculate_testvectors
	//---------------------------------------
	void calculate_testvectors(size_t node)
	{
		if (testvectorsExtern == true)
		{
			global_to_local_testvectors(node);
			//add_additional_testvectors_to_H();
		}
		else
			calculate_EV_testvectors(node);
	}

#if 0
	template<typename T1, typename T2>
	void add_vvT(T1 &mat, double alpha, const T2 &v1, const T2 &v2)
	{
		for (size_t r = 0; r < v1.size(); r++)
		{
			for (size_t c = 0; c < v2.size(); c++)
			{
				mat(r, c) += alpha * v1[r] * v2[c];
			}
		}
	}

	void add_vvT(double &mat, double alpha, double &v1, double &v2)
	{
		mat = alpha * v1 * v2;
	}

	/** add_additional_testvectors_to_H
	 * adds additional factors from the testvectors to H, namely
	 * H += \omega_k t^{(k)} t^{(k)}^T
	 * for the testvectors 1 .. m_testvectors.size()-1
	 * (testvector 0 is approximated exactly)
	 */
	void add_additional_testvectors_to_H()
	{
		//AMG_PROFILE_FUNC();
		// skip first vector (it is approximated exactly)
		// /test if other way round is faster/

		for (size_t r = 0; r < onlyN1.size(); r++)
			for (size_t c = 0; c < onlyN1.size(); c++)
			{
				for (size_t k = 1; k < m_testvectors.size(); k++)
				{
					add_vvT(H(r, c), m_omega[k], localTestvectors[k][r],
							localTestvectors[k][c]);
				}
			}
	}
#endif
private:
	// for speedup purposes, we don't want these arrays to be re-allocated all the time,
	// thats why they are stored here.
	stdvector<bool> bvisited; // used for N2-neighborhood calculation

	// todo: instead of VariableArray2, use ReserveableArray2

	// note: matrices/vectors which are "on" N2 are meant like this: Mlocal(i, j) = M( N2[i], N2[j] )

	/// sqrt( diag(A_ii) )
	smallvec_type aiiSqrt;
	smallmat_type aii;

	/// local matrix S = 1 - wD^{-1} A on N2 = {onlyN1, i, onlyN2}
	DenseMatrix<VariableArray2<smallmat_type> > S;

	/// SF[i] = 1 -wD^{-1} A( N2[i], .)  on N2 = {onlyN1, i, onlyN2}
	stdvector<smallmat_type> SF;

	/// diagonal on {onlyN1, i, onlyN2}
	stdvector<smallmat_type> D;

	/// inverse of diagonal on N2 = {onlyN1, i, onlyN2}
	stdvector<smallmat_type> Dinv;

	/// matrix H = S Y S^T  on {onlyN1}
	DenseMatrix<VariableArray2<smallmat_type> > H;

	/// matrix H = S Y S^T  on {i}
	double Hii;

	// matrix H = S Y S^T  on {i} x {onlyN1}   ( hat_h[j] = H(i, onlyN1[j] ) )
	//stdvector<smallmat_type> hat_h;


	/// testvectors on {onlyN1, i}
	stdvector<DenseVector<stdvector<smallvec_type> > > localTestvectors;

	/// testvector on {i}
	double t_i;

	// for the KKT system

	DenseMatrix<VariableArray2<double> > KKT;
	DenseVector<VariableArray1<double> > rhs;
	DenseVector<VariableArray1<double> > lq;
	stdvector<smallmat_type> q;
//	DenseVector<stdvector<double> > t;

	double F; // calculated F

	/// indices which are in distance 1 of i (exactly), so it doesnt contain i
	stdvector<size_t> onlyN1;

	/// indices which are in distance 2 of i (exactly), so it doesnt contain i or onlyN1.
	stdvector<size_t> onlyN2;

	/// N2 = {onlyN1, i, onlyN2} 		(in this order!)
	stdvector<size_t> N2;

	/// temporary storage of possible parent pairs (elements with theta F < Fmin are put into possible parent pairs)
	stdvector<neighborstruct_var<smallmat_ptype> > i_neighborpairs;

	stdvector<stdvector<neighborstruct_var<double> > > i_bNeighborpairs;

	const TMatrix &A;
	const TMatrix &A_OL2;
	std::vector<double> Fb;
	double m_delta;
	double m_theta;
	double m_damping;
	double m_dProlongationTruncation;
	double m_dHReduceInterpolationNodesParameter;
	int m_iHReduceMinParents;
	stdvector<TVector> &m_testvectors;
	stdvector<double> &m_omega;
	stdvector<double> &fvalues;
	const cAMG_helper &m_amghelper;

	bool testvectorsExtern;
	bool m_bHAsInverse;
	double m_threeParentsPenalty;
	size_t m_iMaxParents;
	bool m_bUseConstantIfFailed;

	size_t m_nrTestvectors;

	const std::vector<levelnode> &debugWatchNodes;
	double m_dDeficitSolveTolerance;

	bool is_watch_node(size_t node)
	{
		size_t level = get_level();
		for(size_t i=0; i<debugWatchNodes.size(); i++)
		{
			if(debugWatchNodes[i].level == level &&
					debugWatchNodes[i].node == node)
				return true;
		}
		return false;
	}

};

// end group plugins_amg_famg
/// \}

} // namespace ug
#endif // __H__LIB_ALGEBRA__FAMG_SOLVER__FAMG_INTERPOLATION_CALCULATOR_H__
