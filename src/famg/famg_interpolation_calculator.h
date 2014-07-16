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
#include "misc.h"

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
#define IFI(s) if(is_watch_node(i)) { s; }

#define UG_LOGI(s) IFI(UG_LOG(s));

#define UG_DLOGI(level, s) \
	{if(ug::GetDebugIDManager().get_debug_level(DBG_FAMG_INTERPOLATION) >= level || is_watch_node(i) ) UG_LOG(s); }



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
						UG_LOG(	"P([" << i << "], [" << it.index() << "])[" << k << "," << k <<"] = " << d << "\n");
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

	bool all_Fb_smaller_delta()
	{
		for(size_t i=0; i<Fb.size(); i++)
			if(Fb[i] > m_delta)
				return false;
		return true;
	}
	bool one_Fb_smaller_delta()
	{
		for(size_t i=0; i<Fb.size(); i++)
			if(Fb[i] < m_delta)
				return true;
		return false;

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
					if(one_Fb_smaller_delta())
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
//		int i_min = -1;
//		double f_min = 1e12;

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

					if(!all_Fb_smaller_delta() || F > m_delta)
						continue;

					F *= m_threeParentsPenalty;
					//check_one_interpolation(i, p0, p1);
					UG_DLOGI(2, "F: " << F << "\n");

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
		UG_DLOGI(2, "node [" << i << "] [" << rating.get_original_index(i) << "]\n");

		UG_COND_THROW(!rating[i].is_valid_rating(),
				"node [" << i << "] has no valid rating " << rating[i].get_rating());

		if (testvectorsExtern && BlockNorm2(m_testvectors[0][i]) < 1e-12)
		{
			UG_DLOGI( 2, "testvector in [" << i << "] nearly 0, set as fine\n");
			rating.set_fine(i);
			return false;
		}

		// see local_eigenvalue_test_computation.h
		//calc_special(i);

		GetNeighborhood(A_OL2, i, onlyN1, onlyN2, bvisited);
		if (onlyN1.size() == 0)
		{
			UG_DLOGI(2, "node [" << i << "] has no connections, set as fine\n");
			rating.set_fine(i);
			return false;
		}

		// get H on onlyN1,i_index,onlyN2.
		get_H(i, rating);


		IF_DEBUG(DBG_FAMG_INTERPOLATION, 5)
		{
			for (size_t n = 0; n < onlyN1.size(); n++)
				UG_DLOGI(4, "[" << onlyN1[n] << "] = " << rating.info(onlyN1[n]) << " ");
			UG_DLOGI(4, "\n");
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
				UG_LOG("Best F is " << get_all_neighbors_interpolation(i) << "\n");
		}


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
			UG_DLOGI(2, std::endl << "[" << i  << "]: UNINTERPOLATEABLE" << std::endl);
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
	 *
	 * NOT USED CURRENTLY
	 */
	template<typename TIndices>
	  void reduceKKTSystem(size_t i, TIndices indices, size_t b)
	{
		const size_t N = indices.size();
		std::vector<bool> strong(N, false);
		markAllStrong(i, indices, strong, m_dHReduceInterpolationNodesParameter, m_iHReduceMinParents, b);

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
	bool check_smoother_sufficient(TIndices &indices, size_t b)
	{
		size_t i_index = onlyN1.size();
		const size_t blockSize=GetRows(H(0,0)); // assert blocksize is FIXED!!!
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
	 * set up matrix KKT and rhs
	 * @param i			index to interpolate
	 * @param indices	indices where to interpolate from
	 * @param b			subblock
	 * @return true
	 */
	template<typename TIndices>
	bool  setup_KKT(size_t i, const TIndices &indices, size_t b)
	{
		// does NOT work!
		if(false && check_smoother_sufficient(indices, b) == true)
		{
			UG_LOGI("SMOOTHER SUFFICIENT\n");
			return false;

		}
		const double aii_bb  = BlockRef(aii, b, b);

		const size_t N = indices.size();
		size_t i_index = onlyN1.size();
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


		// reduce to essential entries
		// todo: do less work here (KKT system is too large)

		// if(indices.size() > 4) reduceKKTSystem(indices, b);


		return true;

	}

	/**
	 * calculates the F-value out of lq and H:
	 * F = <lq, H lq> * (A_ii)_bb
	 * @param i			index to interpolate
	 * @param indices	indices where to interpolate from
	 * @param b			subblock
	 * @return F
	 * throws if F < 0
	 */
	template<typename TIndices>
	double calculate_F(size_t i, const TIndices &indices, size_t b)
	{
		const size_t N = indices.size();
		size_t i_index = onlyN1.size();
		Hii = BlockRef(H(i_index, i_index), b, b);
		double lF=Hii;

		for(size_t k=0; k<N; k++)
			lF -= 2*lq[k]*BlockRef(H(i_index, indices[k]), b, b);
		for(size_t k1=0; k1<N; k1++)
			for(size_t k2=0; k2<N; k2++)
				lF += lq[k1]*BlockRef(H(indices[k1], indices[k2]), b, b)*lq[k2];
		lF *= BlockRef(aii, b, b);

		IFI(
			setup_KKT(i, indices, b);
			UG_LOG("F = " << lF << "\n");
			UG_LOG("i = " << i << " b = " << b << "\n");
			UG_LOG(JuliaString(KKT, "KKT"))
			rhs.maple_print("rhs");
			lq.maple_print("lq");

			UG_LOGN("Fb[" << b << "] = " << lF << ")");
			UG_LOG("\n");
		)

		if(lF<-1e-10)
		{
			return 1000;
			// ERROR: negative F !
			// recalc KKT to display error.
			setup_KKT(i, indices, b);
			UG_LOG("ERROR: negative F values are impossible (i.e., find bug and fix it!)\n");
			UG_LOG("F = " << lF << "\n");
			UG_LOG("i = " << i << " b = " << b << "\n");
			UG_LOG(JuliaString(KKT, "KKT"))
			rhs.maple_print("rhs");
			lq.maple_print("lq");

			UG_LOG("WARNING: Negative entries at node [" << i << "]: ");
			for(size_t l=0; l<indices.size(); l++)
			{	UG_LOG("p" << l << " = " << onlyN1[ indices[l] ] << " " << lq[l] << "\n"); }
			UG_THROW("negative F");
		}
		return lF;
	}

	/**
	 * if an interpolation value in lq is negative, this function will return its index.
	 * values between - and +1e-5 are set to zero.
	 * @param i			index to interpolate
	 * @param indices	indices where to interpolate from
	 * @param b			subblock
	 * @return index of negative lq, else -1
	 */
	template<typename TIndices>
	int check_for_negative_lq(size_t i, const TIndices &indices, size_t b)
	{
		IF_DEBUG(DBG_FAMG_INTERPOLATION, 5) lq.maple_print("lq");
		//bool bHasNegative=false;
		int iNegative=-1;
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
					UG_LOG("WARNING: Negative entries at node [" << i << "]: ");
					for(size_t l=0; l<indices.size(); l++)
					{	UG_LOG("p" << l << " = " << onlyN1[ indices[l] ] << " "); }
					UG_LOG(" solve_KKT: lq[" << k << "] = " << lq[k] << "\n");
				}

				iNegative = k;
			}
			BlockRef(q[k], b, b) = lq[k];
		}
		return iNegative;
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
//		const size_t M = localTestvectors.size();

		TemporaryDebugLevel TDL(DBG_FAMG_INTERPOLATION); IFI(TDL.set(5));
		const size_t blockSize=GetRows(H(0,0)); // assert blocksize is FIXED!!!
		for(size_t k=0; k<indices.size(); k++) q[k] = 0.0;
		F=0;
		Fb.resize(blockSize);
		bool bHasNegative=false;
		for(size_t b=0; b<blockSize; b++)
		{
			UG_LOGI("BLOCK [ " << b << " ]\n");
			if(setup_KKT(i, indices, b) == false) continue;
			if(SolveDeficit(KKT, lq, rhs, m_dDeficitSolveTolerance)==true)
			{
				if(check_for_negative_lq(i, indices, b) != -1)
				{
					UG_DLOGI(2, "[" << i << "], block " << b << " has negatives.\n");
					if(bCalcNegatives==false) return -1;
					bHasNegative = true;
					continue;
				}
				Fb[b] = calculate_F(i, indices, b);
				F += Fb[b];
			}
			else
			{
				Fb[b] = 1e12;
				F += 1e12;
				IF_DEBUG(DBG_FAMG_INTERPOLATION, 2)
				{
				  UG_LOG(" solve_KKT: could not solve KKT (decomposition failed) at node [" << i << "]\n");
				}
			}
		}
		IFI(UG_LOG("SUCCESS, F = " << F << "\n"));
		// TODO: Sinvoll?
		F /= blockSize;

		IFI(UG_LOG("SUCCESS, bCalcNegatives = " << bCalcNegatives << ", bHasNegative = " << bHasNegative << "\n"));
		return (bCalcNegatives==false && bHasNegative == true) ? -1 : 0;
		//return 0;
	}

	/**
	 * same as solve_KKT, only that if there are negative interpolation values, this function returns the index of them via &neg.
	 */
	template<typename TIndices>
	int solve_KKT_get_neg(size_t i, const TIndices &indices, size_t &neg)
	{
		TemporaryDebugLevel TDL(DBG_FAMG_INTERPOLATION); IFI(TDL.set(5));
		const size_t blockSize=GetRows(H(0,0)); // assert blocksize is FIXED!!!
		for(size_t k=0; k<indices.size(); k++) q[k] = 0.0;
		F=0;
		Fb.resize(blockSize);
//		bool bHasNegative=false;
		for(size_t b=0; b<blockSize; b++)
			Fb[b] = -1;
		for(size_t b=0; b<blockSize; b++)
		{
			UG_LOGI("BLOCK [ " << b << " ]\n");
			if(setup_KKT(i, indices, b) == false) continue;
			if(SolveDeficit(KKT, lq, rhs, m_dDeficitSolveTolerance)==true)
			{
				int iNeg = check_for_negative_lq(i, indices, b);
				if(iNeg != -1)
				{
					UG_DLOGI(3, "node [" << i << "] has negative interpolation from " << indices[iNeg] << "\n");
					neg = iNeg;
					return -1;
				}
				Fb[b] = calculate_F(i, indices, b);
				F += Fb[b];
			}
			else
			{
				Fb[b] = 1e12;
				F += 1e12;
				UG_DLOGI(2, "node [" << i << "] could not solve KKT (decomposition failed)\n");
			}
		}
		IFI(UG_LOG("SUCCESS get_neg, F = " << F << "\n"));
		// TODO: Sinvoll?
		F /= blockSize;

		return 0;
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

	/**
	 * check if a interpolation from either p0 or p1 alone is sufficient
	 */
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
//		return BlockRef(H(i_index, j), b,b) / BlockRef(H(j, j), b,b);
		int blockSize = GetRows(H(i_index, j));
		double m=0;
		for(b=0; b<blockSize; b++)
			m=std::max(m, BlockRef(H(i_index, j), b,b) / BlockRef(H(j, j), b,b));
		return m;
	}


	// select only coarse neighbors w/ sufficiently large couplings
	template<typename TIndices>
	void markAllStrong(size_t i, const TIndices &coarseNeighbors,
			   std::vector<bool> &strongMarks, double eps, size_t nnodes, int b=0)
	{

		UG_DLOGI(3, "markAllStrong with sigma = "<<eps <<" for comp=" << b << std::endl);
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
				UG_DLOGI(3, val<< " ")
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
				UG_DLOGI(3, val<< " ")
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
	void markAllStrongAbs(size_t i, const TIndices &coarseNeighbors,
				std::vector<bool> &strongMarks, double eps, int b=0)
	{

		UG_DLOGI(3, "markAllStrongAbs with sigma = "<<eps <<" for comp=" << b << std::endl);
		strongMarks.resize(coarseNeighbors.size());

		// determine maximum of offdiag/diag in i-th row
		double maxVal = 0.0;
		for (size_t j = 0; j < coarseNeighbors.size(); j++)
		{
			const double val = fabs( H_strong_measure(coarseNeighbors[j], b));
			maxVal = std::max(maxVal, val);
			UG_DLOGI(3, val<< " ")
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

	std::vector<double> strongMeasure;
	/**
	 *
	 * @param ind				index in onlyN1 -> onlyN1[ind[i]] is global index
	 * @param coarseNeighbors	subset of indices of onlyN1
	 * @param b
	 */

	void sortByStrength(std::vector<size_t> &ind, int b=-1)
	{
		strongMeasure.resize(onlyN1.size());

		// determine maximum of offdiag/diag in i-th row
//		double maxVal = 0.0;
		size_t blockSize=GetRows(H(0,0));
		// invert with - because bigger is better -> sort is small to big.
		for (size_t j = 0; j < onlyN1.size(); j++)
		{
//			UG_LOG("onlyN1 " << j << " = " << onlyN1[j] << "\n");
			strongMeasure[j] = 0;
			if(b == -1)
				for(size_t bi=0; bi<blockSize; bi++)
					strongMeasure[j] -= fabs(H_strong_measure(j, bi));
			else
				strongMeasure[j] = -fabs(H_strong_measure(j, b));
		}

//		UG_LOG("Strong M easures:\n");
//		PRINT_VECTOR(ind, "ind");
//		for(size_t i=0; i<ind.size(); i++)
//			UG_LOG(onlyN1[ind[i]] << " = " << strongMeasure[ind[i]] << "\n");
		std::sort(ind.begin(), ind.end(), CompareIndicesBy(strongMeasure));
	}

	/**
	 * @param ind				index in onlyN1 -> onlyN1[ind[i]] is global index
	 * @param b
	 */
	void sortByStrength1(std::vector<size_t> &ind, int b=0)
	{
		ind.resize(onlyN1.size());
		for(size_t i=0; i<onlyN1.size(); i++)
			ind[i] = i;
		sortByStrength(ind, b);
	}

	void getAllStrong(size_t i, std::vector<size_t> &res, double tol, int nnodes=1)
	{
		std::vector<size_t> indices;
		for(size_t i=0; i<onlyN1.size(); i++) indices.push_back(i);
		std::vector<bool> strongMarksAll(indices.size(), false);
		std::vector<bool> strongMarks(indices.size(), false);
		size_t blockSize=GetRows(H(0,0));
		for(size_t b=0; b<blockSize; b++)
		{
		  markAllStrong(i, indices, strongMarks, tol, nnodes, b);
		  for(size_t i=0; i<strongMarks.size(); i++)
			  strongMarksAll[i] = strongMarksAll[i] && strongMarks[i];
		}
		res.clear();
		for(size_t i=0; i<strongMarks.size(); i++)
			if(strongMarks[i]) res.push_back(i);

	}



	/*
	 * select only coarse neighbors w/ sufficiently large couplings
	 * CURRENTLY NOT USED:
	 * our current approach is to sort the coarse neighbors depending on their strong-connection-measure
	 * and from this list select 2, 3, and so on.
	 */
#if 0
	void reduceToStrong_H_strong_measure(const std::vector<size_t> &coarseNeighbors,
			std::vector<size_t> &strongNeighbors, double eps, int b=0)
	{

		// determine maximum among coarse neighbors
		double maxVal = H_strong_measure(coarseNeighbors[0], b);
		for (size_t j = 1; j < coarseNeighbors.size(); j++)
			maxVal = std::max(H_strong_measure(coarseNeighbors[j], b), maxVal);

		// select strongly coupled coarse neighbors
		for (size_t j = 0; j < coarseNeighbors.size(); j++)
			if (H_strong_measure(coarseNeighbors[j], b) >= eps * maxVal)
				strongNeighbors.push_back(coarseNeighbors[j]);
	}

	void reduceToStrong_BlockNorm(const std::vector<size_t> &coarseNeighbors,
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
#endif


	/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	/// solve the KKTs from external coarsening
	// key here is
	// 1. we want to interpolate sufficiently well (F < F_max = delta)
	// 2. we don't want to use too many interpolating nodes
	// 3. we don't want negative interpolation values,
	// (unless we are on the last multigrid level: fix this).

	size_t csolve_KKT_non_neg;
	/**
	 * start with all indices
	 * loop:
	 *   try to solve with those indices
	 *   if successfully solved, return
	 *   else remove one index
	 * ATTENTION! This function can have VERY LARGE runtime since it basically checks ALL configurations.
	 */
	int solve_KKT_non_neg(size_t i, std::vector<size_t> &indices, int tvSize, bool bCalcNegatives=false)
	{
		csolve_KKT_non_neg++;
		if(csolve_KKT_non_neg > 20) return -1;
		set_sizes(indices.size(), tvSize);
		IFI(UG_LOG("tvSize = " << tvSize << " (");		for(size_t j=0; j<indices.size(); j++)		{ UG_LOG(onlyN1[indices[j]] << " "); }		UG_LOG(")\n");		)
		int ret = solve_KKT(i, indices);
//		UG_LOGI(" ret = " << ret << ": " << Fb[0] << " " << Fb[1] << "\n");
		if(ret != -1)
			return ret;
		if(indices.size() <= 3) return -1;
		for(size_t j=0; j<indices.size(); j++)
		{
			std::vector<size_t> mI = indices;
			mI.erase(mI.begin()+j);
			ret = solve_KKT_non_neg(i, mI, tvSize);
			if(ret == 0 && all_Fb_smaller_delta())
			{
				indices = mI;
				return ret;
			}
		}
		return -1;
	}


	/**
	 * start with the best 6 indices
	 * loop:
	 *   try to solve with those indices
	 *   if successfully solved, return
	 *   else remove the index with negative values or the last one (with worse H-rating)
	 */
	int solve_KKT_non_neg3(size_t i, std::vector<size_t> &indices, int tvSize, bool bCalcNegatives=false)
	{
		size_t neg;
		size_t K=indices.size();
		if(K > 8) K = 8;

		for(; K>=2; K--)
		{

			IFI(PRINT_VECTOR(indices, "indices"));
			set_sizes(K, tvSize);
			indices.resize(K);
			int ret = solve_KKT_get_neg(i, indices, neg);
			UG_DLOGI(3, "node [" << i << "], K = " << K << " ret = " << ret << " F = " << F << "\n");
			if(ret == 0)
			{
				 if(all_Fb_smaller_delta())
					 return 0;
				 else
					 return -1;
			}
			std::swap(indices[K-1], indices[neg]);

		}
		UG_DLOGI(3, "get_neg -1\n");
		return -1;
	}

	/**
	 * start with the best 2 indices
	 * loop:
	 *   try to solve with those indices
	 *   if successfully solved, return
	 *   else add the next best index
	 *   if an index was negative, remove it from list.
	 */
	int solve_KKT_non_neg2(size_t i, std::vector<size_t> &indices, int tvSize, bool bCalcNegatives=false)
	{
		size_t neg;
		size_t nrNeg=indices.size()+1;

		for(size_t K=0; K<nrNeg;)
		{

			SubVector vec(indices, K);
			if(is_watch_node(i))
				{PRINT_VECTOR(vec, "vec, " << i << " K = " << K); }
			set_sizes(K, tvSize);
			int ret = solve_KKT_get_neg(i, vec, neg);
			if(ret == 0)
			{
				 if(all_Fb_smaller_delta())
				 {
					 if(is_watch_node(i)) UG_LOG("success\n");
					 indices.resize(K);
					 return 0;
				 }
				 if(is_watch_node(i)) UG_LOG("delta\n");
				 K++;
			}
			else
			{
				nrNeg--;
				std::swap(indices[nrNeg-1], indices[neg]);
				if(is_watch_node(i)) UG_LOG("neg\n");
			}
		}
		return -1;
	}


	bool direct_interpolation(size_t i, TPMatrix &P,
			FAMGNodes<TPMatrix> &rating, std::vector<size_t> &coarseNeighbors, bool bIndirect)
	{
//		UG_LOGI("direct interpolation [" << i << "]\n");
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
//		UG_LOG("DIRECT INTERPOLATION " << i << "\n");
		std::vector<size_t> myNeighbors;
		IFI(PRINT_VECTOR(coarseNeighbors, "coarseNeighbors"));
		sortByStrength(coarseNeighbors);

		if(is_watch_node(i))
		{
			UG_LOG("sorted strong: ");
			for(size_t j=0; j<coarseNeighbors.size(); j++) { UG_LOG(onlyN1[coarseNeighbors[j]] << ", "); }
			UG_LOG("\n");
		}
		myNeighbors = coarseNeighbors;
//		PRINT_VECTOR(myNeighbors, "myNeighbors");
		//std::vector<size_t> myNeighbors=coarseNeighbors;
//		reduceToStrong(coarseNeighbors, myNeighbors, 0.7);
		csolve_KKT_non_neg=0;

		bool bSuccess = (solve_KKT_non_neg3(i, myNeighbors, -1) == 0);

		if(bSuccess)
		{
			if(F > m_delta)
			{
				bSuccess = false;
				UG_DLOGI(0, "node [" << i << "] direct interpolation not good enough\n");
			}
		}
		else
		{
			UG_DLOGI(0, "node [" << i << "] direct interpolation not possible\n");
		}


		if(!bSuccess && m_bUseConstantIfFailed)
		{
			bSuccess = solve_KKT_non_neg(i, myNeighbors, 1) == 0;
			if(bSuccess)
			{
				UG_DLOGI(0, "node [" << i << "] direction interpolation with constant successfull\n");
			}
		}

		if(!bSuccess || F > m_delta)
		{
			UG_DLOGI(3, "direct_interpolation, had to set node [" << i << "] coarse: F = " << F);

			rating.external_uninterpolateable(i);
		}
		else
		{
			UG_DLOGI(3, "coarse neighbors, Interpolating from ");
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
					UG_DLOGI(3, "[" << node << "]: " << q[j] << ", ");
					P(i, node) = q[j];
				}
			 }
			else
			{
				// todo: modify dropping algorithm
				// nodes should be dropped if they have (in relation) less influence than others
				if(bIndirect)
				{
					UG_DLOGI(3, "node [" << i << "] indirect interolation\n")
					smallmat_type sm;
					smallmat_type droppedSum=1;
					for(size_t j=0; j<N; j++)
					{
						if(BlockNorm(q[j]) < 1e-2) continue;
						size_t nodej = onlyN1[myNeighbors[j]];
						if(rating[nodej].is_dirichlet()) continue; // because defect is 0 in those nodes.
						if(!rating[nodej].is_coarse())
						{
							for(typename prolongation_matrix_type::row_iterator it = P.begin_row(nodej); it != P.end_row(nodej); ++it)
							{
								sm = q[j] * it.value();
								if(BlockNorm(sm) < 0.02) { droppedSum-=sm; continue;}
								P(i, it.index()) += q[j] * it.value();
							}
						}
						else
							P(i, nodej) += q[j];
						UG_DLOGI(3, "[" << nodej << "]: " << q[j] << ", ");
					}
					Invert(droppedSum);
					for(typename prolongation_matrix_type::row_iterator it = P.begin_row(i); it != P.end_row(i); ++it)
						it.value() *= droppedSum;


				}
				else
				{
					UG_DLOGI(3, "node [" << i << "] directinterolation\n")
					for(size_t j=0; j<N; j++)
					{
						if(BlockNorm(q[j]) < 1e-2) continue;
						size_t nodej = onlyN1[myNeighbors[j]];
						if(rating[nodej].is_dirichlet()) continue; // because defect is 0 in those nodes.
						if(!rating[nodej].is_coarse()) rating.set_coarse(nodej);
						P(i, nodej) = q[j];
						UG_DLOGI(3, "[" << nodej << "]: " << q[j] << ", ");
					}
				}
			}
			check_weights(P, i);
			if(fvalues.size() > i) fvalues[i] = F;
//				UG_LOG("setting [" << i << "] fine.");
			if(bIndirect)
			{
				UG_DLOGI(3, "node [" << i << "] set aggressive fine\n");
				rating.set_aggressive_fine(i);
			}
			else
			{
				UG_DLOGI(3, "node [" << i << "] set fine\n");
				rating.set_fine(i);
			}
		}

		return true;
	}

	/**
	 * this doesn't seem to work too well with linear elasticity, since
	 * the AMG-connection rating is somehow wrong.
	 */
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
			if (rating[conn.index()].is_interpolateable_from())
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

	/**
	 * just put all coarse neighbors into the list
	 * if indirect interpolation is used, also put in "inner" direct fine nodes.
	 */
	void get_coarse_neighborhood(const TMatrix &A,
			const FAMGNodes<TPMatrix> &rating, size_t i, stdvector<size_t> &onlyN1,
			stdvector<size_t> &onlyN2, stdvector<size_t> &coarseNeighbors,
			stdvector<bool> &bvisited, bool bIndirect)
	{

//		bvisited.resize(0);
//		bvisited.resize(A.num_rows(), false);
		GetNeighborhood(A_OL2, i, onlyN1, onlyN2, bvisited);

		for (size_t j = 0; j < onlyN1.size(); j++)
		{
			size_t k=onlyN1[j];
			//if (!rating[onlyN1[j]].is_fine()) // && !rating[onlyN1[j]].is_uninterpolateable())
			if (rating[k].is_interpolateable_from())
				coarseNeighbors.push_back(j);
			if (bIndirect && rating[k].is_fine_direct() && rating.i_must_assign(k))
				coarseNeighbors.push_back(j);
		}
	}

#if 0
	// not used
	double print_strenght_sorted_neighbors(size_t i)
	{
		std::vector<size_t> strong;
		sortByStrength1(strong);
		UG_LOG("sorted strong: ");
		for(size_t j=0; j<strong.size(); j++)
		{UG_LOG(onlyN1[strong[j]] << ", ");}
			UG_LOG("\n");


	}
#endif


	double get_all_neighbors_interpolation(size_t i)
	{
		std::vector<size_t> strong;
		sortByStrength1(strong);

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
			FAMGNodes<TPMatrix> &rating, bool bIndirect)
	{
		//UG_LOG("get_all_neighbors_interpolation [" << i << "]\n");
		AMG_PROFILE_FUNC();
//		UG_DLOG(DBG_FAMG_INTERPOLATION, 3,
		UG_LOGI("get_all_neighbors_interpolation on node " << rating.get_original_index(i) << "\n")

		stdvector<size_t> coarseNeighbors;

		get_coarse_neighborhood(A_OL2, rating, i, onlyN1, onlyN2,
				coarseNeighbors, bvisited, bIndirect);


		IFI(PrintVector(onlyN1, "N1"))
		IFI(PrintVector(coarseNeighbors, "cN1"))
//		UG_LOG("coarseNeighbors: " << coarseNeighbors.size() << "\n");
		if (onlyN1.size() == 0)
		{
			rating.set_dirichlet(i);
			return true;
		}
		else if(coarseNeighbors.size() == 0)
		{
			rating.set_coarse(i);
		}
		return direct_interpolation(i, P, rating, coarseNeighbors, bIndirect);
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
		PROFILE_BEGIN(N2);
		N2 = onlyN1;
		N2.push_back(i);
		N2.insert(N2.end(), onlyN2.begin(), onlyN2.end());
		PROFILE_END();

		PRINT_VECTOR_AMGDEBUG(2, onlyN1, "\nonlyN1");
		PRINT_VECTOR_AMGDEBUG(2, onlyN2, "\nonlyN2");
		PRINT_VECTOR_AMGDEBUG(2, N2, "\nN2");

		// 2. get submatrix in A_OL2 on N2
		//AMG_PROFILE_NEXT(AMG_H_GetLocalMatrix);

		S.resize(N2.size(), N2.size());
		S = 0.0;

		PROFILE_BEGIN(GetLocalMatrix)
		GetLocalMatrix(A_OL2, S, &N2[0], &N2[0]);
		PROFILE_END_(GetLocalMatrix);

		PROFILE_BEGIN(p2)
		size_t i_index = onlyN1.size();
		GetDiagSqrt(aiiSqrt, S(i_index, i_index));
		GetDiag(aii, S(i_index, i_index));
		PROFILE_END();

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
		AMG_PROFILE_FUNC();
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

	void AssignMult3(double &d, double alpha, double &a, double &b)
	{
		d = alpha*a*b;
	}

	template<typename A, typename B, typename C>
	inline void AssignMult3(DenseMatrix<A> &dest, double alpha, const DenseMatrix<B> &mA, const DenseMatrix<C> &mB)
	{
		UG_ASSERT(mA.num_cols() == mB.num_rows(), "");
		dest.resize(mA.num_rows(), mB.num_cols());
		for(size_t r=0; r < mA.num_rows(); r++)
			for(size_t c=0; c < mB.num_cols(); c++)
			{
				dest(r, c) = mA(r, 0)* mB(0, c);
				for(size_t k=1; k < mB.num_rows(); k++)
					dest(r, c) += mA(r, k) * mB(k, c);
				dest(r, c) *= alpha;
			}
	}

	void calculate_SF_from_Dinv_and_local_A(size_t i)
	{
		AMG_PROFILE_FUNC();
		const size_t N = S.num_rows(); // = N2
		const size_t i_index = onlyN1.size();

		if (SF.size() != N)
			SF.resize(N);

		// get SF = 1-wDF^{-1} A  (F-smoothing)
		//AMG_PROFILE_BEGIN(AMG_HA_calculate_SF);

		// bei f-smoothing nie und nimmer damping (arne 3.juni)
		if (BlockNorm(S(i_index, i_index)) < 1e-12)
			UG_LOG("i=[" << i << "] S(i,i) = " << S(i,i) << "\n");

		// todo: only go to onlyN1.size(), since S(i_index, j) = 0 forall j > onlyN1.size()
		// SF[i] = (1 -wD^{-1} A )_( i, N2[i])
		for (size_t j = 0; j < N; j++)
			AssignMult3(SF[j], -1.0, Dinv[i_index], S(i_index, j));
		SF[i_index] += 1.0;

		PRINT_VECTOR_AMGDEBUG(5, SF, "SF");
	}

	void calculate_Sjac_from_Dinv_and_local_A()
	{
		AMG_PROFILE_FUNC();
		const size_t N = S.num_rows(); // = N2
		// get S = 1-w D^{-1} A
		smallmat_type m;
		for (size_t r = 0; r < N; r++)
		{
			// todo: only go to onlyN1.size(), since S(i_index, j) = 0 forall j > onlyN1.size()
			for (size_t c = 0; c < N; c++)
			{
				AssignMult3(m, -m_damping, Dinv[r], S(r, c));
				S(r, c) = m;
				if (r == c)
					S(r, c) += 1.0;
			}
		}

		IF_DEBUG(DBG_FAMG_INTERPOLATION, 5)
			S.maple_print("Sjac");
	}

	void calculate_S_from_SF_and_Sjac()
	{
		AMG_PROFILE_FUNC();
		const size_t N = S.num_rows(); // = N2
		const size_t i_index = onlyN1.size();
		smallmat_type m;

		// get S' = SF S
		// (possible without temporary since SF is mostly Id,
		//  and then SF S is addition of rows of S)

		// S'(i_index, i) = sum_t SF(i_index,t) * S(t, i)
		for (size_t j = 0; j < N; j++)
		{
			//S(i_index, j) = SF[i_index]*S(i_index,j);
			AssignMult(m, S(i_index, j), SF[i_index]);
			S(i_index, j) = m;
			for (size_t t = 0;  t<i_index; t++){
				AddMult(S(i_index, j), SF[t], S(t, j));
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


	inline void AssignMultT(double &dest, const double &mA, const double &mB)
	{
		dest = mA*mB;
	}

	template<typename A, typename B, typename C>
	inline void AssignMultT(DenseMatrix<A> &dest, const DenseMatrix<B> &mA, const DenseMatrix<C> &mB)
	{
		UG_ASSERT(mA.num_cols() == mB.num_rows(), "");
		dest.resize(mA.num_rows(), mB.num_cols());
		for(size_t r=0; r < mA.num_rows(); r++)
			for(size_t c=0; c < mB.num_cols(); c++)
			{
				AssignMultT(dest(r, c), mA(r, 0), MatrixTranspose(mB(c, 0)));
				for(size_t k=1; k < mB.num_rows(); k++)
					AddMult(dest(r, c), mA(r, k), MatrixTranspose(mB(c, k)));
			}
	}


	inline void AssignT(double &dest, const double &mA)
	{
		dest = mA;
	}

	template<typename A, typename B>
	inline void AssignT(DenseMatrix<A> &dest, const DenseMatrix<B> &mA)
	{
		UG_ASSERT(mA.num_cols() == mA.num_rows(), "");
		dest.resize(mA.num_rows(), mA.num_cols());
		for(size_t r=0; r < mA.num_rows(); r++)
			for(size_t c=0; c < mA.num_cols(); c++)
				AssignT(dest(r, c), mA(c, r));
	}



	void calculate_H_from_S_and_Dinv(size_t i)
	{
//		PROFILE_FUNC();
		const size_t N = S.num_rows(); // = N2
		if (H.num_rows() != onlyN1.size() + 1)
			H.resize(onlyN1.size() + 1, onlyN1.size() + 1);

		// get H = S' Y S'^T
		// todo: use symmetric H.

		for (size_t c = 0; c < onlyN1.size() + 1; c++)
			for (size_t r = c; r < onlyN1.size() + 1; r++)
				H(r, c) = 0.0;

		smallmat_type m;
		for (size_t j = 0; j < N; j++)
			for (size_t c = 0; c < onlyN1.size() + 1; c++)
			{
				AssignMultT(m, Dinv[j], S(c, j));
				for (size_t r = c; r < onlyN1.size() + 1; r++)
				{
					// A_ij = B_ik C_kj
					// M_ij = F_ik A_kj = F_ik B_km C_mj

					AddMult(H(r, c), S(r, j), m);

				}
			}

		for (size_t c = 0; c < onlyN1.size() + 1; c++)
			for (size_t r = c; r < onlyN1.size() + 1; r++)
				AssignT(H(c, r), H(r, c));

		//IFI(H.maple_print("H"));
//		IFI(NPrint(H, "H"));

		IF_DEBUG(DBG_FAMG_INTERPOLATION, 6)	{ H.maple_print("H"); }
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
			//IFI(NPrint(S, "locA"));
			calculate_Dinv_from_local_A(i);
			calculate_SF_from_Dinv_and_local_A(i);
			calculate_Sjac_from_Dinv_and_local_A();
			//IFI(NPrint(S, "Sjac"));
			calculate_S_from_SF_and_Sjac();
			//IFI(NPrint(S, "S_SF"));
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
	 */
	void global_to_local_testvectors(size_t node)
	{
		//AMG_PROFILE_FUNC();
		localTestvectors.resize(m_testvectors.size());
		for (size_t k = 0; k < m_testvectors.size(); k++)
		{
			localTestvectors[k].resize(onlyN1.size() + 1);
			localTestvectors[k][onlyN1.size()] = m_testvectors[k][node];

			for (size_t j = 0; j < onlyN1.size(); j++)
				localTestvectors[k][j] = m_testvectors[k][onlyN1[j]];
		}
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

	void calculate_EV_testvectors(size_t i)
	{
		UG_LOG("\n\n\n");
		PROFILE_FUNC_GROUP("debug");
		std::vector<std::vector<size_t> > indicesH;
		GetNeighborhoodHierachy(A_OL2, i, myTmpFAMGParam_LPE_N , indicesH);

		std::vector<size_t> indices;
		PROFILE_BEGIN(A);

		for(size_t k=0; k < indicesH.size(); k++)
			indices.insert(indices.end(), indicesH[k].begin(), indicesH[k].end());

		size_t indices1=indicesH[0].size();
		size_t indices2=indicesH[1].size()+indices1;
		size_t indices3=indicesH[2].size()+indices2;
		size_t indices4=indices.size();
		//indices3 = indices4;

		DenseMatrix<VariableArray2<smallmat_type> > Aloc;

		Aloc.resize(indices.size(), indices.size());

		GetLocalMatrix(A_OL2, Aloc, &indices[0], &indices[0]);

		DenseMatrix<VariableArray2<smallmat_type> > Bloc;
		Bloc.resize(indices4, indices4);
		for(size_t j=0; j<indices4; j++)
			Bloc(j,j) = Aloc(j,j);

		stdvector<DenseVector<VariableArray1<smallvec_type> > > tv;
		tv.resize(m_testvectors.size());
		DenseVector<VariableArray1<smallvec_type> > tmp;

		for(size_t j=0; j<m_testvectors.size(); j++)
			GetLocalVector(m_testvectors[j], tv[j], indices);

		DenseMatrix<VariableArray2<double> > ewA, ewB;
		ewA.resize(tv.size(), tv.size());
		ewB.resize(tv.size(), tv.size());

		for(size_t j1=0; j1<tv.size(); j1++)
		{
			tmp = Aloc*tv[j1];
			for(size_t k=indices3; k<indices4; k++) tmp[k] *= myTmpFAMGParam1;
			for(size_t j2=0; j2<tv.size(); j2++)
				ewA(j1, j2) = VecProd(tmp, tv[j2]);
		}
		for(size_t j1=0; j1<tv.size(); j1++)
		{
			tmp = Bloc*tv[j1];
	//		for(size_t k=indices3; k<indices4; k++) tmp[k] = 0.0;
			for(size_t j2=0; j2<tv.size(); j2++)
				ewB(j1, j2) = VecProd(tmp, tv[j2]);
		}

		DenseMatrix<VariableArray2<double> > ev;
		std::vector<std::complex<double> > lambda;
		ev.resize(ewA.num_rows(), ewA.num_rows());
		lambda.resize(ewA.num_rows());
		GeneralizedEigenvalueProblemComplex(ewA, ev, lambda, ewB, true);


		size_t numberOfReducedTestvectors=1;
		localTestvectors.resize(numberOfReducedTestvectors);


		PRINT_VECTOR(onlyN1, "onlyN1");
		PRINT_VECTOR(indices, "indices");
		for(size_t k=0; k<numberOfReducedTestvectors; k++)
		{
			tmp = 0.0;
			for(size_t c = 0; c<ev.num_rows(); c++)
				tmp += ev(c, k) * tv[c];

			writeSubVec2(Stringify() << m_matrixWritePath << "/PEV_" << i << "_ev_" << k << "__" << lambda[k] << ".vec", Aloc, tmp, indices, indices4);

			localTestvectors[k].resize(onlyN1.size() + 1);
			localTestvectors[k][onlyN1.size()] = tmp[0];

			for (size_t j = 0; j < onlyN1.size(); j++)
			{
				UG_COND_THROW(onlyN1[j] != indices[j+1], "?");
				localTestvectors[k][j] = m_testvectors[k][j+1];
			}
		}
		UG_LOG("\n\n\n");
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
