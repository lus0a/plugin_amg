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


/**
 * FAMGInterpolationCalculator
 *
 * This class calculates the interpolation weights. It should perhaps be merged with FAMGLevelCalculator
 * its only public member functions are
 * void get_possible_parent_pairs(size_t i, stdvector<neighborstruct2> &possible_neighbors, famg_nodes &rating)
 * void get_all_neighbors_interpolation(size_t i, TPMatrix &P, famg_nodes &rating)
 */
template<typename TMatrix, typename TPMatrix, typename TVector>
class FAMGInterpolationCalculator
{
private:
	typedef typename TMatrix::value_type smallmat_type;
	typedef typename TVector::value_type smallvec_type;
	typedef typename TPMatrix::value_type smallmat_ptype;
	FAMGInterpolationCalculator(
			const FAMGInterpolationCalculator<TMatrix, TPMatrix, TVector> &other);

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

	template<typename T>
	FAMGInterpolationCalculator(const TMatrix &A_,
			const TMatrix &A_OL2_, const FAMG<T> &famg,
			stdvector<TVector> &testvectors, stdvector<double> &omega,
			stdvector<double> &_fvalues) :
			A(A_), A_OL2(A_OL2_), m_testvectors(testvectors), m_omega(omega), fvalues(_fvalues),
			m_amghelper(famg.m_amghelper)
	{
		m_delta = famg.get_delta();
		m_theta = famg.get_theta();
		m_bHAsInverse = famg.m_bHAsInverse;
		m_damping =
				famg.get_damping_for_smoother_in_interpolation_calculation();
		m_dProlongationTruncation = famg.get_prolongation_truncation();
		m_dHReduceInterpolationNodesParameter =
				famg.get_H_reduce_interpolation_nodes_parameter();
		m_iMaxParents = famg.m_iMaxParents;

		testvectorsExtern = (m_testvectors.size() > 0);
	}

	void init()
	{
		bvisited.resize(A_OL2.num_rows(), false);
	}

	void check_weights(TPMatrix &P, size_t i)
	{
		IF_DEBUG(LIB_ALG_AMG, 2)
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
	double get_possible_parent_pairs(size_t i, FAMGNodes<TPMatrix> &rating)
	{
		AMG_PROFILE_NEXT(FAMG_getPPP_minCalc);
		int i_min = -1;
		double f_min = 1e12;


		neighborstruct_var<smallmat_ptype> s;
		s.parents.resize(2);

		set_sizes(2);
		//size_t blockSize = GetRows(H(0,0));
		double total_f_min=1e12;

#define IFI(s)
//#define IFI(s) if(i==0) { s; }
#define UG_LOGI(s) IFI(UG_LOG(s));

		IFI(PrintVector(onlyN1, "N1"))
		IFI(UG_LOG(CPPString(H, "H")));

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
				if (solve_KKT(i, two_indices(p0, p1)) != 0)
					continue;
				UG_LOGI("F = " << F << "\n");

				if ((F < 10 || F < 10*m_delta) && F < total_f_min)
					total_f_min = F;

				if (F > m_delta)
					continue;

				//check_one_interpolation(i, p0, p1);
				UG_DLOG(LIB_ALG_AMG, 2, "F: " << F << "\n");

				s.F = F;

				s.parents[0].value = q[0];
				s.parents[1].value = q[1];

				s.parents[0].from = gp0;
				s.parents[1].from = gp1;

				IF_DEBUG(LIB_ALG_AMG, 5)
				{
					UG_LOG("interpolating " << gp0 << " with " << q[0] << "\n");
					UG_LOG("interpolating " << gp1 << " with " << q[1] << "\n");
				}
				if (s.F < f_min)
				{
					f_min = s.F;
					i_min = i_neighborpairs.size();
				}

				//std::cout << s.parents[0].value << "   " << s.parents[1].value << " : " << s.F << "\n";

				i_neighborpairs.push_back(s);

				// F = q^T H q
			}
		}
		if(i_min != -1 && i_neighborpairs.size() > 0 && i_neighborpairs[0].F > i_neighborpairs[i_min].F)
			std::swap(i_neighborpairs[0],  i_neighborpairs[i_min]);
//		UG_LOG("best 2-F is " << total_f_min << "\n");
		return total_f_min != 1e12 ? total_f_min : 0.0;
	}
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
				for (size_t p2 = p1 + 1; p2 < onlyN1.size(); p2++)
				{
					size_t gp2 = onlyN1[p2];

//					UG_LOG("parents " << gp0 << ", " << gp1 << ", " << gp2 << "\n");
					if (solve_KKT(i, three_indices(p0, p1, p2)) != 0)
						continue;
					F*=1.10;

					if ((F < 10 || F < 10*m_delta) && F < total_f_min)
						total_f_min = F;

					if (F > m_delta)
						continue;

					//check_one_interpolation(i, p0, p1);
					UG_DLOG(LIB_ALG_AMG, 2, "F: " << F << "\n");

					s.F = F;

					s.parents[0].value = q[0];
					s.parents[1].value = q[1];
					s.parents[2].value = q[2];
					s.parents[0].from = gp0;
					s.parents[1].from = gp1;
					s.parents[2].from = gp2;

					IF_DEBUG(LIB_ALG_AMG, 5)

					IF_DEBUG(LIB_ALG_AMG, 5)
					{
						UG_LOG("interpolating " << gp0 << " with " << q[0] << "\n");
						UG_LOG("interpolating " << gp1 << " with " << q[1] << "\n");
					}
					if (s.F < f_min)
					{
						f_min = s.F;
						i_min = i_neighborpairs.size();
					}

					//std::cout << s.parents[0].value << "   " << s.parents[1].value << " : " << s.F << "\n";

					i_neighborpairs.push_back(s);
				}
				// F = q^T H q
			}
		}
		if(i_min != -1 && i_neighborpairs.size() > 0 && i_neighborpairs[0].F > i_neighborpairs[i_min].F)
			std::swap(i_neighborpairs[0],  i_neighborpairs[i_min]);
//		UG_LOG("best 3-F is " << total_f_min << "\n");
		return f_min;
	}
#endif

	bool prepare_interpolation_calculation(size_t i, FAMGNodes<TPMatrix> rating)
	{
		UG_DLOG(LIB_ALG_AMG, 2,
					"node " << i << " [" << rating.get_original_index(i) << "]\n");

		UG_ASSERT(rating[i].is_valid_rating(),
				"node " << i << " has no valid rating " << rating[i].rating);

		if (testvectorsExtern && BlockNorm2(m_testvectors[0][i]) < 1e-12)
		{
			UG_DLOG(LIB_ALG_AMG, 2,
					"testvector in " << i << " nearly 0, set as fine\n");
			rating.set_fine(i);
			return false;
		}

		// see local_eigenvalue_test_computation.h
		calc_special(i);

		GetNeighborhood(A_OL2, i, onlyN1, onlyN2, bvisited);
		if (onlyN1.size() == 0)
		{
			UG_DLOG(LIB_ALG_AMG, 2, "node has no connections, set as fine\n");
			rating.set_fine(i);
			return false;
		}

		// get H on onlyN1,i_index,onlyN2.
		get_H(i, rating);


		IF_DEBUG(LIB_ALG_AMG, 5)
		{
			for (size_t n = 0; n < onlyN1.size(); n++)
				UG_DLOG(LIB_ALG_AMG, 4,
						onlyN1[n] << " = " << rating.info(onlyN1[n]) << " ");
			UG_DLOG(LIB_ALG_AMG, 4, "\n");
		}
		return true;
	}
	// get_possible_parent_pairs:
	//---------------------------------------
	/** calculates of an index i a list of interpolating parent pairs.
	 * \param i						node index for which to calculate possible parent pairs
	 * \param possible_neighbors	here goes list of possible parent pairs
	 * \param rating				fine/coarse infos of the nodes
	 */
	double get_possible_parents(size_t i,
			stdvector<neighborstruct_var<smallmat_ptype> > &possible_neighbors, FAMGNodes<TPMatrix> &rating)
	{

		//if(i==0)
		{
//			UG_LOG("\n\n####################################\n\n");

//			UG_LOG("\n\n####################################\n\n");
		}

		AMG_PROFILE_FUNC();
		UG_DLOG(LIB_ALG_AMG, 2, "\n\n\n\n============================\n\n\n");
		if(prepare_interpolation_calculation(i, rating) == false) return 0;
//		UG_LOG("----" << i << "\n");

		i_neighborpairs.clear();

//		double bestF = get_all_neighbors_interpolation(i);
//		UG_LOG("Best F is " << bestF << "\n");
		get_possible_parent_pairs(i, rating);
		if(m_iMaxParents >= 3)
			get_possible_parent_threes(i, rating);

		if(i_neighborpairs.size() == 0)
		{
			UG_DLOG(LIB_ALG_AMG, 2,
			std::endl << rating.get_original_index(i) << ": UNINTERPOLATEABLE" << std::endl);
			//UG_ASSERT(0, "node has no parents :'-(");
			UG_ASSERT(rating.i_must_assign(i), i);
			rating.set_uninterpolateable(i);
			return 0.0;
		}
		else
			return extract_best_parents(i, possible_neighbors);

	}

	double extract_best_parents(size_t i, stdvector<neighborstruct_var<smallmat_ptype> > &possible_neighbors)
	{
		//AMG_PROFILE_NEXT(FAMG_getPPP_cleanup);
		possible_neighbors.clear();

		// minimal element is always first in possible_neighbors[i] list.
		double f_min = i_neighborpairs[0].F;

		for (size_t j = 0; j < i_neighborpairs.size(); j++)
			if (m_theta * i_neighborpairs[j].F <= f_min)
				possible_neighbors.push_back(i_neighborpairs[j]);

		//#ifdef FAMG_PRINT_POSSIBLE_PARENTS
/*		IF_DEBUG(LIB_ALG_AMG, 2)
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



	void set_sizes(size_t N)
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
		lq.resize(N+M);
		rhs.resize(N+M);
		KKT.resize(N+M, N+M);
#endif
	}

	// set KKT matrix
	/*
	 * solve_KKT
	 *
	 * \param i index of the fine node
	 * \param indices indices of the nodes to use out of onlyN1
	 *
	 * asserts H is calculated. \sa get_H
	 * \note since indices is templated, you can select fixed subsets
	 * like \sa two_indices , std::vector or \sa index_range.
	 *
	 */
#if SYS_AS_SYS
	template<typename TIndices>
	bool solve_KKT(size_t i, TIndices indices)
	{
		//AMG_PROFILE_FUNC();

		DenseVector<stdvector<smallvec_type> > &t0 = localTestvectors[0];

		const size_t N = indices.size();
		const size_t blockSize = GetRows(H(0, 0)); // assert blocksize is FIXED!!!

		/* we solve the following system. say i[]=indices[], n=indices.size()-1
		 *  ( H(i_0, i_0) ... H(i_0, i_n)  t[i_0]) ( q[i_0] )   ( h[i_0] )
		 *  (      ...            ...       ...  ) ( ...    )   (  ...   )
		 *  ( H(i_n, i_0) ... H(i_n, i_n)  t[i_n]) ( q[i_n] ) = ( h[i_n] )
		 *  ( t[i_0]^T    ... t[i_n]^T       0.0 ) ( lambda ) = ( t_i    )
		 */

		for (size_t r = 0; r < N; r++)
			for (size_t c = 0; c < N; c++)
				KKT.subassign(blockSize * r, blockSize * c,
						H(indices[r], indices[c]));

		for (size_t j = 0; j < N; j++)
		{
			KKT.subassign(blockSize * j, blockSize * N, t0[indices[j]]);
			KKT.subassign(blockSize * N, blockSize * j,
					te_transpose(t0[indices[j]]));
		}

		KKT(blockSize * N, blockSize * N) = 0;

		for (size_t j = 0; j < N; j++)
			rhs.subassign(blockSize * j, hat_h[indices[j]]);
		rhs[blockSize * N] = t_i;

		IF_DEBUG(LIB_ALG_AMG, 5)
			KKT.maple_print("KKT");
		IF_DEBUG(LIB_ALG_AMG, 5)
			rhs.maple_print("rhs");

		if (InverseMatMult(q, 1.0, KKT, rhs))
		{
			F = Hii - VecProd(rhs, q);
			F /= blockSize;
			return true;
		}
		else
		{
			UG_DLOG(LIB_ALG_AMG, 1, "solve_KKT: could not invert KKT :" << KKT);
			return false;
		}
	}
#else

	// Eliminate rows from non-strong couplings from KKT system
	template<typename TIndices>
	void reduceKKTSystem(TIndices indices, size_t b)
	{
		const size_t N = indices.size();
		std::vector<bool> strong(N, false);
		markAllStrong(indices, strong, m_dHReduceInterpolationNodesParameter,b);

		for(size_t r=0; r<N; r++)
		{
			// fix entries for non-strong
			if (strong[r]) continue;
			for(size_t c=0; c<=N; c++) KKT(r, c)= 0.0;
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
	bool check_smoother_sufficient(size_t b, size_t blockSize, size_t i_index)
	{
		Hii = BlockRef(H(i_index, i_index), b, b);
		const double Hiiscale = BlockRef(aii, b, b)*Hii;

		// quick check: do we really need to solve? (i.e., dirichlet component?)
		UG_DLOG(LIB_ALG_AMG, 5, "Hiiscale="<< Hiiscale << std::endl);

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
	int solve_KKT(size_t i, TIndices indices, bool bCalcNegatives=false)
	{
//		AMG_PROFILE_FUNC();

		const size_t blockSize=GetRows(H(0,0)); // assert blocksize is FIXED!!!
		const size_t N = indices.size();
		const size_t M = localTestvectors.size();
		size_t i_index = onlyN1.size();
		for(size_t k=0; k<indices.size(); k++)
			q[k] = 0.0;
		UG_DLOG(LIB_ALG_AMG, 3, "\n");
		F=0;

		bool bNegative=false;
		for(size_t b=0; b<blockSize; b++)
		{
			UG_LOGI("BLOCK [ " << b << " ]\n");

//			if(check_smoother_sufficient(b, blockSize, i_index) == true) continue;

			Hii = BlockRef(H(i_index, i_index), b, b);
			KKT = 0.0;

			// fill KKT matrix
			for(size_t r=0; r<N; r++)
				for(size_t c=0; c<N; c++)
					KKT(r, c)=BlockRef(H(indices[r], indices[c]), b, b);

			KKT(N, N) = 0;

			for(size_t k=0; k<localTestvectors.size(); k++)
				for(size_t j=0; j < N; j++)
					KKT(N+k, j) = KKT(j, N+k)= BlockRef(localTestvectors[k][indices[j]], b);


			UG_LOGI(JuliaString(KKT, "KKT"));
			for(size_t j=0; j<N; j++)
				rhs[j] = BlockRef(H(i_index, indices[j]), b, b);

			for(size_t k=0; k<M; k++)
				rhs[N+k] = BlockRef(localTestvectors[k][i_index],b);

			IF_DEBUG(LIB_ALG_AMG, 5) KKT.maple_print("KKT");
			IFI(rhs.maple_print("rhs"));
			IF_DEBUG(LIB_ALG_AMG, 5) KKT.maple_print("KKT1");
			IF_DEBUG(LIB_ALG_AMG, 5) rhs.maple_print("rhs1");


			// reduce to essential entries
			// todo: do less work here (KKT system is too large)

			reduceKKTSystem(indices, b);
			IF_DEBUG(LIB_ALG_AMG, 5) KKT.maple_print("KKT2");
			IF_DEBUG(LIB_ALG_AMG, 5) rhs.maple_print("rhs2");

			// solve system
			if(SolveDeficit(KKT, lq, rhs)==true)
			{
				IFI(KKT.maple_print("KKT_After"));
				IFI(rhs.maple_print("rhs_after"));
				IFI(lq.maple_print("lq"));

				IF_DEBUG(LIB_ALG_AMG, 5) lq.maple_print("lq");

				for(size_t k=0; k<indices.size(); k++)
				{
					if(lq[k] < 1e-8 && lq[k] > -1e-8)
						lq[k] = 0.0;
					if(lq[k] < 0)
					{
						IFI(UG_LOG("negative interpolation weights. aborting.\n"))
						IF_DEBUG(LIB_ALG_AMG, 2)
						{
							UG_LOG("WARNING: Negative entries at node " << i << ": ");
							for(size_t l=0; l<indices.size(); l++)
							{	UG_LOG("p" << l << " = " << onlyN1[ indices[l] ] << " "); }
							UG_LOG(" solve_KKT: lq[" << k << "] = " << lq[k] << "\n");
						}
						if(bCalcNegatives==false)
							return -1;
						bNegative=true;
					}
					BlockRef(q[k], b, b) = lq[k];
				}
				IF_DEBUG(LIB_ALG_AMG, 5) rhs.maple_print("rhs");
				//IF_DEBUG(LIB_ALG_AMG, 5)

				double lF=Hii;

				for(size_t k=0; k<N; k++)
					lF -= 2*lq[k]*BlockRef(H(i_index, indices[k]), b, b);
				for(size_t k1=0; k1<N; k1++)
					for(size_t k2=0; k2<N; k2++)
						lF += lq[k1]*BlockRef(H(indices[k1], indices[k2]), b, b)*lq[k2];
				lF *= BlockRef(aii, b, b);
				IFI(UG_LOGN(lF))
				F += lF;
//				UG_ASSERT(lF>-1e-10, "ERROR: negative values lF="<<lF<<" (b="<<b <<") are impossible!");
			}
			else
			{
				IFI(KKT.maple_print("KKT_After"));
				IFI(rhs.maple_print("rhs_after"));
				//IF_DEBUG(LIB_ALG_AMG, 2)
				{
//					UG_LOG("Problem at node " << i << ":");
//					for(size_t l=0; l<indices.size(); l++)
//					{	UG_LOG("p" << l << " = " << onlyN1[l] << " "); }

				}
				IFI(UG_LOG(" solve_KKT: could not solve KKT :" << KKT << "\n"));
				return -2;
			}
		}
		IFI(UG_LOG("SUCCESS, F = " << F << "\n"));
		// TODO: Sinvoll?
		F /= blockSize;

		return bNegative == false ? 0 : -1;
	}

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
#endif

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


	// select only coarse neighbors w/ sufficiently large couplings
	template<typename TIndices>
	void markAllStrong(const TIndices &coarseNeighbors,
				std::vector<bool> &strongMarks, double eps, int b=0)
	{

		UG_DLOG(LIB_ALG_AMG, 3, "markAllStrong with sigma = "<<eps <<" for comp=" << b << std::endl);
		const size_t i_index = onlyN1.size();
		strongMarks.resize(coarseNeighbors.size());

		// determine maximum of offdiag/diag in i-th row
		double maxVal = 0.0;
		for (size_t j = 0; j < coarseNeighbors.size(); j++)
		{
			const double val = fabs(BlockRef(H(i_index, coarseNeighbors[j]), b,b) / BlockRef(H(coarseNeighbors[j], coarseNeighbors[j]), b,b));
			maxVal = std::max(maxVal, val);
			UG_DLOG(LIB_ALG_AMG, 3, val<< " ")
			strongMarks[j]=false;
		}

		// eliminate coarse connections with small contributions
		maxVal *=eps;
		for (size_t j = 0; j < coarseNeighbors.size(); j++)
					if (fabs(BlockRef(H(i_index, coarseNeighbors[j]), b,b))
							>= maxVal * fabs(BlockRef(H(coarseNeighbors[j], coarseNeighbors[j]), b,b)))

						strongMarks[j]=true;

		IF_DEBUG(LIB_ALG_AMG, 3) {
			maple_print(coarseNeighbors, "coarseNghbs");
			maple_print(strongMarks,     "strongMarks");
		}

	}

	template<typename TIndices>
	std::vector<size_t> sortByStrength(const TIndices &coarseNeighbors, int b=-1)
	{
		const size_t i_index = onlyN1.size();
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
					strongMeasure[j] -= (fabs(BlockRef(H(i_index, coarseNeighbors[j]), bi,bi))/fabs(BlockRef(H(coarseNeighbors[j], coarseNeighbors[j]), bi,bi)) );
			else
				strongMeasure[j] = -(fabs(BlockRef(H(i_index, coarseNeighbors[j]), b,b))/fabs(BlockRef(H(coarseNeighbors[j], coarseNeighbors[j]), b,b)) );
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

	void getAllStrong(std::vector<size_t> &res, double eps)
	{
		std::vector<size_t> indices;
		for(size_t i=0; i<onlyN1.size(); i++) indices.push_back(i);
		std::vector<bool> strongMarksAll(indices.size(), false);
		std::vector<bool> strongMarks(indices.size(), false);
		size_t blockSize=GetRows(H(0,0));
		for(size_t i=0; i<blockSize; i++)
		{
			markAllStrong(indices, strongMarks, eps);
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

		UG_DLOG(LIB_ALG_AMG, 3, "reduceToStrong with sigma = "<<eps <<" for comp=" << b << std::endl);
		IF_DEBUG(LIB_ALG_AMG, 3) {
			maple_print(coarseNeighbors, "coarse");
		}
		const size_t i_index = onlyN1.size();

		// determine maximum of offdiag/diag in i-th row
		double maxVal = 0.0;
		for (size_t j = 0; j < coarseNeighbors.size(); j++)
		{
			const double val = BlockRef(H(i_index, coarseNeighbors[j]), b,b) / BlockRef(H(coarseNeighbors[j], coarseNeighbors[j]), b,b);
			maxVal = std::max(maxVal, val);
			UG_DLOG(LIB_ALG_AMG, 3, val<< " ")
		}

		// eliminate coarse connections with small contributions
		maxVal *=eps;
		for (size_t j = 0; j < coarseNeighbors.size(); j++)
		{
			if (BlockRef(H(i_index, coarseNeighbors[j]), b,b)
							>= maxVal * BlockRef(H(coarseNeighbors[j], coarseNeighbors[j]), b,b))
			{
				strongNeighbors.push_back(coarseNeighbors[j]);
			}
		}


		IF_DEBUG(LIB_ALG_AMG, 3) {
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


	bool direct_interpolation(size_t i, TPMatrix &P,
			FAMGNodes<TPMatrix> &rating, std::vector<size_t> &coarseNeighbors)
	{
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
		// global reduce
		/*
		 std::vector<size_t> strongNeighbors;
		reduceToStrong(coarseNeighbors, strongNeighbors, m_dHReduceInterpolationNodesParameter);
		std::vector<size_t> &myNeighbors = strongNeighbors;
        */

		std::vector<size_t> &myNeighbors = coarseNeighbors;
		set_sizes(myNeighbors.size());
		if(solve_KKT(i, myNeighbors) == 0)
		{
//			UG_ASSERT(F>-1e-10, "ERROR: negative F values are impossible (i.e., find bug and fix it!)");

			if(F > m_delta)
			{

				UG_DLOG(LIB_ALG_AMG, 3, "direct_interpolation, had to set node " << i << " coarse: " << F << ">" << m_delta);
							rating.external_set_coarse(i);
			}
			else
			{
				UG_DLOG(LIB_ALG_AMG, 3, "coarse neighbors, Interpolating from ");
				size_t N = myNeighbors.size();
				if(m_dProlongationTruncation > 0.0)
				{
					UG_ASSERT(0, "disabled due to systems");
				}
				if(m_dProlongationTruncation > 0.0)
				 {
					assert(0); // no thanks
					double dmax = BlockNorm(q[0]);
					for(size_t j=1; j<N; j++)
						dmax = std::max(BlockNorm(q[j]), dmax);
					for(size_t j=0; j<N; j++)
					{
						if(BlockNorm(q[j]) < dmax*m_dProlongationTruncation) continue;
						size_t node = onlyN1[myNeighbors[j]];
						UG_DLOG(LIB_ALG_AMG, 3, node << ": " << q[j] << ", ");
						P(i, node) = q[j];
					}
				 }
				else
				{
					for(size_t j=0; j<N; j++)
					{
						size_t nodej = onlyN1[myNeighbors[j]];
						P(i, nodej) = q[j];
						UG_DLOG(LIB_ALG_AMG, 3, nodej << ": " << q[j] << ", ");
					}
				}
				check_weights(P, i);
				if(fvalues.size() > i) fvalues[i] = F;
				rating.set_fine(i);
			}
		}
		else
		{
			rating.external_set_coarse(i);
			UG_DLOG(LIB_ALG_AMG, 3, "direct_interpolation of node " << i << ": could not invert KKT system (coarse neighbors).\n");
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
			//IF_DEBUG(LIB_ALG_AMG, 5)
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
			UG_DLOG(LIB_ALG_AMG, 3, "had to set node " << i << " coarse: ");
			if (F == 0)
			{
				UG_DLOG(LIB_ALG_AMG, 3, "could not invert KKT system.\n");
			}
			else
			{
				UG_DLOG(LIB_ALG_AMG, 3,
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
		GetNeighborhood(A_OL2, i, onlyN1, onlyN2, bvisited);

		for (size_t j = 0; j < onlyN1.size(); j++)
		{
			if (rating[onlyN1[j]].is_coarse())
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
		int ret = solve_KKT(i, ir, true);
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
		AMG_PROFILE_FUNC();
		UG_DLOG(LIB_ALG_AMG, 3,
				"get_all_neighbors_interpolation on node " << rating.get_original_index(i) << "\n")

		stdvector<size_t> coarseNeighbors;

		get_strong_connected_neighborhood2(A_OL2, rating, i, onlyN1, onlyN2,
				coarseNeighbors, bvisited);
		if (onlyN1.size() == 0)
		{
			UG_DLOG(LIB_ALG_AMG, 3, "no neighbors -> F");
			rating.set_fine(i);
			return true;
		}

		//const double &aii = A_OL2(i,i);

		/*if (coarseNeighbors.size() < 2)
			return indirect_interpolation(i, P, rating, coarseNeighbors);
		else if (indirectInterpolationLevel >= 0)
		{
			UG_DLOG(LIB_ALG_AMG, 3,
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

		IF_DEBUG(LIB_ALG_AMG, 5)
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
		IF_DEBUG(LIB_ALG_AMG, 1)
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

		IF_DEBUG(LIB_ALG_AMG, 5)
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

		IF_DEBUG(LIB_ALG_AMG, 5)
		{ S.maple_print("S_SF"); }
	}

	void calculate_H_from_S_and_Dinv()
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

		IF_DEBUG(LIB_ALG_AMG, 2)	{ H.maple_print("H"); }
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
			calculate_S_from_SF_and_Sjac();
			calculate_H_from_S_and_Dinv();
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

	const TMatrix &A;
	const TMatrix &A_OL2;
	double m_delta;
	double m_theta;
	double m_damping;
	double m_dProlongationTruncation;
	double m_dHReduceInterpolationNodesParameter;
	stdvector<TVector> &m_testvectors;
	stdvector<double> &m_omega;
	stdvector<double> &fvalues;
	const cAMG_helper &m_amghelper;

	bool testvectorsExtern;
	bool m_bHAsInverse;
	size_t m_iMaxParents;
};

// end group plugins_amg_famg
/// \}

} // namespace ug
#endif // __H__LIB_ALGEBRA__FAMG_SOLVER__FAMG_INTERPOLATION_CALCULATOR_H__
