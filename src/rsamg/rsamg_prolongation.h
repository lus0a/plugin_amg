/**
 * \file plugins/experimental/amg/src/rsamg/rsamg_prolongation.h
 *
 * \author Martin Rupp
 *
 * \date 06.08.2010
 *
 * Goethe-Center for Scientific Computing 2009-2010.
 */


#ifndef __H__UG__LIB_ALGEBRA__AMG_SOLVER__AMG_RS_PROLONGATION_H__
#define __H__UG__LIB_ALGEBRA__AMG_SOLVER__AMG_RS_PROLONGATION_H__

#include <vector>
#include "rsamg_nodeinfo.h"
#include "rsamg_blockvalues.h"
#include "rsamg_eliminate.h"
#include "../amg_misc.h"

//#define USE_DIRICHLET_AS_INTERPOLATION_NODES
#define SKIPDIRICHLET
namespace ug {

/// \addtogroup plugins_amg_rsamg
/// \{

template<typename TConnection>
bool CompareByValue(const TConnection &i, const TConnection &j)
{
	return dabs(i.dValue) > dabs(j.dValue);
}

template<typename TConnection>
bool CompareByIndex(const TConnection &i, const TConnection &j)
{
	return i.iIndex < j.iIndex;
}






// Truncates entries smaller than (eps x maximum value)
template<typename TProlongationMatrix>
bool TruncateInterpolation(std::vector<typename TProlongationMatrix::connection > &weight, double eps)
{
	// find maximum
	double sum = 0.0;
	double remove = 0.0;
	for(size_t j=0; j<weight.size(); j++)
	{
		sum += weight[j].dValue;
		remove = (weight[j].dValue>remove) ? weight[j].dValue : remove;
	}


	// determine nodes to be removed
	const double threshold=eps*remove;
	remove = 0.0;
	for(size_t j=0; j<weight.size(); j++)
	{
		sum += weight[j].dValue;
		if (weight[j].dValue <threshold) {

			remove += weight[j].dValue;
			weight[j].dValue= 0.0;
		}
	}

	// remove nodes
	if (fabs(remove) < 1e-8) return false;
	const double rescale = 1.0 + remove/(sum-remove);
	for(size_t j=0; j<weight.size(); j++)
	{
			weight[j].dValue *= rescale;
	}

	return true;

};

template<typename TProlongationMatrix>
void printVector (const char* title, size_t i, std::vector<typename TProlongationMatrix::connection > &v)
{
	std::cout << title << i << std::endl;
	for(size_t j=0; j<v.size(); j++)
		{
			std::cout<< "(" << v[j].iIndex << ", " << v[j].dValue << ") ";
		}
	std::cout <<std::endl;
}

template<typename TProlongationMatrix>
inline void SetRSInterpolation(TProlongationMatrix &P, size_t i,
		std::vector<typename TProlongationMatrix::connection > &conPositive, double sumPosNeighbors,
		std::vector<typename TProlongationMatrix::connection > &conNegative, double sumNegNeighbors,
		//double sumPosInterpolatory, double sumNegInterpolatory,
		double diag)
{
	// calculate alpha_i
	// set w_ij = alpha * w'_ij = alpha * a_ij/a_jj.

	// a) positive values
	const size_t numPositive =conPositive.size();
	if (numPositive > 0){
		// interpolation from positive
		printVector<TProlongationMatrix> ("SetRSInterpolation+ (PRE)", i, conPositive);

		double sumPosInterpolatory=0.0;
		for(size_t j=0; j<conPositive.size(); j++)
		{
			sumPosInterpolatory+=conPositive[j].dValue;
		}
		const double alphaPos = (sumPosNeighbors / sumPosInterpolatory) / diag;
		for(size_t j=0; j<conPositive.size(); j++)
		{
			conPositive[j].dValue *= alphaPos;
		}
		printVector<TProlongationMatrix> ("SetRSInterpolation+ (POST)", i, conPositive);
		TruncateInterpolation<TProlongationMatrix> (conPositive, 0.2);
		printVector<TProlongationMatrix> ("SetRSInterpolation+ (POST-TRUNCATE)", i, conPositive);

	} else {
		// no interpolation, thus lump to diagonal
		diag += sumPosNeighbors;
	}


	// b) negative values
	printVector<TProlongationMatrix> ("SetRSInterpolation- (PRE)", i, conNegative);
	double sumNegInterpolatory=0.0;
	for(size_t j=0; j<conNegative.size(); j++)
	{
		sumNegInterpolatory+=conNegative[j].dValue;
	}
	const double alphaNeg = (sumNegNeighbors / sumNegInterpolatory) / diag;
	for(size_t j=0; j<conNegative.size(); j++)
	{
		conNegative[j].dValue *= -alphaNeg;
	}
	printVector<TProlongationMatrix> ("SetRSInterpolation- (POST)", i, conNegative);
	TruncateInterpolation<TProlongationMatrix> (conNegative, 0.2);
	printVector<TProlongationMatrix> ("SetRSInterpolation- (POST-TRUNCATE)", i, conNegative);

	// set w_ij in matrix row P[i] forall j.
	P.set_matrix_row(i, &conNegative[0], conNegative.size());
	if (numPositive > 0) {
		P.add_matrix_row(i, &conPositive[0], conPositive.size());
	}

}





template<typename TProlongationMatrix, typename TMatrix>
inline void SetHensonInterpolation(TProlongationMatrix &P,
		const TMatrix &A, size_t i,
		std::vector<typename TProlongationMatrix::connection > &coarse,
		std::vector<typename TProlongationMatrix::connection > &strongf,
		double diag)
{

	// sort coarse vector w.r.t index
	sort(coarse.begin(), coarse.end(), CompareByIndex<SparseMatrix<double>::connection>);


	typedef typename TMatrix::const_row_iterator const_row_iterator;

	// \sum_{m \in C_i} a_km for all k in N_i
	for(size_t k=0; k<strongf.size(); k++)
	{
		bool found;

		strongf[k].dValue = 0;
		const size_t kindex = strongf[k].iIndex;
		typename TMatrix::value_type sumAkm=0.0;
		for (size_t m=0; m<coarse.size(); m++)
		{
			const_row_iterator connkm = A.get_connection(kindex, coarse[m].iIndex, found);

			if (found) sumAkm += connkm.value();
		}
		strongf[k].dValue = BlockRef(A(i,kindex),0,0) / BlockRef(sumAkm, 0, 0);

		// find kj-connection
		for(size_t j=0; j<coarse.size(); j++)
		{
			const_row_iterator connkj = A.get_connection(kindex, coarse[j].iIndex, found);  // ATTENTION: deg^2*log deg

			if (found) {
				coarse[j].dValue += strongf[k].dValue * BlockRef(connkj.value(),0,0);
			}
		}
	}

	// correct weight by diag ( + sum of weak offdiag)
	for(size_t j=0; j<coarse.size(); j++)
	{coarse[j].dValue /= (-1.0*diag);}

	// set w_ij in matrix row P[i] forall coarse j.
	P.set_matrix_row(i, &coarse[0], coarse.size());
}



inline void SetRSInterpolation_Sparsify(SparseMatrix<double> &P, size_t i,
		std::vector<SparseMatrix<double>::connection > &con,
		double sumPosNeighbors, double sumNegNeighbors,
		double sumPosInterpolatory, double sumNegInterpolatory,
		double diag)
{
	if(con.size()>=1000)
	{
		sort(con.begin(), con.end(), CompareByValue<SparseMatrix<double>::connection>);
		con.resize(6);
		for(size_t j=6; j<con.size(); j++)
		{
			if( con[i].dValue > 0)
				sumPosInterpolatory-=con[i].dValue;
			else
				sumNegInterpolatory-=con[i].dValue;
		}
	}

	// calculate alpha_i
	double alphaPos = - (sumPosNeighbors / sumPosInterpolatory) / diag;
	double alphaNeg = - (sumNegNeighbors / sumNegInterpolatory) / diag;
	// set w_ij = alpha * w'_ij = alpha * a_ij/a_jj.

	for(size_t j=0; j<con.size(); j++)
	{
		if(con[j].dValue < 0)
			con[j].dValue *= alphaNeg;
		else
			con[j].dValue *= alphaPos;
	}

	// set w_ij in matrix row P[i] forall j.
	P.set_matrix_row(i, &con[0], con.size());
}

template<typename TMatrix>
void GetNeighborValues(const TMatrix &A, size_t i, double &minConnValue, double &maxConnValue, double &diag)
{
	diag = 0; maxConnValue=0; minConnValue = 1e12;
	for(typename TMatrix::const_row_iterator conn = A.begin_row(i); conn != A.end_row(i); ++conn)
	{
		if(conn.index() == i)
		{
			diag = amg_diag_value(conn.value());
			continue; // skip diagonal
		}
#ifdef SKIPDIRICHLET
		if(A.is_isolated(conn.index()))
			continue;
#endif
		double connValue = amg_offdiag_value(conn.value());
		if(minConnValue > connValue)
			minConnValue = connValue;
		else if(maxConnValue < connValue)
			maxConnValue = connValue;
	}
}

template<typename TMatrix, typename TNodes>
void GetNeighborValues(const TMatrix &A, TNodes &nodes, size_t i, double &minConnValue, double &maxConnValue,
		double &minCoarseConnValue, double &maxCoarseConnValue, double &sumNegNeighbors, double &sumPosNeighbors, double &diag)
{
	maxConnValue = 0; minConnValue=1e12;
	maxCoarseConnValue = 0; minCoarseConnValue = 1e12;
	sumNegNeighbors = 0; sumPosNeighbors=0;
	diag=0;
	for(typename TMatrix::const_row_iterator conn = A.begin_row(i); conn != A.end_row(i); ++conn)
	{
		if(conn.index() == i)
		{
			diag = amg_diag_value(conn.value());
			continue;
		}
		double connValue = amg_offdiag_value(conn.value());

		if(connValue > 0)	sumPosNeighbors += connValue;
		else				sumNegNeighbors += connValue;

		if(maxConnValue < connValue) maxConnValue = connValue;
		if(minConnValue > connValue) minConnValue = connValue;

		if(nodes[conn.index()].is_coarse() || nodes[conn.index()].is_dirichlet())
		{
			if(maxCoarseConnValue > connValue)	maxCoarseConnValue = connValue;
			if(minCoarseConnValue < connValue)	minCoarseConnValue = connValue;
		}
	}
}
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// CreateProlongation:
//-------------------------
/**
 * Calculates Prolongation P with Matrix_type A and coarse/fine markers in
 * nodes[i].is_/is_coarse by direct interpolation
 * \param 	P				Matrix P: here goes the calculated prolongation
 * \param	A				Matrix A: matrix for which to calculate prolongation on next level
 * \param	nodes			fine/coarse marks of the nodes.
 * \param	theta			\f$\epsilon_{str}\f$.
 * \param   epsilonTruncation if != 0.0, truncates all entries from interpolation with p_ij < epsilonTruncation*max_j p_ij
 */


template<typename TProlongationMatrix, typename TMatrixType>
void CreateRugeStuebenProlongation(TProlongationMatrix &P, const TMatrixType &A,
		AMGNodes &nodes, const cgraph &strong,
		//stdvector<int> &newIndex,
		 double theta, double epsilonTruncation=0.0)
{
	AMG_PROFILE_FUNC();
	//UG_ASSERT(newIndex.size() == nodes.size(), "");
	P.resize_and_clear(A.num_rows(), A.num_rows());

	std::vector<typename TProlongationMatrix::connection> conPositive, conNegative;
	conPositive.reserve(255);
	conNegative.reserve(255);

	// DIRECT INTERPOLATION

	for(size_t i=0; i < A.num_rows(); i++)
	{
		if(!nodes.is_master_or_inner(i))
			continue;
		if(nodes[i].is_coarse())
		{
			// a coarse node
			//UG_ASSERT(newIndex[i] != -1, "coarse node but no new index?");
			//P(i, newIndex[i]) = 1.0;
			P(i,i) = 1.0;
		}
		else if(A.is_isolated(i))
		{
			//P[i].initWithoutDiag(); // boundary values need not to be prolongated
		}
		else if(nodes[i].is_fine_direct())
		{
			// a non-interpolated fine node. calculate interpolation weights

			// calc min off-diag-entry, and sum of Neighbors
			double maxConnValue, minConnValue, diag;
			GetNeighborValues(A, i, minConnValue, maxConnValue, diag);

			// todo: check if it is ok to do it THIS way:
			/*if(epsilonTruncation > 0)  // [AMGKS99] 7.2.4 truncation of interpolation
				barrier = std::min(theta*dmax, epsilonTruncation*maxCoarseConnValue);
			else*/
				//barrier = theta*dmax;
		//	double barrier = theta*std::max(maxConnValue, dabs(minConnValue));





			// Eliminate connections
			TMatrixType Bi;
			Bi.resize_and_clear(A.num_rows(),A.num_cols());
			//StrongOffDiagCheck<typename TMatrixType::value_type> rule(strong);
			PositiveOffDiagCheck<typename TMatrixType::value_type> positive;
			EliminateEntries(positive, A, i, Bi, true);

			//double sumPosInterpolatory=0, sumNegInterpolatory=0;
			double sumPosNeighbors=0.0, sumNegNeighbors=0.0;
			conPositive.clear();
			conNegative.clear();

			// step 1: set w'_ij = a_ij/a_jj for suitable j
			const TMatrixType &C = Bi;
			for(typename TMatrixType::const_row_iterator conn = C.begin_row(i); conn != C.end_row(i); ++conn)
			{
				const size_t j = conn.index() ;
				if(j == i) continue; // skip diagonal

#ifdef SKIPDIRICHLET
				if(A.is_isolated(j))
					continue;
#endif

				double connValue = amg_offdiag_value(conn.value());
				if(connValue > 0) 	sumPosNeighbors += connValue;
				else				sumNegNeighbors += connValue;


				if (!strong.has_connection(i,j))
						continue;
				//if(dabs(connValue) < barrier)
				//	continue;

				// only strong neigbors here
				if(!nodes[j].is_coarse()
#ifdef USE_DIRICHLET_AS_INTERPOLATION_NODES
					&& !nodes[j].is_dirichlet()
#endif
					) continue;

				if(nodes[j].is_dirichlet()) continue;

				// only coarse nodes here
				typename TProlongationMatrix::connection c;
				c.iIndex = j;
				c.dValue = connValue;

				UG_ASSERT(c.iIndex >= 0, "not coarse?");

				// add to interpolatory set and update sums

				// 	con.push_back(c);
				if(connValue > 0)	{conPositive.push_back(c);} //sumPosInterpolatory += connValue;
				else				{conNegative.push_back(c);} //sumNegInterpolatory += connValue;
			}

			if((conPositive.size() + conNegative.size()) > 0)
				SetRSInterpolation(P, i,
						conPositive, sumPosNeighbors,
						conNegative, sumNegNeighbors,
						//sumPosInterpolatory, sumNegInterpolatory,
						diag);
			else
			{
				IF_DEBUG(DebugRSAMGProlongation, 2)
				{
					UG_LOG("node " << i << " was marked direct, but no interpolation found." << std::endl);// barrier = " << barrier <<
					UG_LOG("[\n");
					for(typename TMatrixType::const_row_iterator conn = A.begin_row(i); conn != A.end_row(i); ++conn)
					{
						if(conn.index() == i) continue; // skip diagonal
						if(!nodes[conn.index()].is_coarse())
						{	UG_LOG(" " << conn.index() << " is not coarse\n"); }
						else
						{	UG_LOG(" " << conn.index() << " has value " << amg_offdiag_value(conn.value()) << "\n"); }
					}
					UG_LOG("]\n");
				}
				//UG_ASSERT(0,i);
				// no suitable interpolating nodes for node i,
				// so this node has to be treated by indirect interpolation
				//nodes.set_unassigned_fine_indirect(i);
				nodes.set_coarse(i);
			}
		}
		/*else
		{
			UG_ASSERT(0, "?");
			//unassigned++;
			//UG_ASSERT(aggressiveCoarsening != 0, "no aggressive Coarsening but node " << i << " is fine and indirect??");
		}*/
	}

	if(nodes.get_unassigned())
		UG_DLOG(DebugRSAMGProlongation, 1, nodes.get_unassigned() << " ? ")
	if(nodes.get_unassigned_indirect_fine())
		UG_DLOG(DebugRSAMGProlongation, 1, "Pass 1: " << nodes.get_unassigned_indirect_fine() << " nodes need indirect interpolation.")
}


template<typename TProlongationMatrix, typename Matrix_type>
void CreateHensonProlongation(TProlongationMatrix &P, const Matrix_type &A, AMGNodes &nodes,
		//stdvector<int> &newIndex,
		const cgraph &strong, double theta, double epsilonTruncation=0.0)
{
	AMG_PROFILE_FUNC();
	P.resize(A.num_rows(), A.num_rows());

	std::vector<typename TProlongationMatrix::connection> coarse(255);
	std::vector<typename TProlongationMatrix::connection> strongf(255);

	// DIRECT INTERPOLATION
	for(size_t i=0; i < A.num_rows(); i++)
	{
		if(!nodes.is_master_or_inner(i)) continue;

		if(nodes[i].is_coarse())
		{
			// a coarse node
			P(i,i) = 1.0;
		}
		else if(A.is_isolated(i))
		{
			//P[i].initWithoutDiag(); // boundary values need not to be prolongated
		}
		else if(nodes[i].is_fine_direct())
		{
			// a non-interpolated fine node. calculate interpolation weights
			coarse.clear();
			strongf.clear();


			// step 1: set w'_ij = a_ij/a_jj for suitable j

			double diag = 0.0;
			for(typename Matrix_type::const_row_iterator conn = A.begin_row(i); conn != A.end_row(i); ++conn)
			{
				const size_t j = conn.index();

#ifdef SKIPDIRICHLET
				if(A.is_isolated(j)) continue;
#endif
				if( (j==i) || (!strong.has_connection(i,j)) ) {
					// skip diagonal and lump weak neighbors to diagonal
					diag += conn.value();
					continue;
				}

				typename TProlongationMatrix::connection c;
				c.iIndex = j;
				c.dValue = amg_offdiag_value(conn.value());

				// only strong neigbors here
				if(!nodes[j].is_coarse()
#ifdef USE_DIRICHLET_AS_INTERPOLATION_NODES
						&& !nodes[j].is_dirichlet()
#endif
				)
				{
					if(!nodes[j].is_dirichlet()) strongf.push_back(c);
					continue;
				}

				// add coarse nodes to interpolatory set
				if(!nodes[j].is_dirichlet()) coarse.push_back(c);
				UG_ASSERT(j >= 0, "not coarse?");
			}


			if(coarse.size() > 0)
				SetHensonInterpolation(P, A, i, coarse, strongf, diag);
			else
			{
				IF_DEBUG(DebugRSAMGProlongation, 2)
				{
				 /*UG_LOG("node " << i << " was marked direct, but no interpolation found. barrier = " << barrier << std::endl);
					UG_LOG("[\n");
					for(typename Matrix_type::const_row_iterator conn = A.begin_row(i); conn != A.end_row(i); ++conn)
					{
						if(conn.index() == i) continue; // skip diagonal
						if(!nodes[conn.index()].is_coarse())
						{	UG_LOG(" " << conn.index() << " is not coarse\n"); }
						else
						{	UG_LOG(" " << conn.index() << " has value " << amg_offdiag_value(conn.value()) << "\n"); }
					}
					UG_LOG("]\n");*/
				}
				//UG_ASSERT(0,i);
				// no suitable interpolating nodes for node i,
				// so this node has to be treated by indirect interpolation
				//nodes.set_unassigned_fine_indirect(i);
				nodes.set_coarse(i);
			}
		}
		/*else
		{
			UG_ASSERT(0, "?");
			//unassigned++;
			//UG_ASSERT(aggressiveCoarsening != 0, "no aggressive Coarsening but node " << i << " is fine and indirect??");
		}*/
	}

	if(nodes.get_unassigned())
		UG_DLOG(DebugRSAMGProlongation, 1, nodes.get_unassigned() << " ? ")
	if(nodes.get_unassigned_indirect_fine())
		UG_DLOG(DebugRSAMGProlongation, 1, "Pass 1: " << nodes.get_unassigned_indirect_fine() << " nodes need indirect interpolation.")
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// CreateIndirectProlongation:
//-------------------------
/**
 * Assume Prolongation of all normal fine nodes is already computed, it calculates the Interpolation of
 * fineIndirect nodes with Matrix_type A and Coarse/FineIndirect markers in nodes[i].is_coarse/is_fine_indirect
 *
 * [AMGKS99] 7.2.2 Multi-pass-interpolation
 * Using A1 or A2 aggressive coarsening, it can be shown that we need at most 4 passes.
 * \note CreateIndirectProlongation is starts with pass 2, since pass 1 is CreateRugeStuebenProlongation.
 *
 * \param 	P					Matrix P: here goes the calculated prolongation
 * \param	A					Matrix A: matrix for which to calculate prolongation on next level
 * \param	nodes
 * \param 	posInConnections	array of size A.num_rows() for speedup of neighbor-neighbor-calculation inited with -1.
 * \param	theta
 */
template<typename TProlongationMatrix, typename Matrix_type>
void CreateIndirectProlongation(TProlongationMatrix &P, const Matrix_type &A,
		AMGNodes &nodes, int *posInConnections, cgraph &strong)
{
	AMG_PROFILE_FUNC();
	std::vector<typename TProlongationMatrix::connection> conPositive, conNegative, con;
	conPositive.reserve(255);
	conNegative.reserve(255);
	con.reserve(255);

	std::vector<int> nrOfPaths;
	nrOfPaths.reserve(255);
	typename TProlongationMatrix::connection c;
	//P.print();
	// INDIRECT INTERPOLATION

	size_t oldUnassigned = -1;
	int pass=2;
	while(nodes.get_unassigned_indirect_fine())
	{
#ifdef AMG_PRINT_INDIRECT
		UG_DLOG(LIB_ALG_AMG, 1, std::endl);
#endif
		UG_DLOG(LIB_ALG_AMG, 1, "Pass " << pass << ": ");
		for(size_t i=0; i<A.num_rows() && nodes.get_unassigned_indirect_fine() > 0; i++)
		{
			if(!nodes[i].is_unassigned_fine_indirect() || A.is_isolated(i))
				continue;

			// calculate min offdiag-entry
			double diag, maxConnValue, minConnValue;
			GetNeighborValues(A, i, minConnValue, maxConnValue, diag);

			double sumPosNeighbors=0.0, sumNegNeighbors=0.0;
			conPositive.clear();
			conNegative.clear();

			//double barrier = theta*std::max(maxConnValue, dabs(minConnValue));

			//cout << "indirect interpolating node " << i << endl;

			// look at neighbors of node i, try to interpolate indirectly through them
			for(typename Matrix_type::const_row_iterator conn = A.begin_row(i); conn != A.end_row(i); ++conn)
			{
				size_t indexN = conn.index();
				if(indexN == i) continue; // skip diagonal

				double connValue = amg_offdiag_value(conn.value());
				if(connValue > 0) 	sumPosNeighbors += connValue;
				else				sumNegNeighbors += connValue;

				if (!strong.has_connection(i,indexN))
						continue;

				// we dont want fine nodes which were indirectly interpolated in THIS pass
				// (no gauss-seidel-style)
				if(nodes[indexN].is_fine_indirect_level(pass))
					continue;
				// all interpolate neighbors are now from pass (pass-1)


				UG_ASSERT(!nodes[indexN].is_coarse(), "Node " << i << " indirect, but neighbor " <<  indexN << " coarse?");

				// now we look at the interpolation values of this neighbor (in P !!!)
				// and set w'(i,indexNN) += A(i,indexN) * P(indexN,indexNN)
				for(typename TProlongationMatrix::row_iterator conn2 = P.begin_row(indexN); conn2 != P.end_row(indexN); ++conn2)
				{
					size_t indexNN = conn2.index();
					int pos = posInConnections[indexNN];

					if(pos == -1)
					{
						pos = posInConnections[indexNN] = con.size();
						c.iIndex = indexNN; assert(c.iIndex >= 0);

						AssignMult(c.dValue, connValue, conn2.value());
						con.push_back(c);
					}
					else
						AddMult(con[pos].dValue, connValue, conn2.value());
				}
			}

			for(size_t j=0; j<con.size(); j++)
			{
				if(con[j].dValue > 0)
					conPositive.push_back(con[j]);
				else
					conNegative.push_back(con[j]);
				posInConnections[con[j].iIndex] = -1; 	// reset posInConnections
			}

			if(con.size() == 0)
				continue;

			nodes.set_fine_indirect_level(i, pass);

#ifdef AMG_PRINT_INDIRECT
			cout << i << " ";
#endif
			//cout << endl;

			SetRSInterpolation(P, i,
					conPositive, sumPosNeighbors,
					conNegative, sumNegNeighbors,
					diag);
		}


		if(nodes.get_unassigned_indirect_fine() == oldUnassigned)
		{
			UG_LOG(std::endl << "unassigned nodes left: " << std::endl);
			for(size_t i=0; i<A.num_rows(); i++)
			{
				if(nodes[i].is_unassigned_fine_indirect())
				 UG_LOG(i << " ");
			}
			UG_LOG("\n");
		}
		UG_ASSERT(nodes.get_unassigned_indirect_fine() != oldUnassigned, "Pass " << pass <<
				": Indirect Interpolation hangs at " << nodes.get_unassigned_indirect_fine() <<
				" unassigned indirect fine nodes.");

#ifdef AMG_PRINT_INDIRECT
		UG_DLOG(LIB_ALG_AMG, 1, "calculated, ");
#endif
		UG_DLOG(LIB_ALG_AMG, 1, nodes.get_unassigned_indirect_fine() << " left. ");
		pass++;
		oldUnassigned = nodes.get_unassigned_indirect_fine();
		break;
	}

	UG_ASSERT(nodes.get_unassigned_indirect_fine() == 0, "number of unassigned nodes is still " <<
			nodes.get_unassigned_indirect_fine());

	P.defragment();
	//P.print();
}
// end group plugins_amg_rsamg
/// \}

} // namespace ug

#endif /* __H__UG__LIB_ALGEBRA__AMG_SOLVER__AMG_RS_PROLONGATION_H__ */
