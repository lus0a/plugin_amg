/*
 * famg_external_coarsening.h
 *
 *  Created on: 21.11.2011
 *      Author: mrupp
 */

#include <fstream>

#include "../rsamg/rsamg.h"
#include "../rsamg/rsamg_coarsening.h"

#include "lib_algebra/adapter/scalar_submatrix_adapter.hh"

#include "common/error.h"
#include <boost/lambda/lambda.hpp>

#ifdef UG_PARALLEL
#include "array_binary_operation_scheme.h"
#endif

#ifndef FAMG_EXTERNAL_COARSENING_H_
#define FAMG_EXTERNAL_COARSENING_H_

namespace ug{

/// \addtogroup plugins_amg_famg
/// \{




/**
 * calculate the FAMG interpolation to a given coarse/fine marking
 */
template<typename algebra_type>
void FAMGLevelCalculator<algebra_type>::external_coarsening_calculate_prolongation()
{
	try{
	AMG_PROFILE_FUNC();

	Stopwatch SW;
	if(m_bTiming) { UG_LOG("\ncalculate prolongation... "); SW.start(); }
	size_t N = A.num_rows();
	PROGRESS_START(prog, N, "FAMG::external_coarsening_calculate_prolongation " << N);

//	UG_LINE_LOG();
//	rating.print();
	for(size_t i=0; i<N; i++)
	{
		PROGRESS_UPDATE(prog, i);
		if(rating[i].is_fine_direct() && rating.i_must_assign(i))
		{
			UG_DLOG(DBG_FAMG_INTERPOLATION, 2, "calculator.get_all_neighbors_interpolation " << i << ", " << rating[i] << "\n");
			calculator.get_all_neighbors_interpolation(i, PoldIndices, rating, false);
		}
	}
	PROGRESS_FINISH(prog)


	if(m_famg.m_bIndirectInterpolation)
	{
		PROGRESS_START(prog2, N, "FAMG::external_coarsening_calculate_prolongation, indirect N = " << N);
		for(size_t i=0; i<N; i++)
		{
			PROGRESS_UPDATE(prog2, i);
			if(rating[i].is_uninterpolateable() && rating.i_must_assign(i))
			{
				UG_DLOG(DBG_FAMG_INTERPOLATION, 2, "calculator.get_all_neighbors_interpolation indirect " << i << ", " << rating[i] << "\n");
				calculator.get_all_neighbors_interpolation(i, PoldIndices, rating, true);
			}
		}
	}
	if(m_bTiming) UG_LOG("took " << SW.ms() << " ms");

	}UG_CATCH_THROW_FUNC();
}

/**
 * transfer external coarsening in RS AMGNodes to FAMGNodes.
 */
template<typename T>
inline void AMGnodeinfoToFAMGnodeinfo(AMGNodes &amgNodes, FAMGNodes<T> &famgNodes, cgraph &graphS)
{
	AMG_PROFILE_FUNC();
	size_t N = amgNodes.size();
	for(size_t i=0; i<N; i++)
	{
		if(graphS.is_isolated(i) || amgNodes[i].is_dirichlet())
			famgNodes.set_dirichlet(i);
		else if(amgNodes[i].is_coarse())
			famgNodes.external_set_coarse(i);
		else if(amgNodes[i].is_fine_direct() || amgNodes[i].is_fine())
			famgNodes.set_fine(i);
		else
		{

			famgNodes.set_fine(i);
	//		UG_LOG("nodes " << i << " = " << famgNodes[i] << "\n");
//			UG_LOG("rating " << i << " = " << amgNodes[i] << "\n");
		}
	}

}

/**
 * calculate the external coarsening graph from A_OL2.
 * @param graphS
 */
template<typename algebra_type>
void FAMGLevelCalculator<algebra_type>::get_external_coarsening_graph(cgraph &graphS)
{
	try{
	AMG_PROFILE_FUNC();
	const int submask=m_famg.get_strong_connection_submask();
	const double sigma = m_famg.get_strong_connection_external();

	switch (submask) {
	case 0: // primary sub 0
		{ScalarSubMatrixAdapter<algebra_type, CPUAlgebra, 0, 0>
		sub(const_cast<typename algebra_type::matrix_type&> (A_OL2)); // warning: slow
		CreateStrongConnectionGraph(sub, graphS, sigma);
		break;}

	case 1: // primary sub 1
		{ScalarSubMatrixAdapter<algebra_type, CPUAlgebra, 1, 1>
		sub(const_cast<typename algebra_type::matrix_type&> (A_OL2)); // warning: slow
		CreateStrongConnectionGraph(sub, graphS, sigma);
		break;}

	case 2: // primary sub 2
		{ScalarSubMatrixAdapter<algebra_type, CPUAlgebra, 2, 2>
		sub(const_cast<typename algebra_type::matrix_type&> (A_OL2)); // warning: slow
		CreateStrongConnectionGraph(sub, graphS, sigma);
		break;}

	case -2: //full primary matrix (norm based)
	{// todo: could be replaced by primary matrix!
		CreateStrongConnectionGraph(A_OL2, graphS, sigma);
		break;}

	default:
		// standard scalar case
		CreateStrongConnectionGraph(A_OL2, graphS, sigma);
	}

	}UG_CATCH_THROW_FUNC();
}

static void MyCreateMeasureOfImportancePQ(const cgraph &strong, const cgraph &strongT, nodeinfo_pq_type &PQ, AMGNodes &nodes)
{
	AMG_PROFILE_FUNC();
	UG_ASSERT(nodes.size() == strongT.size(), "");

	PQ.create(nodes.size(), &nodes[0]);
	for(size_t i=0; i < nodes.size(); i++)
	{
		if(strong.is_isolated(i))
			nodes.set_dirichlet(i);
		else if(nodes.is_master_or_inner(i) == false)
			nodes.set_parallel_dont_care(i);
		else if(!nodes[i].is_coarse())
		{
			//UG_ASSERT(graph.iNrOfConnections[i] > 0, "node " << i << " has " << graph.iNrOfConnections[i] << " connections?");
			nodes.set_rating(i, strongT.num_connections(i));
			PQ.insert_item(i);
		}
	}
}

/**
 * make FAMG coarsening consistent by Communicate coarsening information form ol2master -> ol2slave.
 */
template<typename algebra_type>
void FAMGLevelCalculator<algebra_type>::make_coarsening_consistent()
{
#ifdef UG_PARALLEL
	AMG_PROFILE_FUNC();
	StdArrayCommunicationScheme<FAMGNodes<prolongation_matrix_type> > scheme2(rating);
	CommunicateOnInterfaces(PN.comm(), OL2MasterLayout, OL2SlaveLayout, scheme2);


	ConsistencyCheckClass<FAMGNodes<prolongation_matrix_type>, typename FAMGNodes<prolongation_matrix_type>::value_type> scheme(rating, "amgnodes");
	CommunicateOnInterfaces(PN.comm(), PN.master_layout(), PN.slave_layout(), scheme);

	if(!AllProcsTrue(scheme.isOK(), PN.proc_comm()))
	{
		PRINTLAYOUT(PN.proc_comm(), PN.comm(), PN.master_layout(), PN.slave_layout());
		rating.print();
		UG_COND_THROW(!AllProcsTrue(scheme.isOK(), PN.proc_comm()), "amgnodes" << " not consistent!");
	}
#endif

}

/*
template<typename TPMatrix>
class FAMGSendCoarseNodes :
	public CommunicationScheme<FAMGCoarseningCommunicationScheme<TPMatrix>, bool>
{
public:
	FAMGCommunicateCoarsening(FAMGNodes<TPMatrix> &r) : rating(r)
	{

	}
	bool send(int pid, size_t index) const
	{
		int state = rating[index].get_state();
		if(state == FAMG_COARSE_RATING)
			return 1;
		else
			return 0;
	}
	void receive(int pid, size_t index, bool state)
	{
		if(state == 0)
		{
			//rating.external_set_uncalculated_fine(index);
			rating.set_fine(index);
		}
		else if(state == 1)
			rating.external_set_coarse(index);
		else if(state == 3)
		{
			rating.set_dirichlet(index);
		}
		// else ; // nothing

	}

	inline int get_element_size() const
	{
		return sizeof(int);
	}
private:
	FAMGNodes<TPMatrix> &rating;
};
*/

/**
 * will perform a post-coarsening step on all uninterpolateable nodes.
 * ex 2d neumann:
 *
 * U U U U U U U  (U=uninterpolatable)
 * U C F C F C U
 * U F F F F F U
 * U C F C F C U
 * U U U U U U U
 *
 */
template<typename algebra_type>
void FAMGLevelCalculator<algebra_type>::rs_amg_external_postcoarsening()
{

	AMG_PROFILE_FUNC();
	Stopwatch SW;
	if(m_bTiming) { UG_LOG("\nexternal post-coarsening... " );  SW.start(); }
	nodeinfo_pq_type PQ;
	size_t N = A_OL2.num_rows();
	//UG_LOG("RATING:\n"); rating.print();

	// use RSAMG's methods to do standard rs coarsening
#ifdef UG_PARALLEL
	AMGNodes nodes(N, PN);
	UG_COND_THROW(m_famg.m_pParallelCoarsening == NULL, "please set_parallel_coarsening !");
	UG_DLOG(DBG_AMG_COARSENING, 0, m_famg.m_pParallelCoarsening->config_string() << "\n");
#else
	AMGNodes nodes(N);
#endif
	cgraph graphS, graphST, graphS0;
	get_external_coarsening_graph(graphS0);
	SetIsolatedCoarse(graphS0, nodes);


	//UG_COND_THROW(pcl::ProcRank() > 1, "communicate uninterpolateable nodes");

	graphS.resize(graphS0.size());
	// only consider connections between uninterpolateable nodes.
	for(size_t i=0; i<nodes.size(); i++)
		if(rating[i].is_uninterpolateable())
		{
			for(cgraph::row_iterator it = graphS0.begin_row(i); it != graphS0.end_row(i); ++it)
			{
				if(rating[*it].is_uninterpolateable())
					graphS.set_connection(i, *it);
			}
		}
		else
			nodes.set_coarse(i);

	graphST.set_as_transpose_of(graphS);


	MyCreateMeasureOfImportancePQ(graphS, graphST, PQ, nodes);

	for(size_t i=0; i<nodes.size(); i++)
	{
		if(rating[i].is_uninterpolateable() && !nodes[i].is_dirichlet())
		{
			int iR=10;
			for(cgraph::row_iterator it = graphS0.begin_row(i); it != graphS0.end_row(i); ++it)
			{
				if(rating[*it].is_fine())
					iR+=3;
				//if(rating[*it].is_coarse())
					//nodes[i].rating ++;

			}
			nodes[i].rating += iR;
			PQ.update(i);
		}
	}

#ifdef UG_PARALLEL
	m_famg.m_pParallelCoarsening->coarsen(PN, masterLayouts, slaveLayouts, graphS, graphST, PQ, nodes, true, false);
#else
	Coarsen(graphS, graphST, PQ, nodes, true, false);

	if (m_famg.get_prevent_ff_connection_external())
	{
		PreventFFConnections(graphS, graphST, nodes);
	}
#endif
//TODO: Einfluss von fine connected nodes auf Rating!
	// coarsening done. transfer AMG nodeinfo -> FAMG nodeinfo
	std::vector<bool> wasUninterpolatable(nodes.size(), false);

	for(size_t i=0; i<nodes.size(); i++)
		if(rating[i].is_uninterpolateable())
		{
			wasUninterpolatable[i] = true;
			if(nodes[i].is_fine_direct() || nodes[i].is_fine())
				rating.set_fine(i);
			else
				rating.set_coarse(i);
		}

	//----------------------
	if(m_famg.m_writeMatrices)
		write_XC_coarsening(graphS, nodes, "_postC");

	for(size_t i=0; i<N; i++)
	{
		if(wasUninterpolatable[i] && rating[i].is_fine_direct() && rating.i_must_assign(i))
			calculator.get_all_neighbors_interpolation(i, PoldIndices, rating, false);
	}


	// if there are still coarse nodes left, which are not used in the interpolation, try to interpolate them again because
	// there might be more coarse nodes now

	typedef typename prolongation_matrix_type::row_iterator Prow_iterator;
	/*std::vector<bool> interpolatedFrom(N, false);
	for(size_t i=0; i<N; i++)
	{

		if(nodes[i].is_coarse()) continue;
		for(Prow_iterator it = PoldIndices.begin_row(i); it != PoldIndices.end_row(i); ++it)
			interpolatedFrom[it.index()] = true;
	}

	//return;
	SendVecFromTo(PN.comm(), PN.slave_layout(), PN.master_layout(), interpolatedFrom, boost::lambda::_1 || boost::lambda::_2);*/
	for(size_t i=0; i<N; i++)
	{
		//if(interpolatedFrom[i] == false &&
				/*rating[i].is_coarse() ||*/
		if(rating[i].is_uninterpolateable()

				/*&& rating.i_must_assign(i))*/
				&& rating.is_inner_node(i))
		{
			rating.node_reset(i);
			calculator.get_all_neighbors_interpolation(i, PoldIndices, rating, false);
			//UG_LOG("postpost " << i << ", now " << rating[i] << ".");
			if(rating[i].is_fine())
			{
				//UG_LOG(". interpolating from ")
				// adjust interpolatedFrom
				/*
				for(Prow_iterator it = PoldIndices.begin_row(i); it != PoldIndices.end_row(i); ++it)
				{

					//UG_LOG(it.index() << " ");
					interpolatedFrom[it.index()] = true;
				}
				*/
//				UG_LOG("\n");
			}
			else
			{
//				UG_LOG(" not possible.\n");
			}
		}
		else
		{
//			UG_LOG("postpost " << i << ", now " << rating[i] << ", interpolatedFrom = " << TrueFalseString(interpolatedFrom[i]) << "\n");
			}
	}


//	set_uninterpolateable_as_coarse();
}

#ifdef UG_PARALLEL
/**
 * Griebel LE strategy for neumann-coarsening
 * the idea is to detect the neumann boundary by comparing the number of connections of node i with the avg. number of connections of i's neighbors
 * if node i has 'under-average' number of connections, it is considered to be a neumann node and gets a higher rating, so that it is coarsened first.
 * @param graphS
 * @param graphST
 * @param PQ
 * @param nodes
 */
template<typename algebra_type>
void FAMGLevelCalculator<algebra_type>::rs_amg_coarsening_LE_griebel(const cgraph &graphS, const cgraph &graphST, nodeinfo_pq_type &PQ, AMGNodes &nodes)
{
	AMG_PROFILE_FUNC();
	const matrix_type &M = A_OL2;
	typedef typename matrix_type::const_row_iterator const_row_iterator;
	for(size_t i=0; i<A.num_rows(); i++)
	{
		if(nodes[i].is_assigned()) continue;

		size_t sum = 0;
		for(const_row_iterator it = M.begin_row(i); it != M.end_row(i); ++it)
			sum += A_OL2.num_connections(it.index());
		size_t Ni = A_OL2.num_connections(i);
		if(Ni*Ni < sum)
		{
			nodes[i].rating += 10;
			PQ.update(i);
		}
	}
}

/**
 * This algorithm tries to favor nodes to be coarse when they don't have much possible interpolation
 * options. These nodes often get coarse anyways, so we set them coarse so that other nodes can be fine.
 * example: 2d neumann boundary:
 * std coarsening
 * F F F F F F F
 * F C F C F C F
 * F F F F F F F
 * F C F C F C F
 * F F F F F F F
 * will be modified to
 * U U U U U U U  (U=uninterpolatable)
 * U C F C F C U
 * U F F F F F U
 * U C F C F C U
 * U U U U U U U
 *
 * now U can be set to coarse (\sa set_uninterpolateable_coarse)
 * or we can do some post-coarsening on the set of uninterpolated nodes (\sa rs_amg_external_postcoarsening)
 *
 * prefered coarsening
 * C F C F C F C
 * F F F F F F F
 * C F C F C F C
 * F F F F F F F
 * C F C F C F C
 *
 * nrOfInterpolateableFrom
 *
 * 0 1 1 1 1 1 0
 * 1 4 4 4 4 4 1
 * 1 4 4 4 4 4 1
 * 1 4 4 4 4 4 1
 * 0 1 1 1 1 1 0
 *
 * pre-coarsening
 *
 * C F C F C F C
 * F x x x x x F
 * C x x x x x C
 * F x x x x x F
 * C F C F C F C
 *
 *
 * By now, this algorithm does not work 100%.
 *
 * @param graphS
 * @param graphST
 * @param PQ
 * @param nodes
 */
template<typename algebra_type>
void FAMGLevelCalculator<algebra_type>::rs_amg_coarsening_new(const cgraph &graphS, const cgraph &graphST, nodeinfo_pq_type &PQ, AMGNodes &nodes)
{
	try{
	AMG_PROFILE_FUNC();
	//return;
	std::vector<int> mark(nodes.size(), 0);
	std::vector<size_t> N1;
	std::vector<size_t> N2;

	//maxheap<int> myHeap(nodes.size(), mark);

	std::vector<size_t> markAsNotInterpolateable(nodes.size(), 0);
	PROGRESS_START(prog, nodes.size(), "FAMG::rs_amg_coarsening_new " << nodes.size());
	for(size_t i=0; i<nodes.size(); i++)
	{
		PROGRESS_UPDATE(prog, i);
		// reset all ratings.
		// 1. set node i coarse, neighbors fine
		// 2. coarsen (neighbors of fine node)
		// 4. calculate how many C-neighbors the F-neighbors j of i have
		N1.clear();
		N2.clear();

//		UG_LOG("\n\ntry setting [" << i << "] coarse\n");


		const int MARK_UNASSIGNED = 0;
		const int MARK_COARSE=-2;
		const int MARK_FINE=-10;

		mark[i] = MARK_COARSE;

		// mark neighbors fine
		for(cgraph::const_row_iterator conn = graphS.begin_row(i); conn != graphS.end_row(i); ++conn)
		{
			size_t j=*conn;
			mark[j] = MARK_FINE-1;
			N1.push_back(j);
//			UG_LOG("marked [" << j << "] fine\n");
		}

		// get N2
		for(cgraph::const_row_iterator conn = graphS.begin_row(i); conn != graphS.end_row(i); ++conn)
		{
			size_t j=*conn;
			for(cgraph::const_row_iterator conn2 = graphS.begin_row(j); conn2 != graphS.end_row(j); ++conn2)
			{
				size_t k=*conn2;
				if(mark[k] == MARK_UNASSIGNED)
				{
					mark[k] = 1;
					N2.push_back(k);
				}
				else if(mark[k] < MARK_UNASSIGNED) continue;
				else mark[k]++;
			}
		}
//		for(size_t j=0; j<N2.size(); j++) { UG_LOG("N2 " << N2[j] << " has mark " << mark[N2[j]] << "\n"); }

		if(N2.size() > 0)
		for(;;)
		{
			size_t kMax=N2[0];
			/// get max
			for(size_t j=1; j<N2.size(); j++)
			{
				size_t k=N2[j];
				if(mark[k] > mark[kMax])
					kMax=k;
			}

			if(mark[kMax] < MARK_UNASSIGNED) break;

//			UG_LOG(kMax << " gets coarse\n");

			// "set coarse"
			mark[kMax] = MARK_COARSE;
			// set neighbors fine
			for(cgraph::const_row_iterator conn = graphST.begin_row(kMax); conn != graphST.end_row(kMax); ++conn)
			{
				size_t k = *conn;
				if(mark[k] > MARK_UNASSIGNED)
				{
					mark[k] = MARK_FINE - 1;
				}
				else if(mark[k] <= MARK_FINE)
					mark[k]--;
			}
		}



		// 5. if j has 1 neighbor -> markedAsNotInterpolateable[j]++

		//    if j has 2 neighbors ->
		//       Linear: check if on one line . if note -> markedAsNotInterpolateable[j]++
		//       LE: check if LE-interp. if not -> markedAsNotInterpolateable[j]++

		stdvector<MathVector<3> > &pos = m_famg.m_amghelper.positions[level];

		for(size_t j=0; j<N1.size(); j++)
		{
			size_t k= N1[j];
			if(mark[k] == MARK_COARSE) continue;
			int Cneighbors = MARK_FINE - mark[k];

//			UG_LOG(k << " has " << Cneighbors << " coarse neighbors>\n");

			if(Cneighbors == 1)
			{
				markAsNotInterpolateable[k]++;
//				UG_LOG("marked [" << k << "] as not interpolateable (coarse node [" << i << "])\n");
			}
			else if(Cneighbors == 2)
			{
				markAsNotInterpolateable[k]++;
				continue;
				size_t i0=-1, i1;
				for(cgraph::const_row_iterator conn = graphS.begin_row(k); conn != graphS.end_row(k); ++conn)
					if(mark[*conn] == MARK_COARSE)
					{
						if(i0 == (size_t)-1)
							i0 = *conn;
						else { i1 = *conn; break; }

					}
				UG_LOG(k << " " << i0 << " " << i1 << "\n");


				double alpha[3];
				bool bSame[3];
				for(size_t b=0; b<3; b++)
				{
					if(fabs(pos[i1][b]-pos[i0][b])<1e-4)
						bSame[b]=true;
					else
						alpha[b] = (pos[k][b]-pos[i0][b]) /(pos[i1][b]-pos[i0][b]);
				}
				UG_LOG("\n");
				if((bSame[0] || bSame[1] || fabs(alpha[0] - alpha[1]) < 1e-4)
						&& (bSame[0] || bSame[1] || fabs(alpha[0] - alpha[2]) < 1e-4)
						&& (bSame[1] || bSame[2] || fabs(alpha[1] - alpha[2]) < 1e-4))
					continue;
				else
				{
					markAsNotInterpolateable[k]++;
					UG_LOG("marked [" << k << "] as not interpolateable (coarse node [" << i << "])\n");
				}
			}
			else continue;
		}

		for(size_t j=0; j<N2.size(); j++)
		{
			size_t k= N2[j];
			mark[k] = 0;
		}
		for(size_t j=0; j<N1.size(); j++)
		{
			size_t k= N1[j];
			mark[k] = 0;
		}
		mark[i] = 0;
	}

	SendVecFromTo(PN.comm(), PN.master_layout(), PN.slave_layout(), markAsNotInterpolateable);

	for(size_t i=0; i<markAsNotInterpolateable.size(); i++)
	{
		//if(markAsNotInterpolateable[i] > 0)
//			UG_LOG("node [" << i << "] markAsNotInterpolateable[i] = " << markAsNotInterpolateable[i] << "\n");
		if(nodes[i].rating < 0) continue;
		if(markAsNotInterpolateable[i] > 1)
		{
			SetCoarse(i, graphS, graphST, PQ, nodes, false, false);
			PQ.remove(i);
		}
		else if(markAsNotInterpolateable[i] > 0)
		{
			nodes[i].rating += 10;
			//nodes[i].update_time();
			PQ.update(i);
		}
	}
	}UG_CATCH_THROW_FUNC();
}
#endif


/**
 * do the external coarsening
 * steps are:
 * 1. get_external_coarsening_graph
 * 2. SetIsolatedCoarse
 * 3. graphST.set_as_transpose_of(graphS);
 * 4. RSAMG-CreateMeasureOfImportancePQ
 * 5. rs_amg_coarsening_LE_griebel / rs_amg_coarsening_new if enabled
 * 6. (parallel) Coarsen
 * 7. aggressive coarsening if enabeld
 * 8. RSAMG to FAMG coarsening
 */
template<typename algebra_type>
void FAMGLevelCalculator<algebra_type>::rs_amg_external_coarsening()
{
	try{
	AMG_PROFILE_FUNC();
	Stopwatch SW;
	if(m_bTiming) { UG_LOG("\nexternal coarsening... " );  SW.start(); }
	nodeinfo_pq_type PQ;
	size_t N = A_OL2.num_rows();

	// use RSAMG's methods to do standard rs coarsening
#ifdef UG_PARALLEL
	AMGNodes nodes(N, PN);
	UG_COND_THROW(m_famg.m_pParallelCoarsening == NULL, "please set_parallel_coarsening !");
	UG_DLOG(DBG_AMG_COARSENING, 0, m_famg.m_pParallelCoarsening->config_string() << "\n");
#else
	AMGNodes nodes(N);
#endif
	cgraph graphS, graphST;
	get_external_coarsening_graph(graphS);
	SetIsolatedCoarse(graphS, nodes);
	graphST.set_as_transpose_of(graphS);

	CreateMeasureOfImportancePQ(graphS, graphST, PQ, nodes);

	if(!m_famg.m_bBeam)
	{
#ifdef UG_PARALLEL
		if(m_famg.m_bLE_XC_Griebel)
			rs_amg_coarsening_LE_griebel(graphS, graphST, PQ, nodes);
		else if(m_famg.m_bXC_new)
			rs_amg_coarsening_new(graphS, graphST, PQ, nodes);

		m_famg.m_pParallelCoarsening->coarsen(PN, masterLayouts, slaveLayouts, graphS, graphST, PQ, nodes, true,
				m_famg.get_prevent_ff_connection_external() != 0);
#else
		Coarsen(graphS, graphST, PQ, nodes, true, false);

		if (m_famg.get_prevent_ff_connection_external())
		{
			PreventFFConnections(graphS, graphST, nodes);
		}
#endif



	}
	else
	{
		for(size_t i=0; i<nodes.size(); i++)
		{
			if(i == 0 || i == 1 || i == 256 || i == 257)
				//nodes.set_dirichlet(i);
				nodes.set_coarse(i);
			else if( (i/2) % 2 == 0 || i%2 == 0 || i < 100 || i > 130)
				nodes.set_coarse(i);
			else
				nodes.set_fine(i);
		}

	}

	// aggressive coarsening
	if(m_famg.m_bXCaggressive == true && level == 0)
	{
		AMG_PROFILE_BEGIN(rs_amg_external_coarsening_Agressive);
		UG_DLOG(DBG_AMG_COARSENING, 2, std::endl << "building graph2... ");
		UG_DLOG(DBG_AMG_COARSENING, 2, "unassigned = " << nodes.get_unassigned() << "\n");

		//unassigned = 0; ??
		cgraph graphAC(N);
		size_t m_iAggressiveCoarseningNrOfPaths = 2;
		CreateAggressiveCoarseningGraph(graphST, graphAC, nodes, m_iAggressiveCoarseningNrOfPaths);
		CreateMeasureOfImportanceAggressiveCoarseningPQ(graphAC, PQ, nodes);
		if(m_bTiming) { UG_LOG("took " << SW.ms() << " ms"); }
		// coarsen 2
		//------------------

		if(nodes.get_unassigned() == 0)
		{
			if(m_bTiming) { UG_LOG("\nskipping coarsening2: no unassigned nodes... " );  SW.start(); }
		}
		else
		{
			if(m_bTiming) { UG_LOG("\ncoarsening2... " );  SW.start(); }
			Coarsen(graphAC, graphAC, PQ, nodes, true, true);
			//PreventFFConnections(graphS, graphST, nodes);
		}
	}

	UG_LOG("XC coarsening rate: " << (nodes.get_nr_of_coarse()*100.0)/nodes.size());
	// coarsening done. transfer AMG nodeinfo -> FAMG nodeinfo

	AMGnodeinfoToFAMGnodeinfo(nodes, rating, graphS);

	if(m_bTiming) { UG_LOG("took " << SW.ms() << " ms"); }


	//----------------------
	if(m_famg.m_writeMatrices)
		write_XC_coarsening(graphS, nodes, "");
	//if(nodes[i].is_unassigned_fine_indirect())
	//	rating.set_aggressive_fine(i);
	}UG_CATCH_THROW_FUNC();
}

/**
 * write the coarsening out to
 * matrixWritePath/level/XC+addStr.mat
 * @param graphS 	the strong connected graph
 * @param nodes		RS AMGNodes defining fine/coarse
 * @param addStr	additional String (like _post)
 */
template<typename algebra_type>
void FAMGLevelCalculator<algebra_type>::write_XC_coarsening(cgraph &graphS, AMGNodes &nodes, std::string addStr)
{
	PROFILE_FUNC_GROUP("algebra AMG debug");
	if(!m_famg.m_writeMatrices) return;
	size_t N = nodes.size();
	PROGRESS_START(prog, N, "write_XC_coarsening " << N);

	std::string name = std::string("XC") + addStr;
	std::string filename = m_famg.GetAMGFilenamePath(A_OL2, level, name, level, ".mat", true);
	std::string marksname = m_famg.GetAMGFilenamePath(A_OL2, level, name, level, ".marks", true);
	std::fstream file(filename.c_str(), std::ios::out);
	ConnectionViewer::WriteGridHeader(file, m_famg.m_amghelper.positions[level], m_famg.m_amghelper.dimension);
	for(size_t i=0; i < graphS.size(); i++)
	{
		PROGRESS_UPDATE(prog, i);
		for(cgraph::row_iterator it = graphS.begin_row(i); it != graphS.end_row(i); ++it)
			file << i << " " << (*it) << " 1" << std::endl;
	}
	file << "c " << marksname << "\n";



	std::fstream fcoarse(marksname.c_str(), std::ios::out);
	fcoarse << "0 0 1 1 2\n";

	for(size_t i=0; i<N; i++)
	{
		if(nodes[i].is_coarse())
			fcoarse << i << "\n";
	}
}





template<typename algebra_type>
void FAMGLevelCalculator<algebra_type>::write_coarsening()
{
	/*PROFILE_FUNC_GROUP("debug");

	stdvector<MathVector<3> > &vec = m_famg.m_amghelper.positions[level];

	std::ios_base::openmode mode;
#ifdef UG_PARALLEL
	if(pcl::ProcRank()!=0)
		mode = std::ios::out | std::ios::app | std::ios::binary;
	else
#endif
	{

		CreateDirectory((std::string(m_famg.m_writeMatrixPath) + "/coarsening").c_str(), 0777);
		mode = std::ios::out | std::ios::binary;
	}


#ifdef UG_PARALLEL
	for(int i=0; i<pcl::ProcRank(); i++)
		pcl::AllProcsTrue(true);
	{
#endif
		std::fstream file((std::string(m_famg.m_writeMatrixPath) + "/coarsening/rating" + ToString(level)).c_str(), mode);
		for(size_t i=0; i<vec.size(); i++)
		{
			Serialize(file, vec[i]);
			Serialize(file, rating[i].rating);
		}
#ifdef UG_PARALLEL
	}
	for(int i=pcl::ProcRank(); i<pcl::NumProcs(); i++)
		pcl::AllProcsTrue(true);
#endif
*/
}

template<int Ti>
struct MathVectorComp {
  inline bool operator() (const MathVector<Ti> &a, const MathVector<Ti> &b) const
  {
	  for(size_t i=0; i<Ti; i++)
		  if(a[i] != b[i]) return a[i] < b[i];
	  return false;
  }
};

template<typename algebra_type>
void FAMGLevelCalculator<algebra_type>::read_coarsening()
{
	/*
	 PROFILE_FUNC_GROUP("debug");

	stdvector<MathVector<3> > &vec = m_famg.m_amghelper.positions[level];
	typedef std::map<MathVector<3>, size_t, MathVectorComp<3> > MathVecMap;
	MathVecMap m;
	for(size_t i=0; i<vec.size(); i++)
		m[vec[i]] = i;

	std::fstream file((std::string(m_famg.m_writeMatrixPath) + "/coarsening/rating" + ToString(level)).c_str(), std::ios::in | std::ios::binary);
	for(size_t i=0; i<vec.size(); i++)
	{
		MathVector<3> pos;
		double r;
		Deserialize(file, pos);
		Deserialize(file, r);

		MathVecMap::iterator it = m.find(pos);
		if(it != m.end())
		{
			rating[it->second].rating = r;
			UG_LOG("Setting node " << it->second << " to " << rating[it->second]);
		}
	}
	*/
}

// end group plugins_amg_famg
/// \}

}

#endif /* FAMG_EXTERNAL_COARSENING_H_ */
