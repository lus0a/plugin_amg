/*
 * Copyright (c) 2010-2014:  G-CSC, Goethe University Frankfurt
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

///
#include <iostream>

#include "common/assert.h"
#include "common/log.h"
#include "common/progress.h"

using namespace std;

#include "lib_algebra/common/graph/graph.h"
#include "rsamg_nodeinfo.h"

#include "lib_algebra/common/stl_debug.h"
#include "../amg_profiling.h"
//#define AMG_PRINT_COARSEN

#include "lib_algebra/common/heap/boxsort.h"
#include "lib_algebra/algebra_common/unsorted_sparse_vector.h"
#include "common/debug_print.h"

namespace ug
{

extern DebugID DBG_AMG_COARSENING; // SetDebug(debugId.RSAMG.Coarsening, 3)

/// \addtogroup plugins_amg_rsamg
/// \{

time_t global_t=0;


// CreateMeasureOfImportancePQ:
//------------------------------
/**
 * \brief Creates the "measure of importance" (=rating) priority queue.
 * calculate ratings in nodes[i].rating, build up priority queue PQ.
 * rating = unassigned neighbors + 2 * fine neighbors ( neighbors in strongT)
 * \param 	strong		strong[i] = nodes which influence i
 * \param 	strongT 	strongT[i] = nodes which are influenced by i.
 * \param	PQ			maxheap priority queue for sorting of the nodes wrt the rating
 * \param	nodes
 * \sa	CreateStrongConnectionGraph
 */
void CreateMeasureOfImportancePQ(const cgraph &strong, const cgraph &strongT, nodeinfo_pq_type &PQ, AMGNodes &nodes)
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
		else if(!strong.is_isolated(i))
		{
			//UG_ASSERT(graph.iNrOfConnections[i] > 0, "node " << i << " has " << graph.iNrOfConnections[i] << " connections?");
			nodes.set_rating(i, strongT.num_connections(i));
			PQ.insert_item(i);
		}
	}
}


///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// CreateAggressiveCoarseningGraph
//----------------------------------
/**
 * Creates the graph for aggressive coarsening from the strong connectivity graph
 * only among nodes which are marked coarse.
 * two coarse nodes a and b in graph2 are connected, if
 * - there exist at least 2 ways of length 2 from a to b (A2-Coarsening)
 * - there exist a ways of length 2 from a to b (A1-Coarsening)
 * because of the coarsening process, every coarse node has only fine neighbors in graph,
 * that means those ways are from a coarse node over a fine node to a coarse node.
 *
 * \param	graph				old graph of normal strong connectivity (from CreateGraph)
 * \param 	graph2				new graph of "distance-2-strong-connectivity"
 * \param	nodes				info about nodes if nodes are fine or coarse
 * \param	minNrOfPaths		2 or 1, depending on A2- or A1-Coarsening.
 */
void CreateAggressiveCoarseningGraph(const cgraph &graph, cgraph &graph2, const AMGNodes &nodes,
		int minNrOfPaths)
{
	try	{
	AMG_PROFILE_FUNC();
	UG_ASSERT(nodes.size() == graph.size(), "");
	vector<int> connection(255);
	UnsortedSparseVector<int> nrOfPaths(graph.size());

	//graph.print();
	for(size_t i=0; i < graph.size(); i++)
	{
		if(nodes[i].is_fine_direct())
			continue;

		nrOfPaths.clear();
		// first calculate all nodes reachable with paths of length 2

		// ! i is coarse -> has only fine neighbors
		for(cgraph::const_row_iterator conn = graph.begin_row(i); conn != graph.end_row(i); ++conn)
		{
			size_t indexN = (*conn);
			for(cgraph::const_row_iterator connN = graph.begin_row(indexN); connN != graph.end_row(indexN); ++connN)
			{
				size_t indexNN = (*connN);

				if(indexNN == i || nodes[indexNN].is_fine_direct())
					continue;
				nrOfPaths(indexNN)++;
			}
		}

		// then sort out those which were reached #aggressiveCoarseningNrOfPaths (2 or 1) times
		for(UnsortedSparseVector<int>::iterator it = nrOfPaths.begin(); it != nrOfPaths.end(); ++it)
		{
			if(it->value() >= minNrOfPaths)
			{
				// add connection i -> node
				graph2.set_connection(i, it->index());
				// increase rating of i
			}
		}
	}
	}UG_CATCH_THROW_FUNC()
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// CreateMeasureOfImportanceAggressiveCoarseningPQ
//----------------------------------
/**
 * this function could be merged with CreateMeasureOfImportancePQ
 * \param	graphAC			graph created by CreateAggressiveCoarseningGraph
 * \param 	PQ				maxheap priority queue for sorting of the nodes wrt the rating
 * \param	nodes			fine/coarse marks (skip fine nodes)
 */
void CreateMeasureOfImportanceAggressiveCoarseningPQ(const cgraph &graphAC, nodeinfo_pq_type &PQ, AMGNodes &nodes)
{
	try{
	AMG_PROFILE_FUNC();
	UG_ASSERT(nodes.size() == graphAC.size(), "");
	PQ.create(nodes.size(), &nodes[0]);

	for(size_t i=0; i<nodes.size(); i++)
	{
		if(nodes.is_slave(i)) continue;
		if(nodes[i].is_fine_direct() || nodes[i].is_dirichlet())
			continue;

		nodes.set_rating(i, graphAC.num_connections(i));
		PQ.insert_item(i);
	}

	}UG_CATCH_THROW_FUNC();
}


///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// RemoveUnassignedNeighbors:
//-------------------------
void RemoveUnassignedNeighbors(const cgraph &graph, nodeinfo_pq_type &PQ, AMGNodes &nodes, size_t i)
{
	for(cgraph::const_row_iterator conn = graph.begin_row(i); conn != graph.end_row(i); ++conn)
	{
		size_t indexN = (*conn);
		if(nodes[indexN].is_assigned()) continue;
		PQ.remove(indexN);
	}
}

void ChangeRatingOfUnassignedNeighbors(const cgraph &graph, nodeinfo_pq_type &PQ, AMGNodes &nodes, size_t i, int change)
{
	for(cgraph::const_row_iterator conn = graph.begin_row(i); conn != graph.end_row(i); ++conn)
	{
		int indexN = (*conn);
		// TODO: perhaps we could create a the f-f candidate list here

		if(nodes[indexN].is_assigned())
			continue;
		nodes[indexN].rating += change;
		nodes[indexN].update_time();
		PQ.update(indexN);
	}

}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// MarkUnassignedNeighborsFine:
//-------------------------
void MarkUnassignedNeighborsFine(const cgraph &S, const cgraph &ST, nodeinfo_pq_type &PQ, AMGNodes &nodes, size_t i, bool bMarkAsFineIndirect)
{
	for(cgraph::const_row_iterator conn = ST.begin_row(i); conn != ST.end_row(i); ++conn)
	{
		int indexN = (*conn);

		IF_DEBUG(DBG_AMG_COARSENING, 3)
		{
			UG_LOG(indexN << " ");
			if(nodes[indexN].is_assigned())
				UG_LOG((nodes[indexN].is_coarse() ? "=c " : "=f "));
		}

		if(nodes[indexN].is_parallel_dont_care() == false && nodes[indexN].is_assigned()) continue;

		//if(bIndirect) nodes[indexN].setFineIndirect();
		//else

		if(bMarkAsFineIndirect)
			nodes.set_unassigned_fine_indirect(indexN);
		else
			nodes.set_fine_direct(indexN);

		IF_DEBUG(DBG_AMG_COARSENING, 3)
		{
			UG_LOG("->f ");
		}

		ChangeRatingOfUnassignedNeighbors(S, PQ, nodes, indexN, +1);
	}
}

void SetCoarse(int i, const cgraph &S, const cgraph &ST, nodeinfo_pq_type &PQ, AMGNodes &nodes, bool bUnsymmetric, bool bAggressiveCoarsening)
{
	UG_DLOG(DBG_AMG_COARSENING, 3, endl << "set coarse: " << i << " rating " << nodes[i].rating  << ". then fine: ");
	UG_ASSERT(!nodes[i].is_assigned(), "node " << i << " is already assigned??? (rating = " << nodes[i].rating << ", unassigend = " << nodes.get_unassigned() << ")");

	nodes.set_coarse(i);					// mark as coarse/assigned

	// unassigned neighbors will be marked as fine, so remove from PQ
	RemoveUnassignedNeighbors(ST, PQ, nodes, i);

	// now mark unassigned neighbors as fine
	MarkUnassignedNeighborsFine(S, ST, PQ, nodes, i, bAggressiveCoarsening);

	if(bUnsymmetric)
		ChangeRatingOfUnassignedNeighbors(S, PQ, nodes, i, -1);
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Coarsen:
//-------------------------
/**
 * \brief Coarsens the graph with ratings of nodes in nodes[i].rating, set up in a priority queue PQ
 * \param S							S[i] nodes influencing i
 * \param ST						ST[i] = nodes beeing influenced by i
 * \param PQ
 * \param nodes
 * \param bUnsymmetric				set to true if S != ST
 * \param bAggressiveCoarsening
 * \return number of new coarse nodes.
 */
int Coarsen(const cgraph &S, const cgraph &ST, nodeinfo_pq_type &PQ, AMGNodes &nodes, bool bUnsymmetric, bool bAggressiveCoarsening)
{
	try{
	AMG_PROFILE_FUNC();
	PROGRESS_START(prog, nodes.size(), "Coarsening graph size " << nodes.size());
	THROW_IF_NOT_EQUAL_3(S.size(), nodes.size(), PQ.arr_size());

#ifdef UG_DEBUG
	size_t n=0;
//	nodes.print();
	for(size_t i=0; i<nodes.size(); i++)
	{
		if(nodes[i].is_assigned() == false)
		{	UG_COND_THROW(!PQ.is_in(i), i << " ? " << nodes[i]); }
		if(PQ.is_in(i))
		{ 	UG_COND_THROW(nodes[i].is_assigned(), i << " " << nodes[i]); }
		if(nodes[i].is_assigned() == false)
			n++;
	}
	UG_COND_THROW(n != nodes.get_unassigned(), "");
#endif

	while(PQ.height() > 0)
	{
		PROGRESS_UPDATE(prog, PQ.height());
		UG_ASSERT(PQ.height() == nodes.get_unassigned(), PQ.height() << " != " << nodes.get_unassigned())
		// get Node with best rating
		int best = PQ.remove_max();

		SetCoarse(best, S, ST, PQ, nodes, bUnsymmetric, bAggressiveCoarsening);
	}
	//cout << endl;

//	UG_ASSERT(nodes.get_nr_of_coarse() > 0, "no coarse nodes???");
	UG_DLOG(DBG_AMG_COARSENING, 1, nodes.get_nr_of_coarse() << " Coarse Nodes. ");

	IF_DEBUG(DBG_AMG_COARSENING, 3)
	{
		UG_LINE_LOG();
		nodes.print();
	}

	return nodes.get_nr_of_coarse();

	}UG_CATCH_THROW_FUNC();
	return 0;
}


/**
 * \brief
 * \param graphS
 * \param graphST
 * \param nodes
 * \param i
 */
void PreventFFConnection(const cgraph &graphS, const cgraph &graphST, AMGNodes &nodes, size_t i, vector<bool> &marks, int &nrOfFFCoarseNodes)
{
	try{
	if(nodes[i].is_coarse() || graphS.num_connections(i)==0) return;

	// mark coarse nodes interpolating this fine node
	for(cgraph::const_row_iterator it = graphS.begin_row(i); it != graphS.end_row(i); ++it)
	{
		if(nodes[(*it)].is_coarse())
			marks[(*it)] = true;
	}

	// prevent strong F-F connections without common Interpolation node
	for(cgraph::const_row_iterator it = graphS.begin_row(i); it != graphS.end_row(i); ++it)
	{
		if(nodes[*it].is_coarse() || graphS.is_isolated(*it))
			continue;

		cgraph::const_row_iterator it2 = graphS.begin_row(*it);
		cgraph::const_row_iterator it2end = graphS.end_row(*it);
		for(; it2 != it2end; ++it2)
		{
			if(nodes[*it2].is_coarse() && marks[*it2])
				break;
		}
		if(it2 == it2end)
		{
			nodes.set_coarse(i);
			nrOfFFCoarseNodes++;
			break;
		}
	}

	// remove marks
	for(cgraph::const_row_iterator conn = graphS.begin_row(i); conn != graphS.end_row(i); ++conn)
	{
		size_t index = (*conn);
		if(nodes[index].is_coarse())
			marks[index] = false;
	}

	}UG_CATCH_THROW_FUNC();
}


///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// PreventFFConnections:
//-------------------------
/**
 * \brief When two fine nodes are strong connected, they need to have a common interpolating
 * node. Otherwise, we set
 * \param graphS
 * \param graphST
 * \param nodes
 * \return				returns number of new coarse nodes.
 */
void PreventFFConnections(const cgraph &graphS, const cgraph &graphST, AMGNodes &nodes)
{
	try{
	AMG_PROFILE_FUNC();
	size_t N = graphS.size();
	UG_ASSERT(N == graphST.size() && N == nodes.size(), "");

	int nrOfFFCoarseNodes=0;
	vector<bool> marks(N, false);

	for(size_t i=0; i< N; i++)
	{
		if(nodes.is_inner(i) && nodes[i].is_fine())
			PreventFFConnection(graphS, graphST, nodes, i, marks, nrOfFFCoarseNodes);
	}

	if(nrOfFFCoarseNodes)
		UG_DLOG(DBG_AMG_COARSENING, 0, "F-F prevention, now " << nodes.get_nr_of_coarse() << " coarse nodes.\n");

	}UG_CATCH_THROW_FUNC();
}

/*

template<typename matrix_type>
double GetRowMaxAbs(const matrix_type &A, size_t row)
{
	double m=0;
	for(typename matrix_type::const_row_iterator it = A.begin_row(row); it != A.end_row(row); ++it)
	{
		double a = BlockAbs(*it);
		if(a > m) m = a;
	}
	return m;
}

void PreventFFConnectionsMetsch(const matrix_type &mat, const cgraph &graphS, const cgraph &graphST, AMGNodes &nodes)
{
	AMG_PROFILE_FUNC();
	size_t N = graphS.size();
	UG_ASSERT(N == graphST.size() && N == nodes.size(), "");

	int nrOfFFCoarseNodes=0;
	vector<bool> marks(N, false);
	vector<double> maxk(N, -1.0);

	for(size_t r=0; i< N; r++)
		maxk[i] = GetRowMaxAbs(mat);


	for(size_t r=0; i< N; r++)
	{
		if(nodes.is_fine(r)==false)
			continue;

		for(cgraph::const_row_iterator it = graphS.begin_row(r); it != graphS.end_row(r); ++it)
		{
			int c=*it;
			if(!nodes[c].is_fine()) continue;



		}

	}

	if(nrOfFFCoarseNodes)
		UG_DLOG(DBG_AMG_COARSENING, 1, "F-F prevention, now " << nodes.get_nr_of_coarse() << " coarse nodes.\n");
}


*/

/**
 * \brief
 * \param graphS
 * \param graphST
 * \param nodes
 * \param i
 * \param G
 * \param ratings
 */
void GetFFConnectionGraph(const cgraph &graphS, const cgraph &graphST, AMGNodes &nodes, size_t i, vector<bool> &marks,
		cgraph &G, stdvector<size_t> &ratings)
{
	try{
	if(nodes[i].is_coarse() || graphS.num_connections(i)==0) return;

	// mark coarse nodes interpolating this fine node
	for(cgraph::const_row_iterator it = graphST.begin_row(i); it != graphST.end_row(i); ++it)
	{
		if(nodes[(*it)].is_coarse())
			marks[(*it)] = true;
	}

	// prevent strong F-F connections without common Interpolation node
	for(cgraph::const_row_iterator it = graphST.begin_row(i); it != graphST.end_row(i); ++it)
	{
		if(nodes.needs_assignment(*it)==false ) continue;
		if(nodes[*it].is_coarse() || graphS.is_isolated(*it))
			continue;

		cgraph::const_row_iterator it2 = graphST.begin_row(*it);
		cgraph::const_row_iterator it2end = graphST.end_row(*it);
		for(; it2 != it2end; ++it2)
		{
			if(nodes[*it2].is_coarse() && marks[*it2])
				break;
		}
		if(it2 == it2end)
		{
			G.set_connection(i, *it);
			G.set_connection(*it, i);
			ratings[i]++;
			ratings[*it]++;
			break;
		}
	}

	// remove marks
	for(cgraph::const_row_iterator conn = graphST.begin_row(i); conn != graphST.end_row(i); ++conn)
	{
		size_t index = (*conn);
		if(nodes[index].is_coarse())
			marks[index] = false;
	}

	}UG_CATCH_THROW_FUNC();
}

/**
 * \brief
 * \param graphS
 * \param graphST
 * \param nodes
 * \return				returns number of new coarse nodes.
 */
void PreventFFConnectionsMinimal(const cgraph &graphS, const cgraph &graphST, AMGNodes &nodes)
{
	try{
	AMG_PROFILE_FUNC();
	size_t N = graphS.size();
	UG_ASSERT(N == graphST.size() && N == nodes.size(), "");

	int nrOfFFCoarseNodes=0;
	vector<bool> marks(N, false);

	cgraph G(graphST.size());
	stdvector<size_t> ratings(graphST.size(), 0);


	for(size_t i=0; i< N; i++)
	{
		if(nodes.is_inner(i)==false)
			continue;
		GetFFConnectionGraph(graphS, graphST, nodes, i, marks, G, ratings);
	}

	BoxPriorityQueue<size_t> pq(ratings);

	for(size_t i=0; i<ratings.size(); i++)
		if(ratings[i]>0)
			pq.insert_item(i);

	while(pq.height() > 0)
	{
		int best = pq.remove_max();
		nodes.set_coarse(best);
		nrOfFFCoarseNodes++;
		for(cgraph::row_iterator it = G.begin_row(best); it != G.end_row(best); ++it)
		{
			ratings[*it]--;
			if(ratings[*it]==0)	pq.remove(*it);
			else				{ pq.update(*it); }
		}
	}
	if(nrOfFFCoarseNodes)
		UG_DLOG(DBG_AMG_COARSENING, 1, "F-F prevention, now " << nodes.get_nr_of_coarse() << " coarse nodes.\n");

	}UG_CATCH_THROW_FUNC();
}


// end group plugins_amg_rsamg
/// \}

} // namespace ug

