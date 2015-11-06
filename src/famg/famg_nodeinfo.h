/*
 * Copyright (c) 2012-2014:  G-CSC, Goethe University Frankfurt
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

#ifndef __H__UG__LIB_DISC__FAMG_SOLVER__FAMG_NODEINFO_H__
#define __H__UG__LIB_DISC__FAMG_SOLVER__FAMG_NODEINFO_H__

namespace ug {

/// \addtogroup plugins_amg_famg
/// \{

//  structs
// fine is the smallest (see is_fine)
#define FAMG_BIGGEST_FINE_RATING		(-30)
#define FAMG_AGGRESSIVE_FINE_RATING		(-31)
#define FAMG_UNCALCULATED_FINE_RATING	(-32)
#define FAMG_FINE_RATING				(-33)
//      FAMG_INDIRECT_LEVEL[i]          (-33-i)


#define FAMG_UNINTERPOLATEABLE			(-10)
#define FAMG_COARSE_RATING				(-11)
#define FAMG_DIRICHLET_RATING			(-12)
// those are border-nodes which are fine on another processor. we need to calculate ratings for them


//#define FAMG_COARSENING_USE_F_FOR_SORT

/**
 * this class contains the rating and the state of a node
 * rating is set to  min_{possible parents} (number of nodes which have to set coarse additionally)  + (F-value of this pair) / delta
 * depending on FAMG_COARSENING_USE_F_FOR_SORT, the additional F-value is ignore or not
 */
class FAMGNode
{
	template<typename> friend class FAMGNodes;
public:


private:
	/// rating is set to  min_{possible parents} (number of nodes which have to set coarse additionally)  + (F-value of this pair) / delta
	/// depending on FAMG_COARSENING_USE_F_FOR_SORT, the additional F-value is ignore or not
	double rating;

private:
	/// set this node to fine
	inline void set_fine(){rating = FAMG_FINE_RATING;}

	/// set this node to fine indirect level
	inline void set_fine_indirect_level(int level){rating = FAMG_FINE_RATING-level;}

	/// set this node to aggressive fine (a pre-step before fine indirect)
	inline void set_aggressive_fine(){rating = FAMG_AGGRESSIVE_FINE_RATING;}

	/// set this node to coarse
	inline void set_coarse(){rating = FAMG_COARSE_RATING;}

	/// set this node to uninterpolateable
	inline void set_uninterpolateable()	{ rating = FAMG_UNINTERPOLATEABLE;	}

	/// set this node to dirichlet (dirichlet nodes are nodes with A_{ij} = { 1 iff i=j else 0)
	inline void set_dirichlet() { rating = FAMG_DIRICHLET_RATING; }

	// uncalculated: only used for on_demand_coarsening
	inline void set_uncalculated_fine() { rating = FAMG_UNCALCULATED_FINE_RATING; }

public:
	FAMGNode() { rating = 0.0; }
	double get_rating() const { return rating; }


	// is_fine is true if node is fine OR uncalculated_fine OR aggressive_fine

	/// \return true if node is fine (fine, aggressive fine or uncalculated fine)
	inline bool is_fine() const { return rating <= FAMG_BIGGEST_FINE_RATING; }

	/// \return true if node is 'direct' fine, that is a node which is fine and interpolates only from coarse
	inline bool is_fine_direct() const { return rating == FAMG_FINE_RATING; }

	/// \return true if aggressive fine, that is a node which should be interpolated by fine nodes.
	inline bool is_aggressive_fine() const { return rating == FAMG_AGGRESSIVE_FINE_RATING; }


	/// \return true if node is coarse
	inline bool is_coarse() const { return rating == FAMG_COARSE_RATING; }

	/// \return true if node is uninterpolateable (-> coarse)
	inline bool is_uninterpolateable() const { return rating == FAMG_UNINTERPOLATEABLE; }

	/// \return true if nodes are 'dirichlet' nodes, that means A_{ij} = { 1 iff i=j else 0
	inline bool is_dirichlet() const { return rating == FAMG_DIRICHLET_RATING; }

	/// uncalculated: only used for on_demand_coarsening
	inline bool is_uncalculated_fine() const { return rating == FAMG_UNCALCULATED_FINE_RATING; }

	inline bool is_interpolateable_from() const { return is_dirichlet() || is_coarse(); }


public:
	/// \true if we could set this node coarse
	inline bool could_be_coarse() const
	{
		return is_valid_rating() || is_coarse() || is_uninterpolateable(); // is_dirichlet ??
	}

	/// \return true if node has a lambda rating which means it is not assigned
	inline bool is_valid_rating() const
	{
		return rating >= 0;
	}
	/// \return true if node is assigned (assigned coarse/fine or similar)
	inline bool is_assigned()
	{
		return rating < 0;
	}


	/// used in the heap operation to construct the coarse grid
	inline bool operator > (const FAMGNode &other) const
	{
#ifdef FAMG_COARSENING_USE_F_FOR_SORT
		if(rating == other.rating)
			return this < &other;
		else
			return rating < other.rating;
#else
		return floor(rating) < floor(other.rating);
#endif
	}

	// this returns the state as a char for use in parallelization algorithms
	inline char get_state() const
	{
		UG_ASSERT(rating <= 0 && rating == floor(rating), "rating is " << rating << ", but needs to be <= 0 and integer");
		return (char) rating;
	}

	inline double get_val() const
	{
		UG_ASSERT(rating >= 0 && rating < 1000, "rating is " << rating << ", out of bounds [0, 1000]");
		return rating;
	}

	friend std::ostream &operator << (std::ostream &out, const FAMGNode &n)
	{
		out << "Rating: " << n.rating;
		if(n.rating == FAMG_FINE_RATING) out << " (fine)";
		else if(n.rating == FAMG_AGGRESSIVE_FINE_RATING) out << " (AC fine)";
		else if(n.rating == FAMG_UNCALCULATED_FINE_RATING) out << " (fine u)";
		else if(n.rating == FAMG_UNINTERPOLATEABLE) out << " (uninterpolateable)";
		else if(n.is_coarse()) out << " (coarse)";
		else if(n.rating == FAMG_DIRICHLET_RATING) out << " (dirichlet)";
		else if(n.rating == 0.0) out << "(unassigned)";
		return out;
	}

	bool operator != (const FAMGNode &other) const
	{
		return rating != other.rating;
	}
	bool operator == (const FAMGNode &other) const
	{
		return rating == other.rating;
	}
};

//! this structure holds coarse/fine information as well as master/slave
//! for use in amg. it tries to replace the missing "node" object here
template<typename TProlongationMatrix>
class FAMGNodes
{
public:
	typedef FAMGNode value_type;

#ifdef UG_PARALLEL
	FAMGNodes(TProlongationMatrix &_P, size_t level, ug::cAMG_helper &amghelper,
			ParallelNodes &pn)
	: m_level(level), P(_P), m_amghelper(amghelper), PN(pn)
	{
		m_iNrOfCoarse = 0;
	}
#else
	FAMGNodes(TProlongationMatrix &_P, size_t level, ug::cAMG_helper &amghelper)
		: m_level(level), P(_P), m_amghelper(amghelper)
		{
			m_iNrOfCoarse = 0;
		}
#endif


	void create(size_t size)
	{
		m_iNrOfCoarse = 0;
		m_iNrOfFine = 0;
		m_iNrOfAggressiveFine = 0;
		nodes.clear(); 		nodes.resize(size);
#ifdef UG_PARALLEL
		OLtype.clear();		OLtype.resize(size, 0);
#endif
	}

	void resize(size_t size)
	{
		size_t oldSize = nodes.size();
		nodes.resize(size);
		(void) oldSize;
#ifdef UG_PARALLEL
		OLtype.resize(size);
		for(size_t i=oldSize; i<size; i++)
		{

		}
#endif
	}

	FAMGNode &operator [] (size_t i) { return nodes[i]; }
	const FAMGNode &operator [] (size_t i) const { return nodes[i]; }

	/// calculates the rating of this nodes for a list of parent nodes
	template<typename parent_nodes_vec_type>
	int get_rating(const parent_nodes_vec_type &M)
	{
		//AMG_PROFILE_FUNC();
		size_t rating = 0;
		for(size_t i=0; i<M.size(); ++i)
		{
			const FAMGNode &ninfo = nodes[M[i].from];
			if(ninfo.is_fine())
			{
				UG_DLOG(LIB_ALG_AMG, 4, " pair [" << M[0].from << "], [" << M[1].from << "] is invalid, since " << get_original_index(M[i].from) << " is fine. ");
				return -1;
			}
			/*else if(ninfo.is_uninterpolateable())
				rating += 1000;
			else*/
			if(ninfo.is_coarse() || ninfo.is_uninterpolateable()) continue;
			else rating++;
		}
		return rating;
	}

	/**
	 * updates the rating of node <node>
	 * and writes it to rating
	 * possible neighbors are always stored as follows:
	 * - possible_neighbors[i][0] is the best currently known pair for interpolating i
	 * - this can change
	 * a) if a parent node gets fine -> removal of this pair
	 * b) if a parent node gets coarse
	 * @param node
	 * @param PN parent nodes
	 * @return true if value has been updated, otherwise false
	 */
	template<typename neighborstruct>
	bool update_rating(size_t node, stdvector<neighborstruct> &PN, bool onlyRecalcBest=false)
	{
		if(onlyRecalcBest ==false && i_must_assign(node) == false)
		{
			return false;
		}

		//AMG_PROFILE_FUNC();
		UG_DLOG(LIB_ALG_AMG, 4, " update rating of node [" << node << "]... ");
		if(onlyRecalcBest ==false && nodes[node].is_valid_rating() == false)
		{
			UG_DLOG(LIB_ALG_AMG, 4, nodes[node].rating << " is not a valid rating."); return false;
		}

		int mini = -1;
		double minrating = 10000;
		//double minF = 1e12;
		UG_DLOG(LIB_ALG_AMG, 4, "has " << PN.size() << " possible parent nodes. ");
		for(size_t i=0; i<PN.size(); )
		{
			double irating = get_rating(PN[i].parents);

			if(irating < 0)
			{
				std::swap(PN[i], PN.back()); // remove this pair
				PN.resize(PN.size()-1);
				UG_DLOG(LIB_ALG_AMG, 4, " removed pair " << i << " ");
				continue;
			}

			UG_DLOG(LIB_ALG_AMG, 4,PN[i]<<"\nrating=" << irating << "\n");
			irating += PN[i].F/1000;

			// choose the pairs with best minimization F.
			if(irating < minrating) // || (irating == minrating && PN[i].F < minF))
			{
				minrating = irating;
				//minF = PN[i].F;
				mini = i;
			}
			i++;

		}
		UG_DLOG(LIB_ALG_AMG, 4, "now has " << PN.size() << " possible parent nodes. ");
		for(size_t i=0; i<PN.size(); i++)
			UG_DLOG(LIB_ALG_AMG, 4, "[" << PN[i].parents[0].from << "]-[" << PN[i].parents[1].from << "] ");
		UG_DLOG(LIB_ALG_AMG, 4, "\nmini=" << mini << "\n");

		if(mini != -1 && mini != 0) std::swap(PN[0], PN[mini]);
		if(onlyRecalcBest) return true;

		if(mini != -1)
		{
			UG_DLOG(LIB_ALG_AMG, 4, " has rating " << minrating << "\n");

			if(nodes[node].rating != minrating)
			{
				UG_DLOG(LIB_ALG_AMG, 4, " new rating! ");
				nodes[node].rating = minrating;
				return true;
			}
			else return false;
		}
		else
		{
			UG_DLOG(LIB_ALG_AMG, 4, " is uninterpolateable! ");
			if(nodes[node].is_uninterpolateable())
				return false;
			else
			{
				set_uninterpolateable(node);
				return true;
			}
		}
	}


	size_t size() const { return nodes.size(); }



#ifdef UG_PARALLEL
public:
	/// \return true if node is a 'inner' node, i.e. nodes not lying on the boundary
	bool is_inner_node(size_t i) const
	{
		return PN.is_inner(i);
	}

	/// \return true if node is on boundary and master
	bool is_master(size_t i) const
	{
		return PN.is_master(i);
	}

	/// \return true if node is on boundary and slave
	bool is_slave(size_t i) const
	{
		return PN.is_slave(i);
	}

	/// \return 0 for master/inner nodes, else the distance
	/// to master/inner nodes. this is != 0 for overlap nodes
	size_t distance_to_master_or_inner(size_t i) const
	{
		return PN.distance_to_master_or_inner(i);
	}

	/*bool is_master_on(size_t i, int pid)
	{
		if(is_slave(i, 0))
			return m_masterOn[i] == pid;
		else
			return false;
	}*/

	/// \return true if we must assign this node:
	/// this is true if it is an inner node or it is master
	bool i_must_assign(size_t i) const
	{
		return is_inner_node(i) || is_master(i);
	}

	/// \return true if this process can set this node coarse on its on
	/// this is true if we can set it coarse (i.e. not already fine)
	/// and it is master/inner or has distance one to a master or inner
	bool i_can_set_coarse(size_t i) const
	{
		return nodes[i].could_be_coarse() && distance_to_master_or_inner(i) <= 1;
	}

	/// returns the overlap type (inner, master, slave or distanceToMasterOrInner=X)
	std::string OL_type(size_t i) const
	{
		std::stringstream  ss;
		ss << PN.overlap_type(i);
		return ss.str();
	}

	/// returns some info about the node i
	std::string info(size_t i) const
	{
		return Stringify() << "Index [" << i << "]: " << nodes[i] << ", " << PN.overlap_type(i);
	}

	/// calculate the number of unassigned nodes
	void calculate_unassigned()
	{
		m_iUnassigned = 0;
		for(size_t i=0; i<size(); i++)
		{
			if(i_must_assign(i))
				m_iUnassigned++;
		}
	}

#else
public:
	std::string info(size_t i) const
	{
		std::stringstream ss;
		ss << "Index " << i << ": " << nodes[i];
		return ss.str();
	}

	bool is_inner_node(size_t i) const
	{
		return true;
	}

	bool is_master(size_t i) const
	{
		return true;
	}

	bool is_master_on(size_t i, int pid) const
	{
		return pid == 0;
	}

	bool is_slave(size_t i, int OLlevel=0) const
	{
		return false;
	}

	bool i_must_assign(size_t i) const
	{
		return true;
	}

	bool i_can_set_coarse(size_t i) const
	{
		return nodes[i].could_be_coarse();
	}
	void calculate_unassigned()
	{
		m_iUnassigned = size();
	}
#endif
	void print()
	{
		UG_LOG("\n");
		for(size_t i=0; i<size(); i++)
			UG_LOG(info(i) << "\n");
		UG_LOG("\n");
	}

	void node_reset(size_t index)
	{
		if(nodes[index].is_assigned())
			m_iUnassigned++;
		if(nodes[index].is_aggressive_fine())
			m_iNrOfAggressiveFine--;
		if(nodes[index].is_coarse())
			m_iNrOfCoarse--;
		if(nodes[index].is_fine())
			m_iNrOfFine--;
		nodes[index].rating = 0;
	}

	/**
	 * this sets the node i fine.
	 * if it wasn't fine before and i_must_assign(i) is true,
	 * this reduces the number of unassigned nodes.
	 * asserts that i was not set coarse
	 * @param i
	 */
	void set_fine(size_t index)
	{
		UG_ASSERT(nodes[index].is_coarse() == false, "try to set node " << index << " fine, but is coarse???");
		if(nodes[index].is_fine() == false)
		{
			if(i_must_assign(index))
				m_iUnassigned--;
			m_iNrOfFine++;
		}
		nodes[index].set_fine();
		UG_DLOG(DBG_AMG_COARSENING, 4, "set node [" << index << "] fine\n");
	}

	/**
	 * this sets the node i aggressive fine.
	 * if i_must_assign(i) is true,
	 * this reduces the number of unassigned nodes.
	 * asserts that i was not set coarse
	 * @param i
	 */
	void set_aggressive_fine(size_t index)
	{
		//UG_ASSERT(nodes[index].is_fine() == false, "node " << index << " is already fine");
		if(i_must_assign(index))
			m_iUnassigned--;
		m_iNrOfAggressiveFine++;
		nodes[index].set_aggressive_fine();
		UG_DLOG(DBG_AMG_COARSENING, 4, "set node [" << index << "] aggressive fine\n");
	}

	/**
	 * set node to uncalculated fine
	 * @param index
	 */
	void external_set_uncalculated_fine(size_t index)
	{
		UG_ASSERT(!is_inner_node(index), "node " << index << " can not be set uncalculated_fine, since it is an inner node");
		//UG_ASSERT(!nodes[index].is_fine(), "node " << index << " is already fine?");
		if(nodes[index].is_fine() == false)
			m_iNrOfFine++;
		nodes[index].set_uncalculated_fine();
		m_iUnassigned--;
	}

	/**
	 * set node i coarse
	 * this function does not check whether we can
	 * set this node coarse
	 * @param i
	 */
	void external_set_coarse(size_t index)
	{
		if(nodes[index].is_coarse() == false)
		{
			nodes[index].set_coarse();
			m_iNrOfCoarse ++;

			if(i_must_assign(index))
				m_iUnassigned--;
//			P(index, index) = 1.0;
		}
		UG_DLOG(DBG_AMG_COARSENING, 4, "set node [" << index << "] coarse\n");
	}

	void external_uninterpolateable(size_t index)
	{
		nodes[index].set_uninterpolateable();
	}

	/**
	 * asserts that we can set this node coarse in a parallel sense.
	 * calls external_set_coarse
	 * @param index
	 */
	void set_coarse(size_t index)
	{
		UG_ASSERT(i_can_set_coarse(index), "Index " << index << " can not be coarse (" << nodes[index] << ").\n");
		external_set_coarse(index);
	}


	void set_dirichlet(size_t index)
	{
		//set_coarse(index); return;
		if(nodes[index].is_dirichlet() == false && i_must_assign(index))
			m_iUnassigned--;
		nodes[index].set_dirichlet();
		UG_DLOG(DBG_AMG_COARSENING, 4, "set node [" << index << "] dirichlet\n");
	}

	void set_uninterpolateable(size_t index)
	{
		UG_ASSERT(nodes[index].is_valid_rating() || nodes[index].is_uninterpolateable(),
				"node " << index << " is already set to " << nodes[index].rating);
		nodes[index].set_uninterpolateable();
		UG_DLOG(DBG_AMG_COARSENING, 4, "set node [" << index << "] uninterpolateable\n");
	}

	size_t get_nr_of_coarse()
	{
		return m_iNrOfCoarse;
	}

	size_t get_unassigned()
	{
		return m_iUnassigned;
	}

	size_t get_nr_of_fine()
	{
		return m_iNrOfFine;
	}

	size_t get_nr_of_aggressive_fine()
	{
		return m_iNrOfAggressiveFine;
	}

	size_t get_original_index(size_t i)
	{
		UG_ASSERT(m_level != (size_t)-1, "");
		return m_amghelper.GetOriginalIndex(m_level, i);
	}

//private:
public:
	stdvector<FAMGNode> nodes; // !!! this HAS to be a consecutive array, because it is used in the heap

private:
	size_t m_iNrOfCoarse;			// number of coarse nodes so far
	size_t m_iUnassigned;			// number of still unassigned nodes
	size_t m_iNrOfAggressiveFine; 	// number of aggressive fine nodes
	size_t m_iNrOfFine; 			// number of fine nodes

	size_t m_level;

	TProlongationMatrix &P;
	ug::cAMG_helper &m_amghelper;

#ifdef UG_PARALLEL
	stdvector<int> OLtype;
	std::map<size_t, int> m_masterOn;
	ParallelNodes &PN;
#endif
};


template<typename TNeighborstruct, typename prolongation_matrix_type>
void UpdateNeighbors(const cgraph &SymmNeighGraph, size_t node, stdvector<stdvector<TNeighborstruct> > &possible_neighbors,
		FAMGNodes<prolongation_matrix_type> &nodes, maxheap<FAMGNode> &heap)
{
	for(cgraph::const_row_iterator conn = SymmNeighGraph.begin_row(node); conn != SymmNeighGraph.end_row(node); ++conn)
	{
		size_t neigh = (*conn);
		UG_ASSERT(node != neigh, "");
		if(!nodes[neigh].is_valid_rating() || !nodes.i_must_assign(neigh))
			continue;

		if(nodes.update_rating(neigh, possible_neighbors[neigh]))
		{
			if(nodes[neigh].is_uninterpolateable())
			{
				UG_DLOG(LIB_ALG_AMG, 5, " remove from heap! ");
				heap.remove(neigh);
			}
			else
			{
				UG_DLOG(LIB_ALG_AMG, 5, " update in heap! ");
				heap.update(neigh);
			}
		}
	}
}

template<typename neighborstruct, typename heap_type, typename prolongation_matrix_type>
void GetRatings(stdvector<stdvector<neighborstruct> > &possible_neighbors,
		FAMGNodes<prolongation_matrix_type> &nodes, heap_type &heap)
{
	UG_DLOG(LIB_ALG_AMG, 4, "\nGetRatings...\n\n");
	for(size_t i=0; i<nodes.size(); i++)
	{
		UG_DLOG(LIB_ALG_AMG, 4, "node [" << i << "]: ");
		if(nodes.i_must_assign(i) && nodes[i].is_valid_rating())
		{
			nodes.update_rating(i, possible_neighbors[i]);
			if(nodes[i].is_valid_rating())
				heap.insert_item(i);
		}
		UG_DLOG(LIB_ALG_AMG, 4, "rating = " << nodes[i].get_rating() << "\n");
	}
}

// end group plugins_amg_famg
/// \}

}



#endif // __H__UG__LIB_DISC__FAMG_SOLVER__FAMG_NODEINFO_H__
