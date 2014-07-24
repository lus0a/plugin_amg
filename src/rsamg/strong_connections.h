/*
 * strong_connections.h
 *
 *  Created on: 13.05.2013
 *      Author: mrupp
 */

#ifndef STRONG_CONNECTIONS_H_
#define STRONG_CONNECTIONS_H_
#include "lib_algebra/common/graph/graph.h"
#include "lib_algebra/algebra_common/sparse_vector.h"
#include "rsamg_nodeinfo.h"
#include "rsamg_eliminate.h"
#include "rsamg_blockvalues.h"

namespace ug{



//! Adds 'strong negative connections' to graph.
/*! Criterion: $$-a_{ij} \ge \epsilon \max_{a_{ik}<0} |a_{ik}|$$ */
class StrongNegativeConnectionsByTrace
{
protected:
	double m_theta;
public:
	StrongNegativeConnectionsByTrace(double theta) : m_theta(theta){}

	template<typename TRowType>
	void find(const TRowType &Ai, size_t i, cgraph &graph)
	{
		typedef typename TRowType::const_iterator const_iterator;
		typedef typename TRowType::iterator iterator;
		typedef typename TRowType::value_type value_type;


		// rescale by diagonal
		value_type Bii=0.0;
		for(const_iterator conn = Ai.begin(); conn != Ai.end(); ++conn)
		{
			if (conn.index() == i) {
				GetInverse(Bii, conn.value()); break;
			}
		}


		for(const_iterator conn = Ai.begin(); conn != Ai.end(); ++conn)
		{ value_type &Aij = const_cast<value_type &> (conn.value());
			Aij = Bii * Aij; }


		// find maximum
		double dmax = 0.0;
		for(const_iterator conn = Ai.begin(); conn != Ai.end(); ++conn)
		{
			if(conn.index() == i) continue; // skip diag
			if(conn.value() != 0.0 && amg_diag(conn.value()) < dmax)
				dmax = amg_diag(conn.value());
		}
		// UG_LOG("Strong negative( = influencing) connections for " << i<< " (" << dmax <<") :" );

		// set connections
		for(const_iterator conn = Ai.begin(); conn != Ai.end(); ++conn)
		{
			if(conn.index() == i) continue; // skip diag
			const double d = amg_diag(conn.value());
			if( d <= m_theta * dmax) {
				graph.set_connection(i, conn.index());
				//	UG_LOG(conn.index() << "(" << d << "), ")
			}
		}
		//UG_LOG("\n")
	}
};


//! Adds 'strong positive connections' to graph.
/*! Criterion: $$a_{ij} \ge \epsilon^+ \max_{k\neq i} |a_{ik}|$$ */
class StrongPositiveConnections
{
protected:
	double m_theta;
public:
	StrongPositiveConnections(double theta) : m_theta(theta){}

	template<typename TRowType>
	void find(const TRowType &Ci, size_t i, cgraph &graph)
	{
		PROFILE_FUNC();
		double dmax = 0.0;
		typedef typename TRowType::const_iterator const_iterator;

		// find maximum
		for(const_iterator conn = Ci.begin(); conn != Ci.end(); ++conn)
		{
			if(conn->index() == i) continue; // skip diag
			const double dval =  fabs(amg_offdiag_value(conn->value()));
			if(dval > dmax) dmax = dval;
		}

		// set connections
		for(const_iterator conn = Ci.begin(); conn != Ci.end(); ++conn)
		{
			const size_t j =conn->index();
			if(j == i) continue; // skip diag
			const double dval =  amg_offdiag_value(conn->value());
			if( dval >= m_theta * dmax)
				graph.set_connection(i, j);
		}
	}
};


//! Adds 'strong negative connections' to graph.
/*! Criterion: $$-a_{ij} \ge \epsilon \max_{a_{ik}<0} |a_{ik}|$$ */
class StrongNegativeConnectionsByBlockNorm
{
protected:
	double m_theta;
public:
	StrongNegativeConnectionsByBlockNorm(double theta) : m_theta(theta){}

	template<typename TRowType>
	void find(const TRowType &Ci, size_t i, cgraph &graph)
	{
		double dmax = 0.0;
		typedef typename TRowType::const_iterator const_iterator;

		for(const_iterator conn = Ci.begin(); conn != Ci.end(); ++conn)
		{
			if(conn.index() == i) continue; // skip diag
			double d = BlockNorm(conn.value());
//			std::cout << d << std::endl;
			if(d > dmax) dmax = BlockNorm(conn.value());
		}

		for(const_iterator conn = Ci.begin(); conn != Ci.end(); ++conn)
		{
			if(conn.index() == i) continue; // skip diag
			if( BlockNorm(conn.value()) >= m_theta * dmax)
				graph.set_connection(i, conn.index());
		}
	}
};

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// CreateConnectionGraph:
//------------------------------
/**
 * Create graph of strong connections from a matrix A
 * i is strong coupled to j, if \f$ -A_{ij} \ge m_dTheta {\max\limits}_{\phi_1(A_{ik}) < 0} { \phi_2(A_{ik}) }\f$
 * at the moment, we use the function \f$\phi_1 = \phi_2 =\f$ amg_offdiag_value
 * \param	A			matrix A for which to calculate strong connectivity graph
 * \param 	graph		the calculated strong connectivity graph of A
 *						graph is afterwards made up of connections from a node i to j if
 *						j has a strong connection to i
 * \param	m_dTheta	the smaller, the more connections are considered strong e.g. 0.25 is a reasonable default value
 * \note	you can also pre-calculate some connectivity-measurement matrix C, and call this function with C.
 */


template<typename matrix_type>
void CreateStrongConnectionGraph(const matrix_type &A, cgraph &graph, double theta)
{
	AMG_PROFILE_FUNC();
	graph.resize(A.num_rows());


	StrongNegativeConnectionsByTrace strong(theta);

	for(size_t i=0; i< A.num_rows(); i++)
	{

		// Skip isolated node
		if(A.is_isolated(i)){
		//	UG_LOG(i << "is A-isolated!" << std::endl);
			continue;
		}

		// Eliminate positive row entries

		// 1) Extract & modify entries
		SparseVector<typename matrix_type::value_type> Ai_local(A.num_cols());
		PositiveOffDiagCheck<typename matrix_type::value_type> positiveEntries;
		EliminateEntries(positiveEntries, A, i, Ai_local);

		// 2) Find strong connections
		strong.find(Ai_local, i, graph);
	}

#ifdef AMG_PRINT_GRAPH
	graph.print();
#endif
}

template<typename matrix_type>
void CreateStrongConnectionGraphForSystems(const matrix_type &A, cgraph &graph, double theta)
{
	AMG_PROFILE_FUNC();
	graph.resize(A.num_rows());

	for(size_t i=0; i< A.num_rows(); i++)
	{
		// Skip isolated node
		if(A.is_isolated(i)) continue;

		// copy all off-diag entries
		SparseVector<typename matrix_type::value_type> Ri(A.num_cols());
		CopyOffDiagEntries(A, i, Ri);

		// Find strong connections
		StrongNegativeConnectionsByBlockNorm strong(theta);
		strong.find(Ri, i, graph);
	}
#ifdef AMG_PRINT_GRAPH
	graph.print();
#endif
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// CreateConnectionGraph:
//------------------------------
/**
 * Create graph of strong connections from a matrix A
 * i is strong coupled to j, if \f$ -A_{ij} \ge m_dTheta {\max\limits}_{\phi_1(A_{ik}) < 0} { \phi_2(A_{ik}) }\f$
 * at the moment, we use the function \f$\phi_1 = \phi_2 =\f$ amg_offdiag_value
 * \param	A			matrix A for which to calculate strong connectivity graph
 * \param 	graph		the calculated strong connectivity graph of A
 *						graph is afterwards made up of connections from a node i to j if
 *						j has a strong connection to i
 * \param	m_dTheta	the smaller, the more connections are considered strong e.g. 0.25 is a reasonable default value
 * \note	you can also pre-calculate some connectivity-measurement matrix C, and call this function with C.
 */


template<typename matrix_type>
void CreateFullConnectionGraph(const matrix_type &A, cgraph &graph, double theta)
{
	AMG_PROFILE_FUNC();

	const size_t nrows = A.num_rows();
	graph.resize(nrows);


	for(size_t i=0; i< nrows; ++i)
	{

		// skip isolated
		if(A.is_isolated(i)) continue;


		for(typename matrix_type::const_row_iterator conn = A.begin_row(i); conn != A.end_row(i); ++conn)
		{
			graph.set_connection(i, conn.index());
		}

	}

#ifdef AMG_PRINT_GRAPH
	graph.print();
#endif
}

inline void SetIsolatedCoarse(cgraph &graph, AMGNodes &nodes)
{
	for(size_t i=0; i< graph.size(); i++)
		if(graph.num_connections(i) == 0)
			nodes.set_coarse(i);
}

// Updates graph by including large positive connections
template<typename matrix_type>
void UpdateStrongConnectionGraph(const matrix_type &A, cgraph &graph, double theta_pos)
{
	AMG_PROFILE_FUNC();

	for(size_t i=0; i< A.num_rows(); i++)
	{
		// Skip isolated node
		if(A.is_isolated(i))
			continue;

		// Find strong connections
		StrongPositiveConnections strong(theta_pos);
		strong.find(A.get_row(i), i,graph);

	}
#ifdef AMG_PRINT_GRAPH
	graph.print();
#endif
}

}
#endif /* STRONG_CONNECTIONS_H_ */
