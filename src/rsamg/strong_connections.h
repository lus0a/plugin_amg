/*
 * strong_connections.h
 *
 *  Created on: 13.05.2013
 *      Author: mrupp
 */

#ifndef STRONG_CONNECTIONS_H_
#define STRONG_CONNECTIONS_H_


namespace ug{

//! Adds 'strong negative connections' to graph.
/*! Criterion: $$-a_{ij} \ge \epsilon \max_{a_{ik}<0} |a_{ik}|$$ */
template<typename TRowType>
void FindStrongNegativeConnections(const TRowType &Ci, size_t i, double theta, cgraph &graph)
{
	double dmax = 0.0;
	typedef typename TRowType::const_iterator const_iterator;

	for(const_iterator conn = Ci.begin(); conn != Ci.end(); ++conn)
	{
		if(conn.index() == i) continue; // skip diag
		if(conn.value() != 0.0 && amg_trace(conn.value()) < dmax)
			dmax = amg_trace(conn.value());
	}
	UG_LOG("Strong negative( = influencing) connections for " << i<< " (" << dmax <<") :" );
	for(const_iterator conn = Ci.begin(); conn != Ci.end(); ++conn)
	{
		if(conn.index() == i) continue; // skip diag
		const double d = amg_trace(conn.value());
		if( d <= theta * dmax) {
			graph.set_connection(i, conn.index());
			UG_LOG(conn.index() << "(" << d << "), ")
		}
	}
	UG_LOG("\n")
}


//! Adds 'strong positive connections' to graph.
/*! Criterion: $$a_{ij} \ge \epsilon^+ \max_{k\neq i} |a_{ik}|$$ */
template<typename TRowType>
void FindStrongPositiveConnections(const TRowType &Ci, size_t i, double theta, cgraph &graph)
{
	PROFILE_FUNC();
	double dmax = 0.0;
	typedef typename TRowType::const_iterator const_iterator;

	for(const_iterator conn = Ci.begin(); conn != Ci.end(); ++conn)
	{
		if(conn->index() == i) continue; // skip diag
		const double dval =  fabs(amg_offdiag_value(conn->value()));
		if(dval > dmax) dmax = dval;
	}

	for(const_iterator conn = Ci.begin(); conn != Ci.end(); ++conn)
	{
		const size_t j =conn->index();
		if(j == i) continue; // skip diag
		const double dval =  amg_offdiag_value(conn->value());
		if( dval >= theta * dmax)
			graph.set_connection(i, j);
	}
}

//! Adds 'strong negative connections' to graph.
/*! Criterion: $$-a_{ij} \ge \epsilon \max_{a_{ik}<0} |a_{ik}|$$ */
template<typename TRowType>
void FindStrongConnectionsByBlockNorm(const TRowType &Ci, size_t i, double theta, cgraph &graph)
{
	double dmax = 0.0;
	typedef typename TRowType::const_iterator const_iterator;

	for(const_iterator conn = Ci.begin(); conn != Ci.end(); ++conn)
	{
		if(conn.index() == i) continue; // skip diag
		double d = BlockNorm(conn.value());
		std::cout << d << std::endl;
		if(d > dmax) dmax = BlockNorm(conn.value());
	}

	for(const_iterator conn = Ci.begin(); conn != Ci.end(); ++conn)
	{
		if(conn.index() == i) continue; // skip diag
		if( BlockNorm(conn.value()) >= theta * dmax)
			graph.set_connection(i, conn.index());
	}
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
void CreateStrongConnectionGraph(const matrix_type &A, cgraph &graph, double theta, AMGNodes &nodes)
{
	AMG_PROFILE_FUNC();
	graph.resize(A.num_rows());

	for(size_t i=0; i< A.num_rows(); i++)
	{

		// Skip isolated node
		if(A.is_isolated(i))
			continue;

		// Eliminate positive row entries
		SparseVector<typename matrix_type::value_type> Ri(A.num_cols());
		PositiveOffDiagCheck<typename matrix_type::value_type> positiveEntries;
		EliminateEntries(positiveEntries, A, i, Ri);

		// Find strong connections
		FindStrongNegativeConnections(Ri, i, theta, graph);

		if(graph.num_connections(i) == 0)
			nodes.set_coarse(i);
	}

#ifdef AMG_PRINT_GRAPH
	graph.print();
#endif
}

template<typename matrix_type>
void CreateStrongConnectionGraphForSystems(const matrix_type &A, cgraph &graph, double theta, AMGNodes &nodes)
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
		FindStrongConnectionsByBlockNorm(Ri, i, theta, graph);

		if(graph.num_connections(i) == 0)
			nodes.set_coarse(i);
	}

#ifdef AMG_PRINT_GRAPH
	graph.print();
#endif
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
		FindStrongPositiveConnections(A.get_row(i), i, theta_pos, graph);

	}
#ifdef AMG_PRINT_GRAPH
	graph.print();
#endif
}

}
#endif /* STRONG_CONNECTIONS_H_ */
