/*
 * Copyright (c) 2012-2013:  G-CSC, Goethe University Frankfurt
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


#ifndef AMG_MISC_H_
#define AMG_MISC_H_

#include "rsamg/rsamg_impl.h"
namespace ug
{

/// \addtogroup plugins_amg
/// \{

template<typename TMatrix>
void ReduceToStrongConnections(TMatrix &m1, const TMatrix &m2, double theta)
{
	PROFILE_FUNC();
	UG_LOG("Reducing matrix...");
	std::vector<typename TMatrix::connection> con;
	m1.resize_and_clear(m2.num_rows(), m2.num_cols());

#ifdef UG_PARALLEL
	m1.set_layouts(m2.layouts());
	m1.set_storage_type(m2.get_storage_mask());
#endif

	con.reserve(50);
	for(size_t i=0; i<m1.num_rows(); i++)
	{
		con.clear();
		double minConnValue, maxConnValue, diag;
		GetNeighborValues(m2, i, minConnValue, maxConnValue, diag);
		typename TMatrix::connection c;
		double barrier = theta*std::max(maxConnValue, dabs(minConnValue));
		typename TMatrix::value_type offLumpedSum=0;
		int iDiag=-1;
		for(typename TMatrix::const_row_iterator conn = m2.begin_row(i); conn != m2.end_row(i); ++conn)
		{
			if(conn.index() != i)
			{
				if(dabs(amg_offdiag_value(conn.value())) < barrier) // only strong
				{
					offLumpedSum += conn.value();
					continue;
				}
			}
			else
				iDiag = con.size();
			c.dValue = conn.value();
			c.iIndex = conn.index();
			con.push_back(c);
		}
		UG_ASSERT(iDiag != -1, "");
		con[iDiag].dValue += offLumpedSum;
		m1.set_matrix_row(i, &con[0], con.size());
	}
	UG_LOG("done.\n");
}

// end group plugins_amg
/// \}

}

#endif // AMG_MISC_H_
