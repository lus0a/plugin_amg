/**
 * \file amg_misc.h
 *  Created on: 19.12.2011
 *      Author: mrupp
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
