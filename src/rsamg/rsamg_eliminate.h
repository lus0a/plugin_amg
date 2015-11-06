/*
 * Copyright (c) 2013-2014:  G-CSC, Goethe University Frankfurt
 * Author: Arne Nägel
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

#ifndef RSAMG_ELIMINATE_H_
#define RSAMG_ELIMINATE_H_

#include "rsamg_blockvalues.h"
namespace ug
{


template <class T> struct PositiveOffDiagCheck
{
	//! returns true for any positive off-diag connection
	static bool eval(size_t i, size_t j, T Aij) {
	   //return ((amg_offdiag_value(Aij) > 0.0) && (i!=j));
	   return ((amg_diag(Aij) > 0.0) && (i!=j));
	}
};


template <class T> class StrongOffDiagCheck
{
public:
	StrongOffDiagCheck(const cgraph &s) : strong(s) {};

	//! returns true for any strong connection
	 bool eval(size_t i, size_t j, T Aij) const {
		return ((strong.has_connection(i,j)) && (i!=j));
	}
private:
	 const cgraph &strong;
};



template<typename TMatrixType, typename TRowType>
void CopyOffDiagEntries(const TMatrixType &A, size_t i, TRowType &row, bool enforceNew=false)
{
	// copy matrix entries selected by rule
	for(typename TMatrixType::const_row_iterator connij = A.begin_row(i); connij != A.end_row(i); ++connij)
	{
		const size_t j=connij.index();
		row(connij.index()) = (i==j) ? 0.0 : connij.value();
	}
}


template<typename TMatrixType, typename TRowType, typename TElimRule>
void EliminateEntries(const TElimRule &elimRule, const TMatrixType &A, size_t i, TRowType &row, bool enforceNew=false)
{
	// copy matrix entries we wish to keep
	for(typename TMatrixType::const_row_iterator connij = A.begin_row(i); connij != A.end_row(i); ++connij)
	{
		size_t j=connij.index();
		typename TMatrixType::value_type Aij=connij.value();
//		UG_LOG( i << " " << connij.index() << "\n");
		row(connij.index()) = (elimRule.eval(i, j, Aij)) ? 0.0 : connij.value();
	}

	/*
	std::cout << "EliminateEntries (BEGIN) "<< i << "\n";
	row.print();
	*/

	// eliminate all 'positive' row entries
	for(typename TMatrixType::const_row_iterator conn = A.begin_row(i); conn != A.end_row(i); ++conn)
	{
		const size_t j=conn.index();
		typename TMatrixType::value_type Aij=conn.value();


		if (elimRule.eval(i,j,Aij))
		{
			// eliminate entry Aij
			typename TMatrixType::value_type Bjj;
			typename TMatrixType::const_row_iterator connjj = A.get_connection(j,j);
			GetInverse(Bjj, connjj.value());

			typename TMatrixType::value_type Cij=Aij*Bjj;			

			//	std::cout << "Eliminating A (" << i << "," << j<< ")=" << Aij << "using" <<  Bjj << Cij << std::endl;

			for(typename TMatrixType::const_row_iterator connjk = A.begin_row(j); connjk != A.end_row(j); ++connjk)
			{
				const size_t k=connjk.index();
				if (k==j) continue;
				
			
				typename TMatrixType::value_type Aijjk = Cij * connjk.value();
				//	std::cout << "adding  (" << i << "," << k<< ")=" << -Aijjk << " = -Cij*"<< connjk.value() <<" to "<< row(k) << std::endl;

				if(enforceNew || row.has_connection(k))
					row(k) -= Aijjk;
			}
		}

	}


	/*	std::cout << "EliminateEntries (POST) "<< i << "\n";
	row.print();
	std::cout << "EliminateEntries (END) "<< i << "\n";
	*/
}
}
#endif /* RSAMG_ELIMINATE_H_ */
