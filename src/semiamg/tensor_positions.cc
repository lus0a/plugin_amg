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


#include "tensor_positions.h"

namespace ug {

/*
bool TensorPositionProvider::get_positions(std::vector<MathVector<2> >&vec)
{
	vec.resize(m_nvec);
	size_t ind=0;

	for (size_t j=0; j<m_J; j++)
		for (size_t i=0; i<m_I; i++)
					{
						vec[ind][0] = m_tpi;
						vec[ind][1] = m_tpi;
						vec[ind][2] = 0;
						ind++;
					}
	return (ind==m_nvec);
};

bool TensorPositionProvider::get_positions(std::vector<MathVector<3> >&vec)
{
	vec.resize(m_nvec);
	size_t ind=0;
	for (size_t k=0; k<m_K; k++)
		for (size_t j=0; j<m_J; j++)
			for (size_t i=0; i<m_I; i++)
			{
					vec[ind][0] = m_tpi.get_coord(0, i);
					vec[ind][1] = m_tpi;
					vec[ind][2] = m_tpi;
					ind++;
			}
	return (ind==m_nvec);
};
*/


SimpleTensorPositionProvider::SimpleTensorPositionProvider(size_t I, size_t J, size_t K,
			double hI, double hJ, double hK) :
	m_nvec(I*J*K), m_I(I), m_J(J), m_K(K)
	{
		m_h[0]=hI; m_h[1]=hJ; m_h[2]=hK;
		set_coord_offset(0.0, 0.0, 0.0);
	}

SimpleTensorPositionProvider::SimpleTensorPositionProvider(size_t I, size_t J, size_t K,
				double hI, double hJ, double hK,
				double cx, double cy, double cz) :
		m_nvec(I*J*K), m_I(I), m_J(J), m_K(K)
		{
			m_h[0]=hI; m_h[1]=hJ; m_h[2]=hK;
			set_coord_offset(cx, cy, cz);
		}


bool SimpleTensorPositionProvider::get_positions(std::vector<MathVector<2> >&vec)
{
	vec.resize(m_nvec);
	size_t ind=0;

	for (size_t j=0; j<m_J; j++)
		for (size_t i=0; i<m_I; i++)
					{
						vec[ind][0] = i*m_h[0] + m_coord_offset[0];
						vec[ind][1] = j*m_h[1] + m_coord_offset[1];
						vec[ind][2] = m_coord_offset[2];
						ind++;
					}
	return (ind==m_nvec);
};

bool SimpleTensorPositionProvider::get_positions(std::vector<MathVector<3> >&vec)
{
	vec.resize(m_nvec);
	size_t ind=0;
	for (size_t k=0; k<m_K; k++)
		for (size_t j=0; j<m_J; j++)
			for (size_t i=0; i<m_I; i++)
			{
					vec[ind][0] = i*m_h[0] + m_coord_offset[0];
					vec[ind][1] = j*m_h[1] + m_coord_offset[1];
					vec[ind][2] = k*m_h[2] + m_coord_offset[2];
					ind++;
			}
	return (ind==m_nvec);
};


}
