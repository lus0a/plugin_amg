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
