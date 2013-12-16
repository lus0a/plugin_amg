/*
 * tensor_positions.h
 *
 *  Created on: 26.11.2013
 *      Author: anaegel
 */

#ifndef TENSOR_POSITIONS_H_
#define TENSOR_POSITIONS_H_

//#include <cstddef>
//#include <cstring>
#include <vector>


#include "ug.h"
#include "lib_algebra/operator/debug_writer.h"
#include "common/math/math_vector_matrix/math_vector.h"

namespace ug {
class TensorPositionInformation{
public:
	typedef std::vector<double> pos_vec_type;
	typedef std::vector<size_t> ind_vec_type;
	typedef std::map<double, size_t> map_type;

	// indices: 0 ... pos.size()-1
	void set_axis_coord(unsigned int tdim, const pos_vec_type &pos)
	{
		m_coord[tdim] = pos;
		map_type &mymap = m_map[tdim];

		pos_vec_type::const_iterator piter = pos.begin();
		size_t index=0;//offset;
		for (; piter!=pos.end(); piter++, index++)
		{
			mymap[(*piter)] = index;

		}
		m_size[tdim] = pos.size();
	}

	// user supplied indices (positions and indices)
	void set_axis_coord(unsigned int tdim, const pos_vec_type &pos, const ind_vec_type &ind)
	{
		m_coord[tdim].reserve(pos.size());
		map_type &map = m_map[tdim];
		pos_vec_type::const_iterator piter = pos.begin();
		ind_vec_type::const_iterator iiter = ind.begin();
		for (; piter!=pos.end(); piter++, iiter++)
		{
			m_coord[tdim][*iiter] = (*piter);
			map[(*piter)] = (*iiter);
		}

		m_size[tdim] = pos.size();
	}

	/// maps: pos -> index_axis (for axis)
	size_t get_axis_index (unsigned int tdim, double pos) const
	{ return (m_map[tdim].find(pos)->second);}

	/// maps:coord -> (triple value) index
	template <unsigned int dim>
	bool get_node_index(const MathVector<dim> &coord, size_t *index)
	{
		for (unsigned int i=0; i<dim; i++)
		{index[i] = get_axis_index(i, coord[i]);}
		return true;
	}



	/// coordinates
	double get_axis_coord (unsigned int tdim, size_t i) const
	{ return m_coord[tdim][i];}

	template <unsigned int dim>
	bool get_node_coord(const size_t index[], MathVector<dim> &coord)
	{
		for (unsigned int i=0; i<dim; i++)
		{coord[i] = get_axis_coord(i, index[i]);}
		return true;
	}




	size_t size() const
	{ return m_size[0]*m_size[1]*m_size[2]; }

	/*
	/// return (single value) index;
	size_t get_index(size_t *triple)
	{
		size_t index = triple[2]* m_size[1];
		index += triple[1];
		index *= m_size[0];
		index += triple[0];
		return index;
	}

	/// return (single value) index;
	size_t get_index(const MathVector<dim> &coord)
	{
		size_t triple[3];
		triple[0] = get_axis_index(0, coord[0]);
		triple[1] = get_axis_index(1, coord[1]);
		triple[2] = get_axis_index(2, coord[2]);
		return get_index(triple);
	}
*/

private:
	pos_vec_type m_coord[3];
	map_type m_map[3];
	size_t m_size[3];
};



class TensorLevelInformation
: protected TensorPositionInformation
{
public:
	typedef std::vector<size_t> level_vec_type;

	template <unsigned int dim>
	size_t get_node_level(const MathVector<dim> &coord) const
	{
		size_t level = m_level[0][get_axis_index(0, coord[0])];
		for (unsigned int i=1; i<dim; i++)
		{ level = std::min(level, m_level[i][get_axis_index(i, coord[i])]); }
		return level;
	}

	template <unsigned int dim>
	void get_node_levels(const MathVector<dim> &coord, size_t *level) const
	{
		for (unsigned int i=0; i<dim; i++)
		{ level[i] = m_level[i][get_axis_index(i, coord[i])]; }
	}


	void set_axis_level(unsigned int tdim, const level_vec_type &lev, const pos_vec_type &pos, const ind_vec_type &ind)
	{
		TensorPositionInformation::set_axis_coord(tdim, pos, ind);
		m_level[tdim] = lev;
	}

	void set_axis_level(unsigned int tdim, const level_vec_type &lev, const pos_vec_type &pos)
	{
		TensorPositionInformation::set_axis_coord(tdim, pos);
		m_level[tdim] = lev;
	}

private:
	level_vec_type m_level[3];

};


/*
class TensorPositionProvider
: public IPositionProvider<2>, public IPositionProvider<3>
{
public:

	TensorPositionProvider(TensorPositionInformation &tpi) m_tpi(tpi){}
	virtual TensorPositionProvider~() {};
	virtual bool get_positions(std::vector<MathVector<2> >&vec);
	virtual bool get_positions(std::vector<MathVector<3> >&vec);

private:
	TensorPositionInformation & m_tpi;
};
*/
class SimpleTensorPositionProvider
: public IPositionProvider<2>, public IPositionProvider<3>
{
public:

	SimpleTensorPositionProvider(size_t I, size_t J, size_t K,
			double hI=1.0, double hJ=1.0, double hK=1.0) :
	m_nvec(I*J*K), m_I(I), m_J(J), m_K(K)
	{
		m_h[0]=hI; m_h[1]=hJ; m_h[2]=hK;
		set_coord_offset(0.0, 0.0, 0.0);
	}

	SimpleTensorPositionProvider(size_t I, size_t J, size_t K,
				double hI, double hJ, double hK,
				double cx, double cy, double cz) :
		m_nvec(I*J*K), m_I(I), m_J(J), m_K(K)
		{
			m_h[0]=hI; m_h[1]=hJ; m_h[2]=hK;
			set_coord_offset(cx, cy, cz);
		}



	virtual ~SimpleTensorPositionProvider() {};

	//void init(void *g, size_t nvec) {};

	virtual bool get_positions(std::vector<MathVector<2> >&vec);
	virtual bool get_positions(std::vector<MathVector<3> >&vec);

private:
	void set_coord_offset(double cx, double cy, double cz)
	{m_coord_offset[0] = cx; m_coord_offset[1] = cy; m_coord_offset[2] = cz;}

	size_t m_nvec;
	size_t m_I, m_J, m_K;
	double m_h[3];
	double m_coord_offset[3];
};



}
#endif /* TENSOR_POSITIONS_H_ */
