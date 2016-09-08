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
	typedef std::vector<int> ind_vec_type;
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
	{
		std::map<double, size_t>::const_iterator it=m_map[tdim].find(pos);
		UG_ASSERT(it!=m_map[tdim].end(), "ERROR: Map invalid, no element for x="<< pos);
		return it->second;
	}

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


	// determines neighbor coordinates (N1)
	size_t get_node_nbcoord(const MathVector<3> &center, MathVector<3> nbcoord[])
	{
	  // indices of all neighbors
	  ind_vec_type nbtensor[3];
	  const size_t dim =3;
	  for (unsigned int i=0; i < dim; ++i)
	    {
	       int indexi = get_axis_index(i, center[i]);
	       for (size_t ii=std::min(indexi-1,0); ii<m_size[i];++ii)
		 { nbtensor[i].push_back(ii);}
	    }	
	  
	  // loop neighbors
	  size_t a = 0;
	  for (ind_vec_type::const_iterator it0= nbtensor[0].begin(); it0!=nbtensor[0].end(); ++it0)
	    for (ind_vec_type::const_iterator it1= nbtensor[1].begin(); it1!=nbtensor[1].end(); ++it1)
	      for (ind_vec_type::const_iterator it2= nbtensor[2].begin(); it2!=nbtensor[2].end(); ++it2)
		{
		  size_t indexa[3];
		  indexa[0] = *it0;
		  indexa[1] = *it1;
		  indexa[2] = *it2;
		  get_node_coord<3>(indexa, nbcoord[a]);
		  a++;
		}

	  return a;
	  
	}

protected:
	pos_vec_type m_coord[3];
	map_type m_map[3];
	size_t m_size[3];
};



class TensorLevelInformation
: protected TensorPositionInformation
{
public:
	typedef std::vector<int> level_vec_type;

	template <unsigned int dim>
	size_t get_node_level(const MathVector<dim> &coord) const
	{
		int level = m_level[0][get_axis_index(0, coord[0])];
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



	// determines neighbor coordinates (N1)
	size_t get_parents(const MathVector<3> &center,
				int level,
				MathVector<3> pcoord[],
				double pweight[])
	{
	  // find coarse neighbors
	 
	  const size_t dim =3;
	  ind_vec_type nbtensor[dim];
	  int index[dim];      // indices of center
	  
	  int ncoarse[dim];    // number of coarse nodes in this direction
	  double weight[dim];  // weight in this direction
	 
	  for (unsigned int i=0; i < dim; ++i)
	    {
	    int ii = get_axis_index(i, center[i]);
	    index[i] = ii;
	    ncoarse[i]=0;
	    //const int endi = std::min((int)m_size[i]-1, index[i]+1);
	    //for (int nbi=std::max(index[i]-1,0); nbi<=endi; ++nbi)
	    
	    if (m_level[i][ii]>level)
	    {
	      // found coarse node
	      nbtensor[i].push_back(ii); 
	      ncoarse[i]++; 
	    }
	    else
	    {
	      // backward search
	      for (int nbi=ii-1; nbi>=0; --nbi)
	      {
	    	  if (m_level[i][nbi]>level) // found coarse node
	    	  { ncoarse[i]++; nbtensor[i].push_back(nbi); }

	    	  if (m_level[i][nbi]>=level) break;
	      }

	      // forward search
	      for (size_t nbi=ii+1; nbi<m_size[i]; ++nbi)
	      {
	    	  if (m_level[i][nbi]>level)
	    	  { ncoarse[i]++; nbtensor[i].push_back(nbi); }

	    	  if (m_level[i][nbi]>=level) break;
		 
	      }
	    }

	  }	
	  
	  // effective directions (DEBUG)
	  int directions=0; 
	  for (unsigned int i=0; i < dim; ++i)
	  {
	    if (ncoarse[i]>0)
	      directions++;
	  }
	  /*
	  if ((directions<=0)|| directions>3)
	    {
	      std::cerr << "ERROR:\n";
		for (unsigned int i=0; i < dim; ++i)
		{
		  for (int nbi=index[i]-1; nbi>=0; --nbi)
		  {
		    if (m_level[i][nbi]>level) // found coarse node
		    { std::cerr << i << ":" <<nbi <<std::endl;}
		    if (m_level[i][nbi]>=level) break;
		  }
		}
	    }*/
	  UG_ASSERT( ((directions<=3) && (directions>0)), "Huhh: Directions do not match!");
	  
	  // assign weights
	  for (unsigned int i=0; i < dim; ++i)
	  {
	    //	    if (ncoarse[i]==0) {
	    //  nbtensor[i].push_back(index[i]);
	    //}
	    UG_ASSERT(ncoarse[i]>0, "No coarse dim found: "<< i);
	    weight[i] = 1.0/nbtensor[i].size(); 
      }

	  // create nodes
	  size_t indexa[3];
	  size_t a = 0;
	  double sumw =0.0;

	  //std::cerr << "get_parents:" << center << "["<< index[0]<< "," << index[1]<< "," << index[2]<< "]->";
	  for (ind_vec_type::const_iterator it0= nbtensor[0].begin(); it0!=nbtensor[0].end(); ++it0)   
	    for (ind_vec_type::const_iterator it1= nbtensor[1].begin(); it1!=nbtensor[1].end(); ++it1)
	      for (ind_vec_type::const_iterator it2= nbtensor[2].begin(); it2!=nbtensor[2].end(); ++it2)
		{
		  indexa[0] = *it0;
		  indexa[1] = *it1;
		  indexa[2] = *it2;
		  get_node_coord<3>(indexa, pcoord[a]);
		//  std::cerr << pcoord[a] << "["<< indexa[0]<< "," << indexa[1]<< "," << indexa[2]<< "], ";
		  pweight[a] = weight[0]*weight[1]*weight[2];
		  

		  sumw +=pweight[a];
		  a++;		
		}

	  if (sumw>1.0)
	  {
	    std::cerr << sumw << " " << a << " ";

	    for (unsigned int i=0; i < a; ++i)
	    {
	      std::cerr << pweight[i] << " ";
            }

	    for (unsigned int i=0; i < dim; ++i)
	    {
	      std::cerr << weight[i] << " ";
	    }
	    
	    UG_ASSERT(sumw<=1.0, "Sum does not match!");
	  }
	  return a;
	  
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
			double hI=1.0, double hJ=1.0, double hK=1.0);

	SimpleTensorPositionProvider(size_t I, size_t J, size_t K,
				double hI, double hJ, double hK,
				double cx, double cy, double cz);

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
