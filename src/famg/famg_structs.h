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

#ifndef __H__LIB_ALGEBRA__AMG__FAMG_STRUCTS_H__
#define __H__LIB_ALGEBRA__AMG__FAMG_STRUCTS_H__
#include "famg_nodeinfo.h"
namespace ug{

/// \addtogroup plugins_amg_famg
/// \{

template<typename value_type, typename TIndex>
struct s_interpolation
{
	value_type value;
	TIndex from;
};

template<typename value_type> //, typename TIndex=size_t>
struct neighborstruct2
{
	template<typename T>
	void print(FAMGNodes<T> &rating)
	{
		UG_LOG(parents.size() << " parents, F = " << F << ": ");
		for(size_t i=0; i<parents.size(); i++)
			UG_LOG((i>0 ? "," : "") << "p" << i << ": " << parents[i].from << " [" << rating.get_original_index(parents[i].from) << "] -> " << parents[i].value);
		UG_LOG(std::endl);
	}
	FixedArray1<s_interpolation<value_type, size_t>, 2> parents;
	double F;
};

template<typename value_type=double> //, typename TIndex=size_t>
struct neighborstruct_var
{
	template<typename T>
	void print(FAMGNodes<T> &rating)
	{
		UG_LOG(parents.size() << " parents, F = " << F << ": ");
		for(size_t i=0; i<parents.size(); i++)
			UG_LOG((i>0 ? "," : "") << "p" << i << ": " << parents[i].from << " [" << rating.get_original_index(parents[i].from) << "] -> " << parents[i].value);
		UG_LOG(std::endl);
	}
	stdvector<s_interpolation<value_type, size_t> > parents;
	double F;
	std::vector<double> Fb ;
};

template<typename value_type>
std::ostream &operator <<(std::ostream &out, const neighborstruct_var<value_type> &n)
{
	out << n.parents.size() << " parents, F = " << n.F << ": [";
	for(size_t i=0; i<n.Fb.size(); i++)
		out << (i>0 ? "," : "")  << n.Fb[i];
	out << "]. P: ";
	for(size_t i=0; i<n.parents.size(); i++)
	{
		if(BlockNorm(n.parents[i].value) > 1e-8)
			out << " p" << i << ": " << n.parents[i].from << " -> " << n.parents[i].value;
	}
	out << ". ";
	return out;

}

//typedef neighborstruct2_t<double, size_t> neighborstruct2;
//typedef neighborstruct_var_t<double, size_t> neighborstruct_var;

// end group plugins_amg_famg
/// \}

}

#endif
