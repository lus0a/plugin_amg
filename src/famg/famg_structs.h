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
