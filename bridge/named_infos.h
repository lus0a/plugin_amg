/*
 * Copyright (c) 2013-2014:  G-CSC, Goethe University Frankfurt
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

#ifndef NAMED_INFOS_H_
#define NAMED_INFOS_H_

namespace ug{
namespace bridge{
/*
class ParameterDouble
{
	ParameterDouble(std::string name)
	{
		m_name = name;
		m_min = m_max = m_value = NAN;
	}
	ParameterDouble &min(double m)
	{
		m_min = m;
		return *this;
	}
	ParameterDouble &max(double m)
	{
		m_max = m;
		return *this;
	}

	ParameterDouble &value(double m)
	{
		m_value = m;
		return *this;
	}

	operator std::string ()
	{
		Stringify s;
		s << name << "|Double|default";
		if(!isnan(m_min)) 	s << "|min=" << m_min << "D";
		if(!isnan(m_max)) 	s << "|max=" << m_max << "D";
		if(!isnan(m_value))	s << "|value=" << m_value << "D";
	}

	std::string m_name;
	double m_min,m_max, m_value;
};

class ParameterInteger
{
	// ..
};

class ParameterDoubleSlider : public ParameterDouble
{
	// ...
};
class ParameterIntegerSelection
{
 // min, max
};

class ParameterStringFileDialog
{
	void ending(std::string s)
	{
		//..
	}

	void description(std::string s)
	{
		//..
	}

	std::string tostring()
	{
		//..
	}
};
class ParameterString
{
	// ..
};
class ParameterStringLoadFileDialog : public ParameterStringFileDialog
{
public:
	std::string m_name;
	ParameterLoadFileDialog(std::string name) : m_name(name) {}
	std::string operator std::string()
	{
		return Stringify() << m_name << "|load-dialog" << ParameterFileDialog::tostring();
	}
};


class ParameterStringSaveFileDialog : public ParameterStringSaveDialog
{
public:
	std::string m_name;
	ParameterLoadFileDialog(std::string name) : m_name(name) {}
	std::string operator std::string()
	{
		return Stringify() << m_name << "|save-dialog" << ParameterFileDialog::tostring();
	}
};

class ParameterStringSelection
{
	// selection
};

void bla()
{
	bla.parameters(
				ParameterNumber("steps").min(0).max(100).value(2))
			<<  ParameterNumber("steps2").min(0).max(100).value(2))
			<<  ParameterNumber("steps").min(0).max(100).value(2));

}
*/


template<typename TMethod>
class CMethod
{
public:
	CMethod(std::string _name, TMethod _func)
	{
		name = _name;
		func = _func;
		strRetValInfos = "";
		strParameterInfos = "";
		strTooltip = "";
		strHelp = "";
	}
	CMethod &parameter(std::string infos)
	{
		strParameterInfos = infos;
		return *this;
	}
	CMethod &tooltip(std::string tip)
	{
		strTooltip = tip;
		return *this;
	}
	CMethod &returnValue(std::string retval)
	{
		strRetValInfos = retval;
		return *this;
	}
	CMethod &help(std::string h)
	{
		strHelp = h;
		return *this;
	}

	template<typename T>
	ExportedClass<T> &add_to(ExportedClass<T> &c)
	{
		return c.add_method(name, func, strRetValInfos, strParameterInfos, strTooltip, strHelp);
	}
private:
	TMethod func;
	std::string name, strRetValInfos, strParameterInfos, strTooltip, strHelp;
};

template<typename TMethod>
CMethod<TMethod> Method(std::string name, TMethod func)
{
	return CMethod<TMethod>(name, func);
}

}}
// usage like this
// .add_(Method("set_delta", &T::set_delta)
//		.parameter("delta")
//		.tooltip("\"Interpolation quality\" F may not be worse than this (F < m_delta). e.g. 0.5"))
#endif /* NAMED_INFOS_H_ */
