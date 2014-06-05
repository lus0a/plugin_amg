/*
 * named_infos.h
 *
 *  Created on: 29.08.2013
 *      Author: mrupp
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
