/*
 * named_infos.h
 *
 *  Created on: 29.08.2013
 *      Author: mrupp
 */

#ifndef NAMED_INFOS_H_
#define NAMED_INFOS_H_

template<typename TMethod>
class CMethod
{
public:
	CMethod(string _name, TMethod _func)
	{
		name = _name;
		func = _func;
		strRetValInfos = "";
		strParameterInfos = "";
		strTooltip = "";
		strHelp = "";
	}
	CMethod &parameter(string infos)
	{
		strParameterInfos = infos;
		return *this;
	}
	CMethod &tooltip(string tip)
	{
		strTooltip = tip;
		return *this;
	}
	CMethod &returnValue(string retval)
	{
		strRetValInfos = retval;
		return *this;
	}
	CMethod &help(string h)
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
	string name, strRetValInfos, strParameterInfos, strTooltip, strHelp;
};

template<typename TMethod>
CMethod<TMethod> Method(string name, TMethod func)
{
	return CMethod<TMethod>(name, func);
}

// usage like this
// .add_(Method("set_delta", &T::set_delta).parameter("delta").tooltip("\"Interpolation quality\" F may not be worse than this (F < m_delta). e.g. 0.5"))
#endif /* NAMED_INFOS_H_ */
