/*
 * Copyright (c) 2014:  G-CSC, Goethe University Frankfurt
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

#ifndef ARRAY_BINARY_OPERATION_SCHEME_H_
#define ARRAY_BINARY_OPERATION_SCHEME_H_

#include "lib_algebra/parallelization/communication_scheme.h"
#include <boost/lambda/lambda.hpp>

namespace ug{

template<typename TArray, class BinaryOperation>
class ArrayBinaryOperationCommunicationScheme :
		public CommunicationScheme<ArrayBinaryOperationCommunicationScheme<TArray, BinaryOperation>,
			typename TArray::value_type >
{
	typedef typename TArray::value_type value_type;
public:
	ArrayBinaryOperationCommunicationScheme(TArray &t, BinaryOperation op) : m_arr(t), m_op(op)
	{
		boost::ignore_unused_variable_warning(boost::lambda::_1);
		boost::ignore_unused_variable_warning(boost::lambda::_2);
		boost::ignore_unused_variable_warning(boost::lambda::_3);
	}

	inline const value_type &send(int pid, size_t index) const
	{
		return m_arr[index];
	}

	inline void receive(int pid, size_t index, value_type &v)
	{
		m_arr[index] = m_op(m_arr[index], v);
	}

	inline int get_element_size() const
	{
		if(block_traits<value_type>::is_static) return sizeof(value_type);
		else
			return -1;
	}
private:
	TArray &m_arr;
	BinaryOperation m_op;
};


template<typename T>
StdArrayCommunicationScheme<T> GetStdArrayCommunicationScheme(T &vec)
{
	return StdArrayCommunicationScheme<T>(vec);
}


template<typename T, class BinaryOperation>
void SendVecFromTo(pcl::InterfaceCommunicator<IndexLayout> &communicator,
		const IndexLayout &sendingLayout, const IndexLayout &receivingLayout, T &vec, BinaryOperation op)
{
	UG_LOG("CommunicateOnInterfaces");
	UG_LOG("sendingLayout: " << sendingLayout << "\n");
	UG_LOG("receivingLayout: " << receivingLayout << "\n");
	ArrayBinaryOperationCommunicationScheme<T, BinaryOperation> scheme(vec, op);
	CommunicateOnInterfaces(communicator, sendingLayout, receivingLayout, scheme);
}


template<typename T>
void SendVecFromTo(pcl::InterfaceCommunicator<IndexLayout> &communicator,
		const IndexLayout &sendingLayout, const IndexLayout &receivingLayout, T &vec)
{
	SendVecFromTo(communicator, sendingLayout, receivingLayout, vec, boost::lambda::_2);
}


}
#endif /* ARRAY_BINARY_OPERATION_SCHEME_H_ */
