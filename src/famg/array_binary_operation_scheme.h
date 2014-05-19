/*
 * array_binary_operation_scheme.h
 *
 *  Created on: 23.04.2014
 *      Author: mrupp
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
