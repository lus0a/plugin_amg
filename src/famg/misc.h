/*
 * vector_misc.h
 *
 *  Created on: 19.05.2014
 *      Author: mrupp
 */

#ifndef VECTOR_MISC_H_
#define VECTOR_MISC_H_

#include <vector>
#include "lib_algebra/small_algebra/storage/fixed_array.h"
#include "lib_algebra/small_algebra/blocks.h"

namespace ug{

class two_indices : public FixedArray1<size_t, 2>
{
	using FixedArray1<size_t, 2>::at;
public:
	two_indices(size_t i0, size_t i1)
	{
		at(0) = i0;
		at(1) = i1;
	}
};
class three_indices : public FixedArray1<size_t, 3>
{
	using FixedArray1<size_t, 3>::at;
public:
	three_indices(size_t i0, size_t i1, size_t i2)
	{
		at(0) = i0;
		at(1) = i1;
		at(2) = i2;
	}
};

class index_range
{
public:
	index_range(size_t s) :
			m_size(s)
	{
	}
	inline size_t size() const
	{
		return m_size;
	}
	inline size_t operator[](size_t i) const
	{
		return i;
	}
private:
	const size_t m_size;
};


class SubVector
{
public:
	SubVector(std::vector<size_t> &t, size_t subsize) : v(t)
	{
		m_subsize = subsize;
	}

	inline size_t size() const
	{
		return m_subsize;
	}
	inline size_t operator[](size_t i) const
	{
		return v[i];
	}
private:
	size_t m_subsize;
	std::vector<size_t> &v;
};


template<typename T>
void BlockToDouble(DenseMatrix<VariableArray2<double> > &Ad, T &A)
{
	size_t blockSize = GetRows(A(0,0));
	Ad.resize(blockSize*A.num_rows(), blockSize*A.num_cols());

	for(size_t i=0; i<A.num_rows(); i++)
	{
		for(size_t j=0; j<A.num_cols(); j++)
			for(size_t r=0; r<blockSize; r++)
				for(size_t c=0; c<blockSize; c++)
					Ad(i*blockSize+r, j*blockSize+c) = BlockRef(A(i, j), r, c);
	}
}

template<typename T>
void DoubleToBlock(T &A, DenseMatrix<VariableArray2<double> > &Ad, size_t blockSize)
{
	A.resize(Ad.num_rows()/blockSize, Ad.num_cols()/blockSize);

	for(size_t i=0; i<Ad.num_rows(); i++)
		for(size_t j=0; j<Ad.num_cols(); j++)
			BlockRef(A(i/blockSize, j/blockSize), i%blockSize, j%blockSize) = Ad(i, j);
}


template<typename TA, typename TCompareValues>
void SortAlikeB(std::vector<TA> &vA, const TCompareValues &vB)
{
	if(vA.size() != vB.size()) throw "size mismatch!";
	std::vector<size_t> indices(vB.size());
	GetSortedIndices(indices, vB);
	std::vector<TA> vA2(vB.size());
	for(size_t i=0; i<indices.size(); i++)
		vA2[indices[i]] = vA[i];
	vA = vA2;
}
}
#endif /* VECTOR_MISC_H_ */
