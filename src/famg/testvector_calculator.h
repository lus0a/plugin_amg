#ifndef TESTVECTOR_CALCULATOR_H_
#define TESTVECTOR_CALCULATOR_H_

#include <string>
#include <vector>

namespace ug{

template<typename TAlgebra>
class ITestvectorCalculator
{
	typedef typename TAlgebra::matrix_type matrix_type;
	typedef typename TAlgebra::vector_type vector_type;
	typedef MatrixOperator<matrix_type, vector_type> matrix_operator_type;
public:
	virtual void calc(SmartPtr<matrix_operator_type> A, MathVector<3> *pos, std::vector<vector_type> &testvectors) = 0;
	virtual std::string name() const  =0 ;
	virtual std::string config_string() const { return name(); }
	virtual bool needs_positions() const = 0;
	virtual ~ITestvectorCalculator() {}
};

}

#endif
