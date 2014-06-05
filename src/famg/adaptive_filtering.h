/*
 * adaptive_fitlering.h
 *
 *  Created on: 21.05.2014
 *      Author: mrupp
 */

#include "testvector_calculators.h"

#ifndef ADAPTIVE_FITLERING_H_
#define ADAPTIVE_FITLERING_H_

namespace ug{

template<typename TAlgebra>
class AdaptiveFilteringTestvectors : public ITestvectorCalculator<TAlgebra>
{
	typedef typename TAlgebra::matrix_type matrix_type;
	typedef typename TAlgebra::vector_type vector_type;
	typedef MatrixOperator<matrix_type, vector_type> matrix_operator_type;

	SmartPtr<ILinearIterator<vector_type> > m_smoother;
	size_t m_numSmooths;
	size_t m_numTestvectors;

	SmartPtr<ITestvectorCalculator<TAlgebra> > m_startTV;
	std::vector<vector_type> additionalTV;

public:
	AdaptiveFilteringTestvectors(SmartPtr<ITestvectorCalculator<TAlgebra> > startTV )
	{
		m_startTV = startTV;
	}

	virtual void calc(SmartPtr<matrix_operator_type> A, MathVector<3> *pos, std::vector<vector_type> &testvectors)
	{
		testvectors.clear();
		for(size_t i=0; i<additionalTV.size(); i++)
		{
			THROW_IF_NOT_EQUAL(additionalTV[i].size(), A->num_rows());
			testvectors.push_back(additionalTV[i]);
		}
		UG_LOG("additional : " << testvectors.size() << " testvectors.\n");

		std::vector<vector_type> t;
		m_startTV->calc(A, pos, t);
		UG_LOG("AFT: " << testvectors.size() << " testvectors.\n");
		for(size_t i=0; i<t.size(); i++)
		{
			THROW_IF_NOT_EQUAL(t[i].size(), A->num_rows());
			testvectors.push_back(t[i]);
		}

	}

	void reset()
	{
		additionalTV.clear();
	}

	void add(const vector_type &vec)
	{
		additionalTV.push_back(vec);
	}

	virtual bool needs_positions() const { return m_startTV->needs_positions(); }
	virtual std::string name() const { return "AdaptiveFilteringTV"; }
	virtual std::string config_string() const
	{
		std::stringstream ss;
		ss << "AdaptiveFilteringTV. additional TV so far: " << additionalTV.size() << "\n";
		ss << " Start Testvectors: " << ConfigShift(m_startTV->config_string());
		return ss.str();
	}
};

#if 0
template <typename TAlgebra>
void FAMG<TAlgebra>::adaptive_filtering_multiplicative(SmartPtr<ILinearOperator<vector_type, vector_type> > linOp, vector_type &x, SmartPtr<IConvergenceCheck<vector_type> > spConvCheck,
		double expectedRate, SmartPtr<ITestvectorCalculator<TAlgebra> > spTV)
{
	UG_LOG("adaptive filtering\n");
	SmartPtr<MatrixOperator<matrix_type,vector_type,vector_type> > op = linOp.template cast_dynamic<MatrixOperator<matrix_type,vector_type,vector_type> >();


	SmartPtr<AdaptiveFilteringTestvectors<TAlgebra> > spAFT = make_sp(new AdaptiveFilteringTestvectors<TAlgebra>(spTV));
	set_testvectors(spAFT);

	// set solution to random
	SmartPtr<vector_type> spSolution = x.clone_without_values(); vector_type& solution = *spSolution;
	SmartPtr<vector_type> spError = x.clone_without_values(); vector_type& error = *spError;
	SmartPtr<vector_type> spB = x.clone_without_values(); vector_type& b = *spB;

	solution.set_random(-1, 1);

	// compute b = A*solution
	op->apply(b, solution);


	this->set_one_init(true);


	std::string writeMatrixPath = m_writeMatrixPath;

	SmartPtr<LinearIteratorProduct<vector_type, vector_type> > productIterator = new LinearIteratorProduct();
	std::vector<SmartPointer<FAMG<TAlgebra> > > vec;

	for(size_t i=0; i<10; i++)
	{
		UG_LOG("------------ TAKE " << i << " ------------ \n");
		UG_LOG(config_string() << "\n");
		// set x to zero
		x.set(0.0);



		m_writeMatrixPath = Stringify() << writeMatrixPath << "/adaptive_" << i << "/";
		CreateDirectory(m_writeMatrixPath);


		// try to solve A x = b

		vec.push_back(this->clone());
		productIterator->add_iterator(vec[i]);

		LinearSolver<vector_type> ls(productIterator, spConvCheck);
		ls.init(linOp);
		ls.apply(x, b);

		UG_LOG(vec[i]->level_information_str() << "\n");
		UG_LOG("rate = " << spConvCheck->rate() << "\n");
		UG_LOG("avg_rate = " << spConvCheck->avg_rate() << "\n");
		UG_LOG("step = " << spConvCheck->step() << "\n");
		UG_LOG("defect = " << spConvCheck->defect() << "\n");
		UG_LOG("reduction = " << spConvCheck->reduction() << "\n");

		if(spConvCheck->rate() < expectedRate)
		{
			UG_LOG("desired rate achieved.\n");
			break;
		}

		// error = solution - x
		VecScaleAdd(error, 1.0, solution, -1.0, x);
		double n = 1/error.norm();
		for(size_t i=0; i<error.size(); i++)
			error[i] *= n;

		if(m_writeMatrices)
		{
			ConnectionViewer::WriteVector(Stringify() << writeMatrixPath << "filter_" << i << ".vec",
						error, &m_amghelper.positions[0][0], m_dbgDimension);
		}

		// create new testvectors
		SmartPtr<AdaptiveFilteringTestvectors<TAlgebra> > spAFT = make_sp(new AdaptiveFilteringTestvectors<TAlgebra>(spTV));
		set_testvectors(spAFT);
		spAFT->add(error);
	}

	m_writeMatrixPath = writeMatrixPath;
}
#endif

template <typename TAlgebra>
void FAMG<TAlgebra>::adaptive_filtering(SmartPtr<ILinearOperator<vector_type, vector_type> > linOp, vector_type &x, SmartPtr<IConvergenceCheck<vector_type> > spConvCheck,
		double expectedRate, SmartPtr<ITestvectorCalculator<TAlgebra> > spTV)
{
	UG_LOG("adaptive filtering\n");
	SmartPtr<MatrixOperator<matrix_type,vector_type,vector_type> > op = linOp.template cast_dynamic<MatrixOperator<matrix_type,vector_type,vector_type> >();


	SmartPtr<AdaptiveFilteringTestvectors<TAlgebra> > spAFT = make_sp(new AdaptiveFilteringTestvectors<TAlgebra>(spTV));
	set_testvectors(spAFT);

	// set solution to random
	SmartPtr<vector_type> spSolution = x.clone_without_values(); vector_type& solution = *spSolution;
	SmartPtr<vector_type> spError = x.clone_without_values(); vector_type& error = *spError;
	SmartPtr<vector_type> spB = x.clone_without_values(); vector_type& b = *spB;

	solution.set_random(-1, 1);

	// compute b = A*solution
	op->apply(b, solution);

	std::string writeMatrixPath = m_writeMatrixPath;

	this->set_one_init(false);
	for(size_t i=0; i<10; i++)
	{
		UG_LOG("------------ TAKE " << i << " ------------ \n");
		UG_LOG(config_string() << "\n");
		// set x to zero
		x.set(0.0);

		m_writeMatrixPath = Stringify() << writeMatrixPath << "/adaptive_" << i << "/";
		CreateDirectory(m_writeMatrixPath);

		// try to solve A x = b


		UG_LOG("CREATING FAMG...\n");
		SmartPtr<FAMG<TAlgebra> > spFAMG(this, NULL);
		LinearSolver<vector_type> ls(spFAMG, spConvCheck);
		ls.init(linOp);
		ls.apply(x, b);
		UG_LOG(super::level_information_str() << "\n");

		UG_LOG("rate = " << spConvCheck->rate() << "\n");
		UG_LOG("avg_rate = " << spConvCheck->avg_rate() << "\n");
		UG_LOG("step = " << spConvCheck->step() << "\n");
		UG_LOG("defect = " << spConvCheck->defect() << "\n");
		UG_LOG("reduction = " << spConvCheck->reduction() << "\n");

		if(spConvCheck->rate() < expectedRate)
		{
			UG_LOG("desired rate achieved.\n");
			break;
		}

		// error = solution - x
		VecScaleAdd(error, 1.0, solution, -1.0, x);
		double n = 1/error.norm();
		for(size_t i=0; i<error.size(); i++)
			error[i] *= n;

		if(m_writeMatrices)
		{
			ConnectionViewer::WriteVector(Stringify() << writeMatrixPath << "filter_" << i << ".vec",
						error, &m_amghelper.positions[0][0], m_dbgDimension);
		}

		// add error to list of testvectors
		spAFT->add(error);
		set_testvectors(spAFT);

	}

	m_writeMatrixPath = writeMatrixPath;
}

}
#endif /* ADAPTIVE_FITLERING_H_ */
