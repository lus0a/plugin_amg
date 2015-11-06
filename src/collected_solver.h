/*
 * Copyright (c) 2012-2013:  G-CSC, Goethe University Frankfurt
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


#ifndef COLLECTED_SOLVER_H_
#define COLLECTED_SOLVER_H_

/// \addtogroup plugins_amg
/// \{

#ifndef UG_PARALLEL
#error "This only works with a UG_PARALLEL define."
#endif


template <typename TAlgebra>
class CollectSolver: public IMatrixOperatorInverse<	typename TAlgebra::matrix_type,
														typename TAlgebra::vector_type>
{
	public:
	// 	Algebra type
		typedef TAlgebra algebra_type;

	// 	Vector type
		typedef typename TAlgebra::vector_type vector_type;

	// 	Matrix type
		typedef typename TAlgebra::matrix_type matrix_type;

		typedef IMatrixOperatorInverse<	typename TAlgebra::matrix_type, typename TAlgebra::vector_type> moi_type;

	///	Base type
		typedef IMatrixOperatorInverse<matrix_type,vector_type> base_type;

	protected:
		using base_type::convergence_check;

	public:
		CollectSolver()
		{
			sName = "CollectedSolver(uninitialized)";
		};

		virtual const char* name() const {return name.c_str();}

		void set_convergence_check(IConvergenceCheck& convCheck)
		{
			UG_ASSERT(m_pInverse != NULL, "needs Inverse first");
			m_pInverse->set_convergence_check(convCheck);
		}
		IConvergenceCheck* convergence_check()
		{
			UG_ASSERT(m_pInverse != NULL, "needs Inverse first");
			return m_pInverse->convergence_check(convCheck);
		}

		void set_inverse(moi_type *inverse)
		{
			m_pInverse = inverse;
			sName = "CollectedSolver(";
			sName.append(inverse->name());
			sName.append(")");
		}

		moi_type *get_inverse()
		{
			return m_pInverse;
		}


		//	set operator L, that will be inverted
		virtual bool init(MatrixOperator<matrix_type, vector_type>& Op)
		{
		// 	remember operator
			m_pOperator = &Op;

		//	get matrix of Operator
			m_pMatrix = &m_pOperator->get_matrix();

		//	check that matrix exist
			if(m_pMatrix == NULL)
				{UG_LOG("ERROR in CollectSolver::init: No Matrix given,\n"); return false;}


			// create basesolver
			collect_matrix(*m_pMatrix, collectedMatrix, masterColl, slaveColl);

			m_SMO.setmatrix(&collectedMatrix);
			m_pInverse->init(m_SMO);
		}


		virtual bool apply_return_defect(vector_type& u, vector_type& f)
		{
			vector_type collC;
			vector_type collD;
			if(pcl::ProcRank() == 0)
			{
				size_t N = collectedBaseA.num_rows();
				collC.resize(N);
				collC.set(0.0);

				collD = d;
				collD.resize(N);
				for(size_t i=Ah.num_rows(); i<N; i++)
					collD[i] = 0.0;

			}


			// send d -> collD
			ComPol_VecAdd<vector_type > compolAdd(&collD, &d);
			pcl::InterfaceCommunicator<IndexLayout> &com = m_A[level]->communicator();
			com.send_data(slaveColl, compolAdd);
			com.receive_data(masterColl, compolAdd);
			com.communicate();

			if(pcl::ProcRank() == 0)
				m_pInverse->apply_return_defect(collC, collD);

			// send c -> collC
			ComPol_VecCopy<vector_type> compolCopy(&c, &collC);
			com.send_data(masterColl, compolCopy);
			com.receive_data(slaveColl, compolCopy);
			com.communicate();

			if(pcl::ProcRank() == 0)
				for(size_t i=0; i<Ah.num_rows(); i++)
					c[i] = collC[i];
			d.set(0.0);

			return true;
		}

	// 	Destructor
		virtual ~LUSolver() {};

	private:
		std::string sName;
		moi_type *m_pInverse;
};

// end group plugins_amg
/// \}

#endif /* COLLECTED_SOLVER_H_ */
