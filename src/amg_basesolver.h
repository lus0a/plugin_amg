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

namespace ug {

template<typename TAlgebra>
void AMGBase<TAlgebra>::create_direct_solver(size_t level)
{
#ifdef UG_PARALLEL
	m_baseLevel = level;
#endif
	try{
	AMG_PROFILE_BEGIN(amg_createDirectSolver);

	UG_ASSERT(block_traits< typename vector_type::value_type >::is_static, "dynamic not yet implemented");

	AMGLevelType &L = *levels[level];
	matrix_operator_type &A = *L.pA; (void) A;
	size_t static_nrUnknowns = block_traits< typename vector_type::value_type >::static_size;
	AMG_LOG(2, "Creating direct solver on level " << level << " (" << A.num_rows() << " nodes, total "
			<< A.num_rows()*static_nrUnknowns << " unknowns)" << std::endl << "Using Direct Solver on Matrix "
			<< A.num_rows()*static_nrUnknowns << "x" << A.num_rows()*static_nrUnknowns << ". ");
	(void) static_nrUnknowns;
	Stopwatch SW; SW.start();

#ifdef UG_PARALLEL
	const pcl::ProcessCommunicator &pc = A.layouts()->proc_comm(); //!
	if(m_basesolver->supports_parallel() == false && pc.size()>1)
	{

		L.bHasBeenMerged = true;
		L.bLevelHasMergers=true;
		// create basesolver
//		TemporaryDebugLevel tdl(LIB_ALG_AMG, 5);
//		TemporaryDebugLevel tdl2(LIB_ALG_MATRIX, 5);

		if(m_bChecks) { TESTLAYOUTS(A.layouts()); }

		CollectMatrixOnOneProc(A, *L.collectedA, L.agglomerateMasterLayout, agglomerateSlaveLayout);
		//UG_LOG("AGGLOMERATE MASTER LAYOUT\n\n" << L.agglomerateMasterLayout);
		//UG_LOG("AGGLOMERATE SLAVE LAYOUT\n\n" << agglomerateSlaveLayout);
		if(pcl::ProcRank() == pc.get_proc_id(0))
		{
			if(m_writeMatrices)
			{
				m_amghelper.receive_agglomerate_positions(level, A.layouts()->comm(), L.agglomerateMasterLayout, L.collectedA->num_rows());

				std::string name = get_output_directory(level) + "collectedA.mat";
				std::vector<MathVector<3> > vec = m_amghelper.positions[level];
				ConnectionViewer::WriteMatrix(name.c_str(), *L.collectedA, &vec[0], m_amghelper.dimension);
			}

			L.spAgglomeratedLayouts->proc_comm() = pc.create_sub_communicator(true);
			// setting empty layouts to collectedA
			L.collectedA->set_layouts(SmartPtr<AlgebraLayouts>(new AlgebraLayouts));
			m_basesolver->init(L.collectedA);

			L.corr.resize(L.collectedA->num_rows());
			L.collC.resize(L.collectedA->num_rows());
			L.collD.resize(L.collectedA->num_rows());
			SetParallelVectorAsMatrix(L.corr, *L.collectedA, PST_ADDITIVE);
			SetParallelVectorAsMatrix(L.collC, *L.collectedA, PST_ADDITIVE);
			SetParallelVectorAsMatrix(L.collD, *L.collectedA, PST_ADDITIVE);
		}
		else
		{
			if(m_writeMatrices)
				m_amghelper.send_agglomerate_positions(level, A.layouts()->comm(), agglomerateSlaveLayout);

			m_agglomerateLevel = level;
			L.spAgglomeratedLayouts->proc_comm() = pc.create_sub_communicator(false);
		}
		L.pAgglomeratedA = L.collectedA;
	}
	else
	{
		L.bLevelHasMergers=false;
		m_basesolver->init(L.pA);
		L.pAgglomeratedA = L.pA;
		L.corr.resize(L.pAgglomeratedA->num_rows());
		SetParallelVectorAsMatrix(L.corr, *L.pAgglomeratedA, PST_ADDITIVE);
	}
#else
	m_basesolver->init(L.pA);
	L.pAgglomeratedA = L.pA;
#endif

	m_dTimingCoarseSolverSetupMS = SW.ms();
	UG_DLOG(LIB_ALG_AMG, 1, "Coarse Solver Setup took " << m_dTimingCoarseSolverSetupMS << "ms." << std::endl);

	}UG_CATCH_THROW_FUNC();
}



template<typename TAlgebra>
bool AMGBase<TAlgebra>::solve_on_base(vector_type &c, vector_type &d, size_t level)
{
	AMG_PROFILE_FUNC();
	m_basesolver->apply_return_defect(c, d);
	return true;
}

}
