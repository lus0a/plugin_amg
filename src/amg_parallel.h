/*
 * Copyright (c) 2013-2014:  G-CSC, Goethe University Frankfurt
 * Author: Arne Nägel
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

#ifndef AMG_PARALLEL_H_
#define AMG_PARALLEL_H_
namespace ug{

#ifdef UG_PARALLEL


// Assumption: Matrix has a unique distribution
template<typename matrix_type>
bool RedistributeMatrixRows(ParallelMatrix<matrix_type> &_mat, ParallelMatrix<matrix_type> &newMat)
{

	// for each master layout M_p linked to an upstream slave S_q:
	// a) send row entries to slave on upstream processor q>p
	// b) remove obsolete slave layout S_q



	// for each slave layout with a downstream master:
	// a) recv row entries from their master (downstream only!)
	// b) remove obsolete master layout M_q

}

template<typename matrix_type>
bool MakeFullRowsMatrix(const ParallelMatrix<matrix_type> &_mat, ParallelMatrix<matrix_type> &newMat,
		ParallelNodes &PN,
		stdvector<IndexLayout> &vMasterLayouts,
		stdvector<IndexLayout> &vSlaveLayouts,
		IParallelCoarsening *pParallelCoarsening=NULL)
{
	PROFILE_FUNC();

	// pcl does not use const much
	ParallelMatrix<matrix_type> &mat = const_cast<ParallelMatrix<matrix_type> &> (_mat);

	IndexLayout totalMasterLayout, totalSlaveLayout;

	GenerateOverlapClass<ParallelMatrix<matrix_type> > c(mat, newMat,
			totalMasterLayout, totalSlaveLayout, vMasterLayouts, vSlaveLayouts, PN);

	if (pParallelCoarsening==NULL) {
		c.m_overlapDepthMaster = 1;
		c.m_overlapDepthSlave = 0;
	} else{
		c.m_overlapDepthMaster = pParallelCoarsening->overlap_depth_master();
		c.m_overlapDepthSlave = pParallelCoarsening->overlap_depth_slave();
	}

	c.m_masterDirichletLast = false;
	c.m_slaveDirichletLast = false;
	bool b = c.calculate();

	//PRINTLAYOUT(mat.layouts()->proc_comm(), mat.communicator(), mat.master_layout(), mat.slave_layout());
	//PRINTLAYOUT(mat.layouts()->proc_comm(), mat.communicator(), vMasterLayouts[0], vSlaveLayouts[0]);

	//overlapSize = c.m_overlapSize;
	return b;
}
#endif

}
#endif /* AMG_PARALLEL_H_ */
