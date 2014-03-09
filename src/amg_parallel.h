/*
 * amg_parallel.h
 *
 *  Created on: 23.11.2013
 *      Author: anaegel
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
