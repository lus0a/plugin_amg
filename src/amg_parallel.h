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
template<typename matrix_type>
bool MakeFullRowsMatrix(const ParallelMatrix<matrix_type> &_mat, ParallelMatrix<matrix_type> &newMat,
		ParallelNodes &PN, /* IParallelCoarsening *pParallelCoarsening,*/
		stdvector<IndexLayout> &vMasterLayouts, stdvector<IndexLayout> &vSlaveLayouts)
{
	PROFILE_FUNC();
	// pcl does not use const much
	//UG_ASSERT(overlap_depth > 0, "overlap_depth has to be > 0");
	ParallelMatrix<matrix_type> &mat = const_cast<ParallelMatrix<matrix_type> &> (_mat);

	IndexLayout totalMasterLayout, totalSlaveLayout;

	GenerateOverlapClass<ParallelMatrix<matrix_type> > c(mat, newMat,
			totalMasterLayout, totalSlaveLayout, vMasterLayouts, vSlaveLayouts, PN);
	//c.m_overlapDepthMaster = pParallelCoarsening->overlap_depth_master();
	//c.m_overlapDepthSlave = pParallelCoarsening->overlap_depth_slave();

	c.m_overlapDepthMaster = c.m_overlapDepthSlave = 1;
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
