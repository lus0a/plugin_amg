/*
 * Copyright (c) 2012-2014:  G-CSC, Goethe University Frankfurt
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

#ifndef UG_PARALLEL
#error "This only works with a UG_PARALLEL define."
#endif


#ifndef RSAMG_PARALLEL_COARSENING_H_
#define RSAMG_PARALLEL_COARSENING_H_

#include "lib_algebra/common/graph/graph.h"
#include "rsamg_nodeinfo.h"

#include "pcl/pcl.h"
#include "lib_algebra/parallelization/parallel_index_layout.h"


namespace ug
{

/// \addtogroup plugins_amg_rsamg
/// \{

class IParallelCoarsening
{
public:
	virtual ~IParallelCoarsening() {}
	virtual void coarsen(ParallelNodes &PN, stdvector<IndexLayout> vMasterLayouts, stdvector<IndexLayout> vSlaveLayouts,
				const cgraph &graphS, const cgraph &graphST, nodeinfo_pq_type &PQ, AMGNodes &nodes, bool bUnsymmetric, bool bPreventFF) = 0;
	virtual int overlap_depth_master() = 0;
	virtual int overlap_depth_slave() = 0;
	virtual const char *config_string() = 0;
};


IParallelCoarsening *GetFullSubdomainBlockingCoarsening();
IParallelCoarsening *GetColorCoarsening();
IParallelCoarsening *GetRS3Coarsening();
IParallelCoarsening *GetCLJPCoarsening();
IParallelCoarsening *GetFalgoutCoarsening();
IParallelCoarsening *GetMinimumSubdomainBlockingCoarsening();
IParallelCoarsening *GetCoarseGridClassificationCoarsening();
IParallelCoarsening *GetSimpleParallelCoarsening();

// end group plugins_amg_rsamg
/// \}

}
#endif /* RSAMG_PARALLEL_COARSENING_H_ */
