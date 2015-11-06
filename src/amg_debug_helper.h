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


#ifndef __H__UG__LIB_DISC__AMG_SOLVER__AMG_DEBUG_HELPER_H__
#define __H__UG__LIB_DISC__AMG_SOLVER__AMG_DEBUG_HELPER_H__

#include "common/common.h" // MathVector<3>
#include "common/math/ugmath.h"
#include "lib_algebra/common/stl_debug.h" // stdvector

#ifdef UG_PARALLEL
#include "pcl/pcl.h"
#include "lib_algebra/parallelization/parallel_index_layout.h"
#include "lib_algebra/parallelization/communication_policies.h"
#endif

#include "amg_profiling.h"

namespace ug{

/// \addtogroup plugins_amg
/// \{

/**
 * \brief helper struct for debug output of matrices in AMG
 * This class knows the position to each node on each level of
 * the algebraic grid.
 * \sa WriteAMGGraphToFile, AMGWriteToFile, AMGWriteToFilePS
 */
struct cAMG_helper
{

	typedef MathVector<3> position_type;

	cAMG_helper() : dimension(0), parentIndex(NULL) {	}

	//! returns to the node i on level 'level' the corresponding node on the level 0
	int GetOriginalIndex(int level, int i) const
	{
		return GetOriginalIndex(level, 0, i);
	}

	//! returns to the node i on level 'level' the corresponding node on the level 0
	int GetOriginalIndex(int level, int onlevel, int i) const
	{
		while(level != onlevel)
		{
			if(i >= (int)(*parentIndex)[level].size())
				return 0;
			i = (*parentIndex)[level--][i];
		}
		return i;
	}

	//! returns to the node i on level 'level' the position
	position_type GetPosForIndexAtLevel(int level, int i) const
	{
		return positions[level][i];
	}

#ifdef UG_PARALLEL
	void update_overlap_positions(int level,
			pcl::InterfaceCommunicator<IndexLayout> &communicator, IndexLayout &overlapMaster, IndexLayout &overlapSlave,
			size_t newSize)
	{
		//AMG_PROFILE_FUNC();
		std::vector<MathVector<3> > &vec2 = positions[level];
		vec2.resize(newSize);

		ComPol_VecCopy<std::vector<MathVector<3> > >	copyPol(&vec2);
		communicator.send_data(overlapMaster, copyPol);
		communicator.receive_data(overlapSlave, copyPol);
		communicator.communicate();
	}

	void receive_agglomerate_positions(int level,
			pcl::InterfaceCommunicator<IndexLayout> &communicator,
			IndexLayout &agglomerateMaster, size_t newSize)
	{
		std::vector<MathVector<3> > &vec2 = positions[level];
		vec2.resize(newSize);

		ComPol_VecCopy<std::vector<MathVector<3> > >	copyPol(&vec2);
		communicator.receive_data(agglomerateMaster, copyPol);
		communicator.communicate();
	}

	void send_agglomerate_positions(int level,
			pcl::InterfaceCommunicator<IndexLayout> &communicator,
			IndexLayout &agglomerateSlave)
	{
		std::vector<MathVector<3> > &vec2 = positions[level];

		ComPol_VecCopy<std::vector<MathVector<3> > >	copyPol(&vec2);
		communicator.send_data(agglomerateSlave, copyPol);
		communicator.communicate();
	}
#endif

	void make_coarse_level(size_t level, const stdvector<int> &parentIndex)
	{
		//AMG_PROFILE_FUNC();
		positions.resize(level+1);
		positions[level].resize(parentIndex.size());
		for(size_t i=0; i < parentIndex.size(); i++)
			positions[level][i] = positions[level-1][parentIndex[i]];
	}
	bool has_positions() const
	{
		return dimension != 0;
	}

	position_type *get_positions(size_t level)
	{
		return &positions[level][0];
	}

	int get_dimension()
	{
		return dimension;
	}

public:

	stdvector<stdvector<MathVector<3> > > positions; ///< positions on the AMG grid 0
	int dimension;					///< dimension (2 or 3)
	stdvector< stdvector<int> > *parentIndex;				///< (*parentIndex)[L][i] is the index of i on level L-1
};

// end group plugins_amg
/// \}

} // namespace ug
#endif // __H__UG__LIB_DISC__AMG_SOLVER__AMG_DEBUG_HELPER_H__
