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

/**
 * Goethe-Center for Scientific Computing 2011.
 */

#include <vector>
#include <map>

namespace ug
{

/// \addtogroup plugins_amg
///\{

// Agglomeration
//---------------------------------------------------------------------------
/**
 * serial load distribution algorithm
 *
 * @param sizes 		the sizes of the processing units
 * @param connections 	connections[i][j] = connection size between unit i and unit j
 * @param mergeWith		(returned) : if mergeWith[i].size() == 1, the processor to send matrix to (if == pcl::ProcRank, no mergin),
 * 						otherwise, list of nodes which are merged, with mergeWith[i][0] = i.
 * @param minimalSize	force all merged nodes to be greater than minimalSize
 * @param preferredSize
 *
 * Example:
 * minimalSize = 400
 * preferredSize = 4000
 * this setting ensures that if we have to agglomerate on one level, we try to do
 * as much agglomeration as possible, so that we have only 1 or 2 levels with many agglomerations and all other levels
 * are free of agglomeration (we want to prevent levels where only 2 processors agglomerate)
 *
 * mergeWith[0] = {0} : no merging, processor 1 stays
 * mergeWith[1] = {1} : no merging, processor 1 stays
 * mergeWith[2] = {5} : processor 2 send matrix to processor 5
 * mergeWith[3] = {5} : processor 3 send matrix to processor 5
 * mergeWith[4] = {5} : processor 4 send matrix to processor 5
 * mergeWith[5] = {5, 3, 4, 2}. : processor 5 receives matrix from processors 3, 4 and 2.
 *
 * @note this is a very easy serial agglomeration algorithm
 * @todo add agglomeration interface and other implementations like METIS
 */
void EasyAgglomeration(const std::vector<size_t> sizes,
		const std::vector<std::map<int, size_t> > connections,
		std::vector<std::vector<int> > &mergeWith,
		size_t minimalSize, size_t preferredSize);

// end group plugins_amg
/// \}

}
