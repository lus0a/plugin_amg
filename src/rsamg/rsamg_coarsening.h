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

#ifndef __H__LIB_ALGEBRA__AMG__AMG_COARSENING_H__
#define __H__LIB_ALGEBRA__AMG__AMG_COARSENING_H__

#include "lib_algebra/common/graph/graph.h"
#include "rsamg_nodeinfo.h"

namespace ug
{

/// \addtogroup plugins_amg_rsamg
/// \{


void CreateMeasureOfImportancePQ(const cgraph &strong, const cgraph &strongT, nodeinfo_pq_type &PQ, AMGNodes &nodes);

void CreateAggressiveCoarseningGraph(const cgraph &graph, cgraph &graph2, const AMGNodes &nodes,
		int nrOfPaths);


void CreateMeasureOfImportanceAggressiveCoarseningPQ(const cgraph &graphAC, nodeinfo_pq_type &PQ, AMGNodes &nodes);

int Coarsen(const cgraph &S, const cgraph &ST, nodeinfo_pq_type &PQ, AMGNodes &nodes, bool bUnsymmetric, bool bAggressiveCoarsening);

void PreventFFConnections(const cgraph &graphS, const cgraph &graphST, AMGNodes &nodes);


void RemoveUnassignedNeighbors(const cgraph &graph, nodeinfo_pq_type &PQ, AMGNodes &nodes, size_t i);
void MarkUnassignedNeighborsFine(const cgraph &S, const cgraph &ST, nodeinfo_pq_type &PQ, AMGNodes &nodes, size_t i, bool bMarkAsFineIndirect);
void ChangeRatingOfUnassignedNeighbors(const cgraph &graph, nodeinfo_pq_type &PQ, AMGNodes &nodes, size_t i, int change);

void PreventFFConnection(const cgraph &graphS, const cgraph &graphST, AMGNodes &nodes, size_t i, std::vector<bool> &marks, int &nrOfFFCoarseNodes);
int SetCoarse(int i, const cgraph &S, const cgraph &ST, nodeinfo_pq_type &PQ, AMGNodes &nodes, bool bUnsymmetric, bool bAggressiveCoarsening);
// end group plugins_amg_rsamg
/// \}

} // namespace ug

#endif // __H__LIB_ALGEBRA__AMG__AMG_COARSENING_H__
