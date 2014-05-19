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
