/*
 * debug_ids.cpp
 *
 *  Created on: 14.06.2013
 *      Author: mrupp
 */

#include "common/debug_id.h"

namespace ug{

DebugID DBG_FAMG("AMG.FAMG");
DebugID DBG_FAMG_CON_GRAPH("AMG.FAMG.CreateSymmConnectivityGraph");
DebugID DBG_FAMG_INTERPOLATION("AMG.FAMG.Interpolation");
DebugID DBG_FAMG_COLORING("AMG.FAMG.Coloring");
DebugID DBG_FAMG_OVERLAP("AMG.FAMG.Overlap");

double myTmpFAMGParam1=1.0;
int myTmpFAMGParam_LPE_N=5;
int myTmpFAMGParam_LE_N=5;
}






