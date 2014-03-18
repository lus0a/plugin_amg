/*
 * debug_ids.h
 *
 *  Created on: 19.06.2013
 *      Author: mrupp
 */

#ifndef DEBUG_IDS_H_
#define DEBUG_IDS_H_


#include "common/debug_id.h"

namespace ug{

extern DebugID DBG_FAMG_CON_GRAPH, DBG_FAMG_INTERPOLATION, DBG_FAMG, DBG_FAMG_COLORING, DBG_FAMG_OVERLAP;

// TemparyDebugLevel can change the debug level of a debug id an
// ensure that it is re-set to the original level at the end of the scope.
class TemporaryDebugLevel
{
	DebugID &did;
	int m_prevDebugLevel;
public:
	TemporaryDebugLevel(DebugID &debugID) : did(debugID)
	{
		m_prevDebugLevel = did.get_debug_level();
	}

	TemporaryDebugLevel(DebugID &debugID, int temporaryDebugLevel) : did(debugID)
	{
		m_prevDebugLevel = did.get_debug_level();
		did.set_debug_level(temporaryDebugLevel);
	}

	TemporaryDebugLevel(DebugID &debugID, DebugID &debugSrc, int diff=0) : did(debugID)
	{
		m_prevDebugLevel = did.get_debug_level();
		did.set_debug_level(debugSrc.get_debug_level()+diff);
	}

	~TemporaryDebugLevel()
	{
		did.set_debug_level(m_prevDebugLevel);
	}

	void set(int temporaryDebugLevel)
	{
		did.set_debug_level(temporaryDebugLevel);
	}
};

extern double myTmpFAMGParam1;
extern int myTmpFAMGParam_LPE_N;
extern int myTmpFAMGParam_LE_N;
}


#endif /* DEBUG_IDS_H_ */
