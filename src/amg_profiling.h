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


#ifndef __H__LIB_ALGEBRA__AMG_PROFILING_H_
#define __H__LIB_ALGEBRA__AMG_PROFILING_H_

/// \addtogroup plugins_amg
/// \{

#define PROFILE_AMG
#ifdef PROFILE_AMG
	#include  "common/profiler/profiler.h"
	#define AMG_PROFILE_FUNC()			PROFILE_FUNC_GROUP("algebra AMG")
	#define AMG_PROFILE_BEGIN(name)		PROFILE_BEGIN_GROUP(name, "algebra AMG")
	#define AMG_PROFILE_NEXT(name)		PROFILE_END(); PROFILE_BEGIN_GROUP(name, "algebra AMG")
	#define AMG_PROFILE_END()			PROFILE_END()
	#define AMG_BPROFILE_BEGIN(name, string) PROFILE_BEGIN_GROUP(name, "algebra AMG"); Stopwatch SW; const char *profileString = string; SW.start();
	#define AMG_BPROFILE_END() PROFILE_END(); UG_DLOG(LIB_ALG_AMG, 1, string << " took " << SW.ms() << " ms\n");
	#define AMG_BPROFILE_NEXT(name, string) AMG_BPROFILE_END(); PROFILE_BEGIN_GROUP(name, "algebra AMG"); SW.start();

#else
	#define AMG_PROFILE_FUNC()
	#define AMG_PROFILE_BEGIN(name)
	#define AMG_PROFILE_NEXT(name)
	#define AMG_PROFILE_END()
	#define AMG_BPROFILE_BEGIN(name, string)
	#define AMG_BPROFILE_NEXT(name, string)
	#define AMG_BPROFILE_END()
#endif


#define AMG_LOG(lvl, msg) { if(m_verbose >= lvl) UG_LOG(msg); }
#define AMG_LOGN(lvl, msg) AMG_LOG(lvl, msg << "\n");


#define AMG_START_TIMING(msg) Stopwatch AMG_TIMING_SW; if(bTiming) { m_ssTiming << msg; SW.start(); }
#define AMG_END_TIMING(msg) if(bTiming) { m_ssTiming << " took " << AMG_TIMING_SW.ms() << " ms\n"; }

// end group plugins_amg
/// \}

#endif /* __H__LIB_ALGEBRA__AMG_PROFILING_H_ */
