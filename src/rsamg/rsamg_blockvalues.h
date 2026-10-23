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


#ifndef __H__UG__LIB_ALGEBRA__RSAMG_BLOCKVALUES_H_
#define __H__UG__LIB_ALGEBRA__RSAMG_BLOCKVALUES_H_
#include "lib_algebra/small_algebra/small_algebra.h"
#include "lib_algebra/small_algebra/additional_math.h"

/// \addtogroup plugins_amg_rsamg
/// \{

inline double amg_diag_value(const double &d) { return d; }
inline double amg_offdiag_value(const double &d) { return d; }
inline double amg_signed_value(const double &d) { return d; }

template<typename T> inline double amg_diag_value(const T &d) { return BlockNorm(d); }
template<typename T> inline double amg_offdiag_value(const T &d) { return -BlockNorm(d); }

// BlockNorm is appropriate for measuring connection strength, but it destroys
// the sign of interpolation weights.  Use the trace to recover the sign while
// retaining the norm as a robust magnitude for non-scalar blocks.
template<typename T> inline double amg_signed_value(const T &d)
{
	const double trace = ug::matTrace(d);
	return trace < 0.0 ? -BlockNorm(d) : BlockNorm(d);
}

template<typename T> inline double amg_trace(const T &mij) {return ug::matTrace(mij);}
template<typename T> inline double amg_diag(const T &mij) {return ug::matDiagMax(mij);}

// end group plugins_amg_rsamg
/// \}

#endif /* RSAMG_BLOCKVALUES_H_ */
