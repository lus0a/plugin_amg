/*
 * rsamg_blockvalues.h
 *
 *  Created on: Jun 12, 2012
 *      Author: rupp4
 */

#ifndef __H__UG__LIB_ALGEBRA__RSAMG_BLOCKVALUES_H_
#define __H__UG__LIB_ALGEBRA__RSAMG_BLOCKVALUES_H_
#include "lib_algebra/small_algebra/small_algebra.h"
#include "lib_algebra/small_algebra/additional_math.h"

/// \addtogroup plugins_amg_rsamg
/// \{

inline double amg_diag_value(const double &d) { return d; }
inline double amg_offdiag_value(const double &d) { return d; }

template<typename T> inline double amg_diag_value(const T &d) { return BlockNorm(d); }
template<typename T> inline double amg_offdiag_value(const T &d) { return -BlockNorm(d); }

template<typename T> inline double amg_trace(const T &mij) {return ug::matTrace(mij);}

// end group plugins_amg_rsamg
/// \}

#endif /* RSAMG_BLOCKVALUES_H_ */
