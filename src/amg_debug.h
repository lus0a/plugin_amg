/**
 * \file amg_debug.h
 * \author Martin Rupp
 *
 * \date 16.06.10
 *
 * Goethe-Center for Scientific Computing 2009-2010.
 */



#ifndef __H__UG__LIB_DISC__AMG__AMG_DEBUG_H__
#define __H__UG__LIB_DISC__AMG__AMG_DEBUG_H__

#include <fstream>
#include "amg_debug_helper.h"
#include "lib_algebra/common/graph/graph.h"
#include "lib_algebra/lib_algebra.h"
#include "amg_profiling.h"
#include "postscript.h"

#ifdef UG_PARALLEL
#include "pcl/pcl_layout_tests.h"
#endif
#include "lib_algebra/common/connection_viewer_output.h"
#include "common/util/file_util.h"
#include "common/debug_print.h"

namespace ug {

/// \addtogroup plugins_amg
/// \{

extern DebugID DBG_AMG_COARSENING, DBG_RSAMG_PROLONGATION, DBG_AMG_AGGLOMERATION, DBG_AMG_COMMUNICATE_PROLONGATION;

// WriteToFile
//--------------------------------------------------
//! writes to a file in somewhat SparseMatrix-market format (for connection viewer)
template<typename T>
void AMGWriteToFile(const T &A, int fromlevel, int tolevel, std::string filename, const cAMG_helper &h)
{
	try{
	PROFILE_FUNC_GROUP("debug");
	if(h.has_positions() == false)
	{
		UG_LOG("AWriteToFile not possible: no positions available.")
		return;
	}

	std::fstream file(filename.c_str(), std::ios::out);

	int minlevel = std::min(fromlevel, tolevel);
	size_t N = std::max(A.num_rows(), A.num_cols());
	if(fromlevel != tolevel)
	{
		for(size_t i=0; i < A.num_rows(); i++)
			for(typename T::const_row_iterator conn = A.begin_row(i); conn != A.end_row(i); ++conn)
			{
				N = std::max(N, (size_t)h.GetOriginalIndex(tolevel, minlevel, i)+1);
				N = std::max(N, (size_t)h.GetOriginalIndex(fromlevel, minlevel, conn.index())+1);
			}
	}
	ConnectionViewer::WriteGridHeader(file, h.positions[minlevel], N, h.dimension);
	PROGRESS_START(prog, A.num_rows(), "AMGWriteToFile name = " << filename << " level = " << fromlevel << " size = " << A.num_rows());
	for(size_t i=0; i < A.num_rows(); i++)
	{
		PROGRESS_UPDATE(prog, i);
		for(typename T::const_row_iterator conn = A.begin_row(i); conn != A.end_row(i); ++conn)
			if(conn.value() != 0.0)
				file << h.GetOriginalIndex(tolevel, minlevel, i) << " " << h.GetOriginalIndex(fromlevel, minlevel, conn.index()) << " " << conn.value() << "\n";
	}
	}UG_CATCH_THROW_FUNC();
}

// WriteToFile
//--------------------------------------------------
//! writes to a file in somewhat SparseMatrix-market format (for connection viewer)
template<typename T>
void AMGWriteVecToFile(const T &v, int level, std::string filename, const cAMG_helper &h, const T *pCompare=NULL)
{
	PROFILE_FUNC_GROUP("debug");
	if(h.has_positions() == false)
	{
		UG_LOG("AWriteVecToFile not possible: no positions available.")
		return;
	}
	std::fstream file(filename.c_str(), std::ios::out);
	ConnectionViewer::WriteGridHeader(file, h.positions[level], h.dimension);
	if(pCompare == NULL)
		for(size_t i=0; i<v.size(); i++)
			file << i << " " << i << " " << v[i] << "\n";
	else
		for(size_t i=0; i<v.size(); i++)
			file << i << " " << i << " " << v[i] - (*pCompare)[i] << "\n";
}


// writeToFile
//--------------------------------------------------
//! writes to a file in somewhat SparseMatrix-market format (for connection viewer)
template<typename T>
void AMGWriteToFilePS(const T &A, int fromlevel, int tolevel, std::string filename, const cAMG_helper &h)
{
	PROFILE_FUNC_GROUP("debug");
	if(h.has_positions() == false)
	{
		UG_LOG("AWriteToFilePS not possible: no positions available.")
		return;
	}

	postscript ps;
	ps.create(filename.c_str());
	int minlevel = std::min(fromlevel, tolevel);
	for(size_t i=0; i < A.num_rows(); i++)
	{
		int from = h.GetOriginalIndex(tolevel, minlevel, i);
		ps.move_to(h.positions[minlevel][from].x(), h.positions[minlevel][from].y());
		ps.print_text( std::string("0") + ToString(i));

		for(typename T::const_row_iterator conn = A.begin_row(i); conn != A.end_row(i); ++conn)
		{
			if(conn.value() != 0.0)
			{
				if(conn.index() != i)
				{
					int to = h.GetOriginalIndex(fromlevel, minlevel, conn.index());
					ps.move_to(h.positions[minlevel][from].x(), h.positions[minlevel][from].y());
					ps.line_to(h.positions[minlevel][to].x(), h.positions[minlevel][to].y());

				}
			}
		}
	}

	//std::cout << std::endl;
}

// could be in cpp
inline void WriteAMGGraphToFile(const cgraph &G, std::string filename, const cAMG_helper &h, int level)
{
	try{
	// attention! this does only work serial.
	PROFILE_FUNC_GROUP("debug");
	if(h.has_positions() == false)
	{
		UG_LOG("WriteAMGGraphToFile not possible: no positions available.")
		return;
	}

	std::fstream file(filename.c_str(), std::ios::out);

	ConnectionViewer::WriteGridHeader(file, h.positions[level], h.dimension);
	PROGRESS_START(prog, G.size(), "WriteAMGGraphToFile name = " << filename << " level = " << level << " size = " << G.size());
	for(size_t i=0; i < G.size(); i++)
	{
		PROGRESS_UPDATE(prog, i);
		for(cgraph::const_row_iterator it = G.begin_row(i); it != G.end_row(i); ++it)
			file << i << " " << (*it) << " 1" << std::endl;
	}
	}UG_CATCH_THROW_FUNC();
}


#define PRINT_VECTOR_AMGDEBUG(dbgLevel, v, msg) {IF_DEBUG(LIB_ALG_AMG, 5){std::stringstream ss; ss << msg; PrintVector(v, ss.str());} }

#ifdef UG_PARALLEL
void MyPrintLayout(const pcl::ProcessCommunicator &pc,
				   pcl::InterfaceCommunicator<IndexLayout> &communicator,
				   const IndexLayout &layout1, const IndexLayout &layout2,
				   const char *name1, const char *name2);
#endif


/**
 * write the error made in prolongating the Linear elasticity testvectors
 */
template<typename TProlongationMatrix, typename vector_type>
void Write_LE_Error(const TProlongationMatrix &PoldIndices, const stdvector<MathVector<3> > &v, std::string writeMatrixPath, std::string name, int level)
{
	try{
	size_t N = PoldIndices.num_rows();
	vector_type tvLin(N);
	vector_type tvCon(N);

	vector_type Rres(N);
	vector_type Cres(N);
#ifdef UG_PARALLEL
	SetParallelVectorAsMatrix(Rres, PoldIndices, PST_CONSISTENT);
	SetParallelVectorAsMatrix(Cres, PoldIndices, PST_CONSISTENT);
#endif
	int dim=2;
	for(size_t i=0; i<N; i++)
	{
		Rres[i] = 0.0;
		Cres[i] = 0.0;
		if(PoldIndices.num_connections(i) <= 1)
			continue;
		const MathVector<3> &p0 = v[i];
		for(typename TProlongationMatrix::const_row_iterator it = PoldIndices.begin_row(i); it != PoldIndices.end_row(i); ++it)
		{
			size_t j = it.index();
			BlockRef(tvLin[j], 0) = p0.x() - v[j].x();
			BlockRef(tvLin[j], 1) = p0.y() - v[j].y();
			if(GetSize(tvLin[j]) == 3)
			{
				dim = 3;
				BlockRef(tvLin[j], 2) = p0.z() - v[j].z();
			}

			tvCon[j] = 1;
		}

		typename vector_type::value_type sumC = -1.0;
		typename vector_type::value_type sumR = 0.0;
		for(typename TProlongationMatrix::const_row_iterator it = PoldIndices.begin_row(i); it != PoldIndices.end_row(i); ++it)
		{
			sumC += it.value() * tvCon[it.index()];
			sumR += it.value() * tvLin[it.index()];
		}
		Rres[i] = sumR;
		Cres[i] = sumC;
	}
	std::string path=std::string("/level") + ToString(level) + "/";
	CreateDirectory(writeMatrixPath + path);
	ConnectionViewer::WriteVectorPar(writeMatrixPath + path + name + std::string("_Lin_err.vec"), Rres, &v[0], dim);
	ConnectionViewer::WriteVectorPar(writeMatrixPath + path + name + std::string("_Con_err.vec"), Cres, &v[0], dim);

	}UG_CATCH_THROW_FUNC();
}



///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// writevec
//--------------------------------
template<typename TAlgebra>
template<typename TVector>
bool AMGBase<TAlgebra>::writevec(std::string filename, const TVector &const_d, size_t level, const TVector *solution)
{
	size_t line=__LINE__;
	try{
	AMG_PROFILE_FUNC();
	UG_ASSERT(m_writeMatrices, "Huhh: writevec without m_writeMatrices?");

	line=__LINE__;

	bool bDrawConnections=false;

#ifdef UG_PARALLEL
	vector_type solution2;
	vector_type d;
	line=__LINE__;
	CloneVector(d, const_d);
	line=__LINE__;
	d = const_d;
	line=__LINE__;
	d.change_storage_type(PST_CONSISTENT);
	line=__LINE__;
	size_t pid = pcl::ProcRank();

	line=__LINE__;
	if(solution)
	{
		CloneVector(solution2, *solution);
		solution2 = *solution;
		solution2.change_storage_type(PST_CONSISTENT);
	}
#else
	vector_type solution2= *solution;
	const vector_type &d = const_d;
	size_t pid = 0;
#endif
	line=__LINE__;

	const char *ext = bDrawConnections ? ".mat" : ".vec";
	std::string name = get_output_directory(level) + std::string("AMG_L") + ToString(level) + "_" + filename + ext;

	line=__LINE__;
#ifdef UG_PARALLEL
	line=__LINE__;
	name = ConnectionViewer::GetParallelName(name, levels[level]->pA->layouts()->proc_comm(), true);
	line=__LINE__;
	/*if(levels[level]->bHasBeenMerged && m_agglomerateLevel != level)
		AMGWriteToFile(levels[level]->uncollectedA, level, level, name.c_str(), m_amghelper);
	else*/
#endif
	line=__LINE__;
	if(bDrawConnections)
	{
		line=__LINE__;
		AMGWriteToFile(*levels[level]->pA, level, level, name.c_str(), m_amghelper);
		std::fstream f(name.c_str(), std::ios::out | std::ios::app);
		//UG_LOG("writevec to " << name.c_str() << "\n");

		line=__LINE__;
		std::string name2 = get_output_directory(level) + std::string("AMG_L") + ToString(level) + "_" + filename + "_p"+ ToString(pid) + ".values";
	#ifdef UG_PARALLEL
		name2 = ConnectionViewer::GetParallelName(name2, levels[level]->pA->layouts()->proc_comm(), false);
	#endif
		line=__LINE__;
		f << "v " << name2 << "\n";

		line=__LINE__;
		std::fstream file(name2.c_str(), std::ios::out);
		if(solution)
			for(size_t i=0; i<d.size(); i++)
				file << i << " " << d[i] - solution2[i] << std::endl;
		else
			for(size_t i=0; i<d.size(); i++)
				file << i << " " << d[i] << std::endl;
	}
	else
	{
		line=__LINE__;
		if(solution)
			AMGWriteVecToFile(d, level, name, m_amghelper, &solution2);
		else
			AMGWriteVecToFile(d, level, name, m_amghelper);
	}
	line=__LINE__;
	}UG_CATCH_THROW(PRETTY_FUNCTION << "failed at line " << line );
	return true;
}


// end group plugins_amg
/// \}

}



#endif // __H__UG__LIB_DISC__AMG_SOLVER__AMG_DEBUG_H__
