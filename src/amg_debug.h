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
	if(m_writeMatrices==false)return true;
	size_t line=__LINE__;
	try{
	AMG_PROFILE_FUNC();


	line=__LINE__;

	bool bDrawConnections=false;

#ifdef UG_PARALLEL
	vector_type solution2;
	vector_type d;
	line=__LINE__;
	CloneVector(d, const_d);
	line=__LINE__;
	VecScaleAssign(d, 1.0, const_d);
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
	}UG_CATCH_THROW(PRETTY_FUNCTION << "failed at line " << line << " " << *const_d.layouts());
	return true;
}



template<typename TAlgebra>
std::string AMGBase<TAlgebra>::GetAMGFilenamePath(const matrix_type &m, int levelDir, std::string name, int level, std::string extension, bool bWriteHeader)
{
#ifdef UG_PARALLEL
	return ConnectionViewer::GetParallelName(get_output_directory(levelDir) + name + "_L" + ToString(level) + extension, m.layouts()->proc_comm(), bWriteHeader);
#else
	return get_output_directory(levelDir) + name + "_L" + ToString(level) + extension;
#endif
}

template<typename TAlgebra>
std::string AMGBase<TAlgebra>::GetAMGFilename(const matrix_type &m, std::string name, int level, std::string extension, bool bWriteHeader)
{
#ifdef UG_PARALLEL
	return ConnectionViewer::GetParallelName(name + "_L" + ToString(level) + extension, m.layouts()->proc_comm(), bWriteHeader);
#else
	return name + "_L" + ToString(level) + extension;
#endif
}

template<typename TAlgebra>
template<typename TNodeType>
void AMGBase<TAlgebra>::write_debug_matrix_markers
	(size_t level, const TNodeType &nodes)
{
	if(!m_writeMatrices) return;
	try{
	matrix_type &A = *levels[level]->pAgglomeratedA;;
	// todo: replace some day with something like nodes.get_mark_count(), nodes.mark_name(i), nodes.is_marked(i)


	std::string path = get_output_directory(level);

	AMG_PROFILE_FUNC();
	std::fstream ffine((path + GetAMGFilename(A, "_fine", level, ".marks")).c_str(), std::ios::out);
	ffine << "1 0 0 1 0\n";
	std::fstream ffine2((path + GetAMGFilename(A, "_aggfine", level, ".marks")).c_str(), std::ios::out);
	ffine2 << "1 0.2 1 1 0\n";
	std::fstream fcoarse((path + GetAMGFilename(A, "_coarse", level, ".marks")).c_str(), std::ios::out);
	fcoarse << "0 0 1 1 2\n";
	std::fstream fother((path + GetAMGFilename(A, "_other", level, ".marks")).c_str(), std::ios::out);
	fother << "1 1 0 1 0\n";
	std::fstream fdirichlet((path + GetAMGFilename(A, "_dirichlet", level, ".marks")).c_str(), std::ios::out);
	fdirichlet << "0 1 1 1 0\n";
	for(size_t i=0; i < A.num_rows(); i++)
	{
		//int o = m_amghelper.GetOriginalIndex(level, i);
		int o=i;
		if(nodes[i].is_aggressive_fine()) ffine2 << o << "\n";
		else if(nodes[i].is_fine_direct()) ffine << o << "\n";
		else if(nodes[i].is_coarse()) fcoarse << o << "\n";
		else if(nodes[i].is_dirichlet()) fdirichlet << o << "\n";
		else fother << o << "\n";
	}
	}UG_CATCH_THROW("AMGBase::write_debug_matrix_markers failed, level = " << level);
}




template<typename TAlgebra>
template<typename TMatrix>
void AMGBase<TAlgebra>::
	write_debug_matrix(TMatrix &mat, size_t fromlevel, size_t tolevel, const char *name, bool bAddMakers)
{
	if(!m_writeMatrices) return;
	try{
	AMG_PROFILE_FUNC();

	int level = std::min(fromlevel, tolevel);

	std::string filename = GetAMGFilenamePath(mat, level, name, fromlevel, ".mat", true);
	AMGWriteToFile(mat, fromlevel, tolevel, filename, m_amghelper);
	if(bAddMakers)
	{
		ConnectionViewer::AddMarkers(filename, GetAMGFilename(mat, "_fine", level, ".marks"));
		ConnectionViewer::AddMarkers(filename, GetAMGFilename(mat, "_coarse", level, ".marks"));
		ConnectionViewer::AddMarkers(filename, GetAMGFilename(mat, "_other", level, ".marks"));
		ConnectionViewer::AddMarkers(filename, GetAMGFilename(mat, "_dirichlet", level, ".marks"));
	}
	}UG_CATCH_THROW("AMGBase::write_debug_matrix failed, fromlevel = " << fromlevel << " tolevel = " << tolevel <<
			" name = " << name);
}

/*template<typename TAlgebra>
void AMGBase<TAlgebra>::save(size_t level, const vector_type &u, const vector_t &b)
{
	AMGLevelType &L = *levels[level];
	matrix_operator_type &A = *L.pA;

	std::fstream f;
	SerializeParallelData(f, A);
	SerializeUniquePart(f, A);
	SerializeUniquePart(f, u);
	SerializeUniquePart(f, b);
	bool b = m_amghelper.has_positions();
	if(b)
		Serialize(f, m_amghelper.positions[level]]);

}

template<typename TAlgebra>
void AMGBase<TAlgebra>::load(size_t level, const vector_type &u, const vector_t &b)
{
	AMGLevelType &L = *levels[level];
	matrix_operator_type &A = *L.pA;

	std::fstream f;
	SerializeParallelData(f, A);
	SerializeUniquePart(f, A);
	SerializeUniquePart(f, u);
	SerializeUniquePart(f, b);
	bool b = m_amghelper.has_positions();
	if(b)
		Serialize(f, m_amghelper.positions[level]]);

}
*/


template<typename TAlgebra>
void AMGBase<TAlgebra>::write_debug_matrices(matrix_type &AH, prolongation_matrix_type &R, const matrix_type &A,
		prolongation_matrix_type &P, size_t level)
{
	if(!m_writeMatrices) return;
	try{
	AMG_PROFILE_FUNC();
	AMG_LOG(1, "writing 4 matrices to " << m_writeMatrixPath << ": A ");

	//if(A.num_rows() < AMG_WRITE_MATRICES_MAX)
	{
		write_debug_matrix(A, level, level, "A", true);		AMG_LOG(1, "P ");
		write_debug_matrix(P, level+1, level, "P", true);	AMG_LOG(1, "R ");
		write_debug_matrix(R, level, level+1, "R", true);	AMG_LOG(1, "AH ");
	}

	AMGWriteToFile(AH, level+1, level+1,
			GetAMGFilenamePath(AH, level, "AH", level+1, ".mat", true), m_amghelper);

	AMG_LOG(1, " done.\n");
	}UG_CATCH_THROW("AMGBase::write_debug_matrices failed, level = " << level);
}


template<typename TAlgebra>
std::string
AMGBase<TAlgebra>::get_output_directory(int level)
{
	std::string path = Stringify() << m_writeMatrixPath << "/level" << level << "/";
	CreateDirectory(path);
	return path;
}

// end group plugins_amg
/// \}

}



#endif // __H__UG__LIB_DISC__AMG_SOLVER__AMG_DEBUG_H__
