/*
 * Copyright (c) 2013-2014:  G-CSC, Goethe University Frankfurt
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

#ifndef LOCAL_EIGENVALUE_TEST_COMPUTATION_H_
#define LOCAL_EIGENVALUE_TEST_COMPUTATION_H_


#if 0
void calc_special(size_t i)
{

}
#else
template<typename T>
void print_matlab(DenseMatrix<T> &m, std::string name)
{

	UG_LOG("\n" << name << reset_floats << "= zeros(" << m.num_rows() << ", " << m.num_cols() << ")\n");
	for(size_t r=0; r<m.num_rows(); r++)
	{
		for(size_t c=0; c<m.num_cols(); c++)
		{
			if(m(r, c) != 0)
			{UG_LOG(name << "[" << r+1 << ", " << c+1 << "] = " << m(r,c) << "\n");}
		}
	}

	UG_LOG("\n");
}

void print_matlab_positions(std::vector<size_t> indices)
{
	size_t level = m_amghelper.positions.size()-1;
	UG_LOG("pos = Array(Pos3d, " << indices.size() << ")\n");
	for(size_t i=0; i<indices.size(); i++)
	{
		UG_LOG("pos[" << i+1 << "] = Pos3d(" << m_amghelper.positions[level][indices[i]][0] << ", "
				<< m_amghelper.positions[level][indices[i]][1] << ", " << m_amghelper.positions[level][indices[i]][2] << ")\n");
	}
	UG_LOG("\n")
}

template<typename T>
inline void ToDouble(DenseMatrix<VariableArray2<double> > &A2, const T &A)
{
	size_t blockSize = GetRows(A(0,0));
	A2.resize(blockSize*A.num_rows(), blockSize*A.num_cols());

	for(size_t i=0; i<A.num_rows(); i++)
	{
		for(size_t j=0; j<A.num_cols(); j++)
			for(size_t r=0; r<blockSize; r++)
				for(size_t c=0; c<blockSize; c++)
					A2(i*blockSize+r, j*blockSize+c) = BlockRef(A(i, j), r, c);
	}
}
/*
template<typename A_type, typename B_type, typename vector_type>
int GeneralizedEigenvalueProblemSys(A_type &A, std::vector<vector_type> &X,
		B_type &lambda, A_type &B, bool bSortEigenvalues=false)
{
	size_t blockSize = GetRows(A(0,0));
	DenseMatrix<VariableArray2<double> > A2, B2, X2;
	size_t N = A.num_rows();
	size_t Nb = N*blockSize;
	std::vector<double> lambda;
	ToDouble(A2, A);
	ToDouble(B2, B);
	lambda.resize(Nb);
	X2.resize(Nb, Nb);

	GeneralizedEigenvalueProblem(A2, X2, lambda, B2, bSortEigenvalues);
	X.resize(Nb);
	for(size_t i=0; i<Nb; i++)
	{
		X[i].resize(N);
		for(size_t j=0; j<N; j++)
			for(size_t k=0; k<blockSize; k++)
				BlockRef(X[i][j], k) = X2(j*blockSize+k, i);
	}
}


void calc_special22(size_t i)
{
	if(!is_watch_node(i))
		return;
	std::vector<std::vector<size_t> > indicesSysH;
	GetNeighborhoodHierachy(A_OL2, i, 3, indicesH);

	std::vector<size_t> indices;
	for(size_t k=0; k <= 3; k++)
		indices.insert(indices.end(), indicesH[k].begin(), indicesH[k].end());

	size_t indices1=indicesH[0].size();
	size_t indices2=indicesH[1].size()+indices1;
	size_t indices3=indicesH[2].size()+indices2;
	size_t indices4=indices.size();
	size_t N=indices4;

	DenseMatrix<VariableArray2<smallmat_type> > Aloc;

	Aloc.resize(indices.size(), indices.size());

	GetLocalMatrix(A_OL2, Aloc, &indices[0], &indices[0]);

	DenseMatrix<VariableArray2<smallmat_type> > Dinv;
	Dinv.resize(N, N);
	for(size_t r=0; r<N; r++)
	{
		Dinv(r,r) = S(r,r);
		Invert(Dinv(r,r));
	}

	DenseMatrix<VariableArray2<smallmat_type> > diagm;
	diagm.resize(N, N);
	for(size_t r=0; r<N; r++)
		diagm(r,r)=1.0;

	DenseMatrix<VariableArray2<smallmat_type> > S = diagm - Dinv*A;
	DenseMatrix<VariableArray2<smallmat_type> > SF = diagm;
	for(size_t c=0; c<N; c++)
		SF(0, c) = S(0, c);

	S = SF*S;
	writeSubMatrix("/Users/mrupp/matrices/S.mat", S, indices, indices4);

	DenseVector<VariableArray1<smallvec_type> > q;
	DenseVector<VariableArray1<smallvec_type> > Sq;

	DenseMatrix<VariableArray2<double> > ev;
	std::vector<double> lambda;
	size_t blockSize = GetRows(Aloc(0,0));
	size_t Nb = N*blockSize;

	DenseMatrix<VariableArray2<double> > A2, Q2;
	A2.resize(Nb, Nb);
	ev.reize(Nb, Nb);
	lambda.resize(Nb);

	ToDouble(A2, Aloc);


	for (size_t p0 = indices1; p0 < indices2; p0++)
	{
		size_t gp0 = indices[p0];
		for (size_t p1 = p0+1; p1 < indices2; p1++)
		{

			double w1=0.5;
			double w2=0.5;
			double eps=0.5;

			//while(eps > 0.001)
			{
				q*=0.0;
				q[0] = 1.0;
				q[p0] = -w1;
				q[p1] = -w2;

				Sq = transpose(S)*q;
				DenseMatrix<VariableArray2<smallmat_type> > Q;
				for(size_t r=0; r<N; r++)
					for(size_t c=0; c<N; c++)
						Q(r, c) = Sq[r]*Sq[c];

				ToDouble(Q2, Q);


				// bei blockalgorithmen entstehen hier N*blockSize eigenwerte
				GeneralizedEigenvalueProblem(Q2, ev, lambda, A2, true);

				writeSubVec(std::string("/Users/mrupp/matrices/ev_")+ToString(gp0)+"_"+ToString(gp0)+"_"+ToString(lambda[N-1])+"_.vec",
						A, ev, N-1, indices, indices4);
			}

		}
	}




	// select q
	// calc Sq=S*q
	// calc Q = Sq*Sq^T
	// (lambda, v)=eig(Q, A)
	// iMax = argmax(lambda)
	// emax = v[:, iMax]
	// lmax = lambda[iMax]
	// return lmax



}*/


void calc_special(size_t i)
{
	if(!is_watch_node(i))
	return;

	//calc_special_local_eigenvalues(i);
	calc_special_local_projected_eigenvalues(i);
	calc_special_local_eigenvalues(i);
}

template<typename T, typename T2>
void GetLocalVector(const T &tsrc, T2 &tdest, std::vector<size_t> &indices)
{
	tdest.resize(indices.size());
	for(size_t i=0; i<indices.size(); i++)
		tdest[i] = tsrc[indices[i]];
}

void calc_special_local_projected_eigenvalues(size_t i)
{
	PROFILE_FUNC_GROUP("debug");
	UG_LOG("\n--------------------------------------------------\n");
	std::vector<std::vector<size_t> > indicesH;
	GetNeighborhoodHierachy(A_OL2, i, myTmpFAMGParam_LPE_N , indicesH);

	std::vector<size_t> indices;
	PROFILE_BEGIN(A);

	for(size_t k=0; k < indicesH.size(); k++)
		indices.insert(indices.end(), indicesH[k].begin(), indicesH[k].end());

	size_t indices1=indicesH[0].size();
	size_t indices2=indicesH[1].size()+indices1;
	size_t indices3=indicesH[2].size()+indices2;
	size_t indices4=indices.size();
	//indices3 = indices4;

	DenseMatrix<VariableArray2<smallmat_type> > Aloc;

	Aloc.resize(indices.size(), indices.size());

	GetLocalMatrix(A_OL2, Aloc, &indices[0], &indices[0]);

	DenseMatrix<VariableArray2<smallmat_type> > Bloc;
	Bloc.resize(indices4, indices4);
	for(size_t j=0; j<indices4; j++)
		Bloc(j,j) = Aloc(j,j);

	stdvector<DenseVector<VariableArray1<smallvec_type> > > tv;
	tv.resize(m_testvectors.size());
	DenseVector<VariableArray1<smallvec_type> > tmp;

	for(size_t j=0; j<m_testvectors.size(); j++)
		GetLocalVector(m_testvectors[j], tv[j], indices);

	DenseMatrix<VariableArray2<double> > ewA, ewB;
	ewA.resize(tv.size(), tv.size());
	ewB.resize(tv.size(), tv.size());

	for(size_t j1=0; j1<tv.size(); j1++)
	{
		tmp = Aloc*tv[j1];
		for(size_t k=indices3; k<indices4; k++) tmp[k] *= myTmpFAMGParam1;
		for(size_t j2=0; j2<tv.size(); j2++)
			ewA(j1, j2) = VecProd(tmp, tv[j2]);
	}
	for(size_t j1=0; j1<tv.size(); j1++)
	{
		tmp = Bloc*tv[j1];
//		for(size_t k=indices3; k<indices4; k++) tmp[k] = 0.0;
		for(size_t j2=0; j2<tv.size(); j2++)
			ewB(j1, j2) = VecProd(tmp, tv[j2]);
	}
	ewA.maple_print("A");
	ewB.maple_print("B");


	DenseMatrix<VariableArray2<double> > ev;
	std::vector<std::complex<double> > lambda;
	ev.resize(ewA.num_rows(), ewA.num_rows());
	lambda.resize(ewA.num_rows());
	GeneralizedEigenvalueProblemComplex(ewA, ev, lambda, ewB, true);

	for(size_t k=0; k<ev.num_rows(); k++)
	{
		tmp = 0.0;
		for(size_t c = 0; c<ev.num_rows(); c++)
			tmp += ev(c, k) * tv[c];

		writeSubVec2(Stringify() << m_matrixWritePath << "/PEV2_" << i << "_ev_" << k << "__" << lambda[k] << ".vec", Aloc, tmp, indices, indices2);
		writeSubVec2(Stringify() << m_matrixWritePath << "/PEV3_" << i << "_ev_" << k << "__" << lambda[k] << ".vec", Aloc, tmp, indices, indices3);
		writeSubVec2(Stringify() << m_matrixWritePath << "/PEV_" << i << "_ev_" << k << "__" << lambda[k] << ".vec", Aloc, tmp, indices, indices4);
	}
}

void writeSubVec2(std::string filename, const DenseMatrix<VariableArray2<smallmat_type> > &mat,
		DenseVector<VariableArray1<smallvec_type> > &v, const std::vector<size_t> &indices, size_t N)
{
	if(m_amghelper.has_positions() == false) return;
	size_t level = m_amghelper.positions.size()-1;

	std::fstream file(filename.c_str(), std::ios::out);

	ConnectionViewer::WriteGridHeader(file, m_amghelper.positions[level], m_amghelper.dimension);

	size_t blockSize = GetRows(mat(0,0));

//		for(size_t i=0; i < m_amghelper.positions[level].size(); i++)
//			file << i << " " << i << " 0.0\n";

	for(size_t i=0; i < N; i++)
	{
		if(blockSize == 1)
			file << indices[i] << " " << indices[i] << " " << v[i] << "\n";
		else
		{
			file << indices[i] << " " << indices[i] << " [ ";
			for(size_t j=0; j<blockSize; j++)
				file << BlockRef(v[i], j) << " ";
			file << "]\n";
		}

	}

	/*for(size_t i=0; i < N; i++)
		for(size_t j=0; j < N; j++)
			if(i != j && BlockNorm(mat(i,j)) > 1e-9)
			{
				file << indices[i] << " " << indices[j] << " " << mat(i,j) << "\n";
			}*/
}

void calc_special_local_eigenvalues(size_t i)
{
	PROFILE_FUNC_GROUP("debug");
	std::vector<std::vector<size_t> > indicesH;
	GetNeighborhoodHierachy(A_OL2, i, myTmpFAMGParam_LE_N, indicesH);

	std::vector<size_t> indices;

	for(size_t k=0; k < indicesH.size(); k++)
		indices.insert(indices.end(), indicesH[k].begin(), indicesH[k].end());

	size_t indices1=indicesH[0].size();
	size_t indices2=indicesH[1].size()+indices1;
	size_t indices3=indicesH[2].size()+indices2;
	size_t indices4=indices.size();
	indices3 = indices4;

	DenseMatrix<VariableArray2<smallmat_type> > Aloc;

	Aloc.resize(indices.size(), indices.size());

	GetLocalMatrix(A_OL2, Aloc, &indices[0], &indices[0]);

	UG_LOG("\nindices = [" << indices1 << ", " << indices2 << ", " << indices3 << ", " << indices4 << "]\n");
	print_matlab(Aloc, "A");
	print_matlab_positions(indices);

	// N2-N3 connections
	DenseMatrix<VariableArray2<smallmat_type> > Bloc;
	Bloc.resize(indices4, indices4);
	for(size_t j=0; j<indices4; j++)
		Bloc(j,j) = Aloc(j,j);
//#define FAMG_EV_ATDA
#define FAMG_EV_DIRCHLET
//#define FAMG_EV_DIRICHLET_NEUMANN
#ifdef FAMG_EV_DIRCHLET
//		for(size_t i2=indices2; i2<indices3; i2++)
//			for(size_t j3=indices3; j3<indices4; j3++)
//				Aloc(i2, i2) += Aloc(i2, j3);
#endif

#ifdef FAMG_EV_DIRICHLET_NEUMANN
	for(size_t i2=indices2; i2<indices3; i2++)
		for(size_t j3=indices3; j3<indices4; j3++)
		{
			if(Aloc(i2, j3) == 0.0) continue;
			Aloc(i2, i2) += Aloc(i2, j3);
			//Bloc(i2, i2) += Bloc(j3, j3);
		}
#endif

	//		for(size_t j=indices2; j<indices3; j++)	Bloc(j,j) *= 1000.0;
#ifdef FAMG_EV_ATDA
	DenseMatrix<VariableArray2<smallmat_type> > Aloc2;
	Aloc2.resize(indices3, indices3);
	// A' = A^T D A   mit D = diag(A) auf N1 !
	// A_{ij} = sum A_{ki} D_{kk} A_{kj}
	for(size_t i=0; i<indices3; i++)
		for(size_t j=0; j<indices3; j++)
		{
			Aloc2(i, j) = 0.0;
			for(size_t k=0; k<indices2; k++)
				Aloc2(i, j) += Aloc(k, i) / Aloc(k,k) * Aloc(k, j);
		}
	Aloc = Aloc2;

#endif


	size_t blockSize = GetRows(Aloc(0,0));

	DenseMatrix<VariableArray2<double> > A2;
	A2.resize(indices3*blockSize, indices3*blockSize);

	DenseMatrix<VariableArray2<double> > B2;
	B2.resize(indices3*blockSize, indices3*blockSize);

	for(size_t i=0; i<indices3; i++)
	{
		for(size_t j=0; j<indices3; j++)
			for(size_t r=0; r<blockSize; r++)
				for(size_t c=0; c<blockSize; c++)
					A2(i*blockSize+r, j*blockSize+c) = BlockRef(Aloc(i, j), r, c);
		for(size_t r=0; r<blockSize; r++)
			B2(i*blockSize+r, i*blockSize+r) = BlockRef(Bloc(i, i), r, r);
	}


	PrintVector(indicesH[0], "N0");
	PrintVector(indicesH[1], "N1");
	PrintVector(indicesH[2], "N2");
	A2.maple_print("A2");

	writeSubMatrix(m_matrixWritePath+std::string("/Aloc.mat"), Aloc, indices, indices3);
	writeSubMatrix(m_matrixWritePath+std::string("/Bloc.mat"), Bloc, indices, indices3);

	DenseMatrix<VariableArray2<double> > ev;
	std::vector<std::complex<double> > lambda;
	ev.resize(A2.num_rows(), A2.num_rows());
	lambda.resize(A2.num_rows());
	GeneralizedEigenvalueProblemComplex(A2, ev, lambda, B2, true);

	for(size_t k=0; k<ev.num_rows(); k++)
	{
		writeSubVec(Stringify() << m_matrixWritePath << "/LEV_" << i << "_ev_" << k << "__" << lambda[k] << ".vec", Aloc, ev, k, indices, indices3);
	}

}
void writeSubMatrix(std::string filename, const DenseMatrix<VariableArray2<smallmat_type> > &mat, const std::vector<size_t> &indices, size_t N)
{
	if(m_amghelper.has_positions() == false) return;
	size_t level = m_amghelper.positions.size()-1;

	std::fstream file(filename.c_str(), std::ios::out);

	ConnectionViewer::WriteGridHeader(file, m_amghelper.positions[level], m_amghelper.dimension);
	for(size_t r=0; r < N; r++)
		for(size_t c=0; c < N; c++)
			if(mat(r, c) != 0.0)
				file << indices[r] << " " << indices[c] << " " << mat(r, c) << "\n";
}

void writeSubVec(std::string filename, const DenseMatrix<VariableArray2<smallmat_type> > &mat, const DenseMatrix<VariableArray2<double> > &v,
		size_t K, const std::vector<size_t> &indices, size_t N)
{
	if(m_amghelper.has_positions() == false) return;
	size_t level = m_amghelper.positions.size()-1;

	std::fstream file(filename.c_str(), std::ios::out);

	ConnectionViewer::WriteGridHeader(file, m_amghelper.positions[level], m_amghelper.dimension);

	size_t blockSize = GetRows(mat(0,0));

//		for(size_t i=0; i < m_amghelper.positions[level].size(); i++)
//			file << i << " " << i << " 0.0\n";

	for(size_t i=0; i < N; i++)
	{
		if(blockSize == 1)
			file << indices[i] << " " << indices[i] << " " << v(i, K) << "\n";
		else
		{
			file << indices[i] << " " << indices[i] << " [ ";
			for(size_t j=0; j<blockSize; j++)
				file << v(i*blockSize+j, K) << " ";
			file << "]\n";
		}

	}

	/*for(size_t i=0; i < N; i++)
		for(size_t j=0; j < N; j++)
			if(i != j && BlockNorm(mat(i,j)) > 1e-9)
			{
				file << indices[i] << " " << indices[j] << " " << mat(i,j) << "\n";
			}*/
}
#endif

#endif /* LOCAL_EIGENVALUE_TEST_COMPUTATION_H_ */
