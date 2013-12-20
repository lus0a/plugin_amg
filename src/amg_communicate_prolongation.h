/**
 * \file amg_communicate_prolongation.h
 *  Created on: 27.04.2011
 *      Author: mrupp
 */

#ifndef AMG_COMMUNICATE_PROLONGATION_H_
#define AMG_COMMUNICATE_PROLONGATION_H_

#ifdef UG_PARALLEL
#include "pcl/pcl.h"
#include "lib_algebra/parallelization/row_sending_scheme.h"
#include "lib_algebra/parallelization/communication_scheme.h"
#endif

// note: documentation is in amg_base.h

namespace ug
{

/// \addtogroup plugins_amg
/// \{

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// create_fine_marks
//------------------------------
// docu see amg_base.h
template<typename TAlgebra>
template<typename TAMGNodes>
void AMGBase<TAlgebra>::create_fine_marks(int level, TAMGNodes &amgnodes, size_t N)
{
	AMG_PROFILE_FUNC();
	stdvector<bool> &vFine = levels[level]->is_fine;
	vFine.resize(N);

	for(size_t i=0; i < N; i++)
		vFine[i] = amgnodes[i].is_fine();
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// create_parent_index
//------------------------------
// docu see amg_base.h
template<typename TAlgebra>
void AMGBase<TAlgebra>::create_parent_index(int level, stdvector<int> newIndex, size_t nrOfCoarse)

{
	AMG_PROFILE_FUNC();
	m_parentIndex.resize(level+2);
	m_parentIndex[level+1].resize(nrOfCoarse);
	for(size_t i=0; i < nrOfCoarse; i++) m_parentIndex[level+1][i] = -1;
	for(size_t i=0; i < newIndex.size(); i++)
		if(newIndex[i] != -1)
			m_parentIndex[level+1][ newIndex[i] ] = i;
}


///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// create_new_index
//------------------------------
// docu see amg_base.h
template<typename TAlgebra>
template<typename TAMGNodes>
void AMGBase<TAlgebra>::create_new_indices(prolongation_matrix_type &PoldIndices, prolongation_matrix_type &PnewIndices,
		size_t N, TAMGNodes &amgnodes, stdvector<int> &newIndex, double dEpsilonTruncation)

{
	try{
	AMG_PROFILE_FUNC();
//	create new indices
	newIndex.clear();
	newIndex.resize(PoldIndices.num_cols(), -1);
	size_t nrOfCoarse=0;
	for(size_t r=0; r<N; r++)
	{
		for(typename prolongation_matrix_type::row_iterator it = PoldIndices.begin_row(r);
				it != PoldIndices.end_row(r); ++it)
		{
			size_t c = it.index();
			UG_ASSERT(amgnodes[c].is_coarse(), c);
			if(newIndex[c] == -1)
				newIndex[c] = nrOfCoarse++;
		}
	}

//	P to new indices
	PnewIndices.resize_and_clear(N, nrOfCoarse);
	for(size_t r=0; r<N; r++)
	{
		double maxCon=0;
		for(typename prolongation_matrix_type::row_iterator conn = PoldIndices.begin_row(r);
			conn != PoldIndices.end_row(r); ++conn)
		{
			double con = BlockNorm(conn.value());
			maxCon = std::max(con, maxCon);
		}

		for(typename prolongation_matrix_type::row_iterator conn = PoldIndices.begin_row(r);
				conn != PoldIndices.end_row(r); ++conn)
		{
			if(BlockNorm(conn.value()) < maxCon*dEpsilonTruncation)
				continue;
			size_t c = conn.index();
			UG_ASSERT(amgnodes[c].is_coarse(), "node " << c << " is not even coarse.");
			UG_ASSERT(newIndex[c] != -1, amgnodes[c]);
			// assure we are only interpolating Master0 and Slave0 nodes from Master0, Slave0 or Slave1 nodes.
			//UG_ASSERT(PN.distance_to_master_or_inner(c) <= 1
				//	|| (PN.is_slave(r) && PN.distance_to_master_or_inner(c) == 2),
					//c << " = " << PN.overlap_type());
			PnewIndices(r, newIndex[c]) = conn.value();
		}
	}

#ifdef UG_DEBUG
	/*for(size_t r=N; r<overlapN; r++)
		UG_ASSERT(amgnodes[r].is_fine() == false || amgnodes[r].is_uncalculated_fine(), amgnodes.info(r));*/
#endif
	//UG_DLOG(LIB_ALG_AMG, 0, "\namgnodes.get_nr_of_coarse() = " << amgnodes.get_nr_of_coarse() << ", nrOfCoarse = " << nrOfCoarse);
	PnewIndices.defragment();
#ifdef UG_PARALLEL
	PnewIndices.set_storage_type(PST_CONSISTENT);
#endif
	}UG_CATCH_THROW("AMGBase::create_new_indices failed");
}


#ifndef UG_PARALLEL

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// serial_process_prolongation:
//--------------------------------
// docu see amg_base.h
template<typename TAlgebra>
template<typename TAMGNodes>
void AMGBase<TAlgebra>::serial_process_prolongation(prolongation_matrix_type &PoldIndices,
		prolongation_matrix_type &PnewIndices, double dEpsilonTruncation, int level, TAMGNodes &amgnodes, stdvector<int> &newIndex)
{
	try{
	AMG_PROFILE_FUNC();

	create_new_indices(PoldIndices, PnewIndices, PoldIndices.num_rows(), amgnodes, newIndex, dEpsilonTruncation);
	create_parent_index(level, newIndex, PnewIndices.num_cols());

	if(needs_positions())
		m_amghelper.make_coarse_level(level+1, m_parentIndex[level+1]);
	create_fine_marks(level, amgnodes, PoldIndices.num_rows());

	}UG_CATCH_THROW("AMGBase::serial_process_prolongation failed, level = " << level);
}

#else


///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// CheckMatrixLayouts
//------------------------------
/**
 * this function checks:
 * - if ( P(i,j) != 0.0) and is_slave(j)), that layouts.slave().interface(masterProc(j)).contains(j)
 * - tests that master/slave layouts are double-ended and matching.
 *
 * @param PN		ParallelNodes used for local->global indexing
 * @param mat		checked matrix. note that we check column entries only (the 'Preimage' of mat).
 * @param layouts	used layouts
 */
template<typename TMatrix>
void CheckMatrixLayouts(ParallelNodes &PN, const TMatrix &mat, AlgebraLayouts &layouts)
{
#ifndef NDEBUG
	AMG_PROFILE_FUNC();
	// check slave layout
	// go to every connection. if node not master on this processor, check that slave interface exists
	bool bEverySlaveIsInLayout=true;
	for(size_t r = 0; r<PN.get_original_size(); r++)
	{
		for(typename TMatrix::const_row_iterator it = mat.begin_row(r); it != mat.end_row(r); ++it)
		{
			int c = it.index();
			const AlgebraID &globalID = PN.local_to_global(c);
			int mPID = globalID.master_proc();
			if(mPID == pcl::GetProcRank()) continue;
			if(layouts.slave().interface_exists(mPID) == false)
			{
				UG_LOG("Node " << c << "(gid " << globalID << ") has no connection to master (interface to master proc "
						<< mPID << " nonexistent)\n");
				bEverySlaveIsInLayout = false;
			}
			else if(!IsInInterface(layouts.slave().interface(mPID), c))
			{
				UG_LOG("Node " << c << "(gid " << globalID << ") has no connection to master (not in interface to master proc "
						<< mPID << ")\n");
				bEverySlaveIsInLayout = false;
			}
		}
	}

	if(!bEverySlaveIsInLayout)
	{
		PRINTPC(layouts.proc_comm());
		PrintLayouts(layouts);
		mat.print();
	}
	UG_ASSERT(bEverySlaveIsInLayout, "Error in Interfaces.");

	// check if slave interfaces match the master interfaces.
	TESTLAYOUTS(layouts);
#endif
}
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// parallel_process_prolongation:
//--------------------------------
// docu see amg_base.h
template<typename TAlgebra>
template<typename TAMGNodes>
void AMGBase<TAlgebra>::parallel_process_prolongation(prolongation_matrix_type &PoldIndices,
		prolongation_matrix_type &PnewIndices, double dEpsilonTruncation, int level,
		TAMGNodes &amgnodes, ParallelNodes &PN, bool bCreateNewNodes,
		AlgebraLayouts &nextLevelLayouts, stdvector<int> &newIndex)
{
	try{
	AMG_PROFILE_FUNC();
	//PRINTPC(PoldIndices.layouts()->proc_comm());
	//PRINTPC(PnewIndices.layouts()->proc_comm());

	// send rows of PoldIndices from master to slave, update P, PN and layouts
	communicate_prolongation_additive_to_consistent(PN, PoldIndices, bCreateNewNodes);


	postset_coarse(PN, PoldIndices, amgnodes);
	create_minimal_layout_for_prolongation(PN, PoldIndices, nextLevelLayouts);

	CheckMatrixLayouts(PN, PoldIndices, nextLevelLayouts);

	create_new_indices(PoldIndices, PnewIndices, PN.get_original_size(), amgnodes, newIndex,
			dEpsilonTruncation);
	PnewIndices.set_storage_type(PST_CONSISTENT);

	ReplaceIndicesInLayout(nextLevelLayouts.master(), newIndex);
	ReplaceIndicesInLayout(nextLevelLayouts.slave(), newIndex);

	TESTLAYOUTS(nextLevelLayouts);

	create_parent_index(level, newIndex, PnewIndices.num_cols());

	if(needs_positions())
	{
		m_amghelper.update_overlap_positions(level, PN.communicator(), PN.get_total_master_layout(),
				PN.get_total_slave_layout(), PN.local_size());
		m_amghelper.make_coarse_level(level+1, m_parentIndex[level+1]);
	}
	create_fine_marks(level, amgnodes, PN.get_original_size());

	}UG_CATCH_THROW("AMGBase::parallel_process_prolongation failed, level = " << level);
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// communicate_prolongation_additive_to_consistent:
//------------------------------
// docu see amg_base.h
template<typename TAlgebra>
void AMGBase<TAlgebra>::communicate_prolongation_additive_to_consistent(ParallelNodes &PN, prolongation_matrix_type &PoldIndices, bool bCreateNewNodes)
{
	try{
	AMG_PROFILE_FUNC();

	RowSendingScheme<prolongation_matrix_type> rowSendingScheme(PoldIndices, PN);

	rowSendingScheme.set_create_new_nodes(bCreateNewNodes);
	rowSendingScheme.issue_send(PN.communicator(), PN.master_layout(), PN.slave_layout());
	PN.communicator().communicate();

	rowSendingScheme.process(PN.slave_layout());
	rowSendingScheme.set_rows_in_matrix(PoldIndices);

	// we then have a consistent PoldIndices
	PoldIndices.set_storage_type(PST_CONSISTENT);

	}UG_CATCH_THROW("AMGBase::communicate_prolongation_additive_to_consistent failed");

}


///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// postset_coarse:
//------------------------------
// docu see amg_base.h
template<typename TAlgebra>
template<typename TAMGNodes>
void AMGBase<TAlgebra>::postset_coarse(ParallelNodes &PN, prolongation_matrix_type &PoldIndices, TAMGNodes &amgnodes)
{
	try{
	AMG_PROFILE_FUNC();

	amgnodes.resize(PoldIndices.num_cols());

	// set coarse flags of new rows
	for(IndexLayout::const_iterator iter = PN.slave_layout().begin();
			iter != PN.slave_layout().end(); ++iter)
	{
		const IndexLayout::Interface &interface = PN.slave_layout().interface(iter);
		for(IndexLayout::Interface::const_iterator iter = interface.begin(); iter != interface.end(); ++iter)
		{
			size_t i = interface.get_element(iter);
			UG_ASSERT(PN.is_slave(i), i);

			if(!PoldIndices.is_isolated(i) && amgnodes[i].is_fine() == false && amgnodes[i].is_dirichlet() == false)
			{
				amgnodes.set_fine(i);
				UG_DLOG(LIB_ALG_AMG, 0, "post-setted " << i << " fine.\n");
			}

			// this is a prolongation
			for(typename prolongation_matrix_type::row_iterator it = PoldIndices.begin_row(i);
				it != PoldIndices.end_row(i); ++it)
			{
				// so all connections must be coarse
				size_t j = it.index();
				if(amgnodes[j].is_coarse() == false)
				{
					amgnodes.external_set_coarse(j);
					UG_DLOG(LIB_ALG_AMG, 4, "post-setted " << j << " coarse.\n");
				}
			}
		}
	}

	}UG_CATCH_THROW("AMGBase::postset_coarse failed");
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// CreateMinimalLayoutCommunicationScheme:
//-----------------------------------------
/** used to delete unused elements of interfaces on boths sides
 */
class CreateMinimalLayoutCommunicationScheme
	: public CommunicationScheme<CreateMinimalLayoutCommunicationScheme, bool>
{
public:
	CreateMinimalLayoutCommunicationScheme(stdvector<bool> &bUsedSlave) : m_bUsedSlave(bUsedSlave)
	{
		m_newMasterLayout.clear();
		m_newSlaveLayout.clear();
	}

	inline char send(int pid, size_t index)
	{
		if(index >= m_bUsedSlave.size() || m_bUsedSlave[index] == false)
			return false;
		m_newSlaveLayout.interface(pid).push_back(index);
		return true;
	}

	inline void receive(int pid, size_t index, char b)
	{
		if(b) m_newMasterLayout.interface(pid).push_back(index);
	}

	inline int get_element_size() const
	{
		return sizeof(char);
	}

public:
	void get_new_layouts(IndexLayout &master, IndexLayout &slave)
	{
		master.clear();
		AddLayout(master, m_newMasterLayout);
		slave.clear();
		AddLayout(slave, m_newSlaveLayout);
	}

private:
	stdvector<bool> &m_bUsedSlave;
	IndexLayout m_newMasterLayout, m_newSlaveLayout;
};




///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// create_minimal_layout_for_prolongation
//------------------------------------------
// docu see amg_base.h
template<typename TAlgebra>
void AMGBase<TAlgebra>::create_minimal_layout_for_prolongation(ParallelNodes &PN, prolongation_matrix_type &P,
		AlgebraLayouts &newLayouts)
{
	try{
	AMG_PROFILE_FUNC();
	stdvector<bool> bUsedSlave(P.num_cols(), false);

	for(size_t r=0; r<PN.get_original_size(); r++)
	{
		for(typename prolongation_matrix_type::row_iterator it = P.begin_row(r);
				it != P.end_row(r); ++it)
		{
			size_t c = it.index();
			if(PN[c].is_slave())
				bUsedSlave[c] = true;
		}
	}

	/*for(size_t i=0; i<bUsedSlave.size(); i++)
	{
		UG_LOG(i << " used = " << bUsedSlave[i] << "\n");
	}*/
	CreateMinimalLayoutCommunicationScheme cmlcs(bUsedSlave);
	CommunicateOnInterfaces(PN.communicator(), PN.get_total_slave_layout(), PN.get_total_master_layout(), cmlcs);
	cmlcs.get_new_layouts(newLayouts.master(), newLayouts.slave());

	}UG_CATCH_THROW("AMGBase::create_minimal_layout_for_prolongation failed");

}

#endif

// end group plugins_amg
/// \}

}
#endif /* FAMG_COMMUNICATE_PROLONGATION_H_ */
