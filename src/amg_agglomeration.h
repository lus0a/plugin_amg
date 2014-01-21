/**
 * \file amg_agglomeration.h
 * \author Martin Rupp
 *
 * \date 03.08.2011
 *
 * Implementation of agglomeration of cores for AMG.
 *
 * Goethe-Center for Scientific Computing 2011-2012.
 */

#include "agglomeration.h"
#ifdef UG_PARALLEL
#include "lib_algebra/parallelization/global_layout.h"
#include "pcl/pcl_tostring.h"
#endif

namespace ug
{

extern DebugID AMGAgglomerationDebug; // SetDebugLevel(debugID.RSAMG.Agglomeration, 0)
/// \addtogroup plugins_amg
/// \{

#ifdef UG_PARALLEL
//////////////////////////////////////////////////////////////////////////////////////////////////////
template<typename TAlgebra>
bool AMGBase<TAlgebra>::gather_vertical(vector_type &vec, vector_type &collectedVec, size_t level,
		ParallelStorageType type)
{
	AMGLevel &L = *levels[level];
	matrix_operator_type &A = *L.pA;
	pcl::InterfaceCommunicator<IndexLayout> &com = A.layouts()->comm();

	if(!levels[level]->bHasBeenMerged)
	{
		//UG_LOG("gather_vertical: not merging at level " << level << "\n");
		if(&vec == &collectedVec) return true;
		else
		{
			for(size_t i=0; i<vec.size(); i++)
				collectedVec[i] = vec[i]; // don't use that
		}
	}
	else
	{
		bool bRoot = m_agglomerateLevel != level;
		GatherVectorOnOne(L.agglomerateMasterLayout, agglomerateSlaveLayout, com,
				collectedVec, vec, type, bRoot);
	}
	return true;
}

//////////////////////////////////////////////////////////////////////////////////////////////////////
template<typename TAlgebra>
bool AMGBase<TAlgebra>::broadcast_vertical(vector_type &vec, vector_type &collectedVec, size_t level,
		ParallelStorageType type)
{
	AMGLevel &L = *levels[level];
	matrix_operator_type &A = *L.pA;
	pcl::InterfaceCommunicator<IndexLayout> &com = A.layouts()->comm();

	if(!levels[level]->bHasBeenMerged)
	{
		//UG_LOG("broadcast_vertical: not merging at level " << level << "\n");
		if(&vec == &collectedVec) return true;
		else
		{
			for(size_t i=0; i<vec.size(); i++)
				vec[i] = collectedVec[i]; // don't use that
		}
	}
	else
	{

		bool bRoot = m_agglomerateLevel != level;
		BroadcastVectorFromOne(L.agglomerateMasterLayout, agglomerateSlaveLayout, com,
				vec, collectedVec, type, bRoot);
	}
	return true;
}

//////////////////////////////////////////////////////////////////////////////////////////////////////
template<typename TAlgebra>
bool AMGBase<TAlgebra>::isMergingSlave(size_t level)
{
	return levels[level]->bHasBeenMerged && m_agglomerateLevel == level;
}

//////////////////////////////////////////////////////////////////////////////////////////////////////
template<typename TAlgebra>
bool AMGBase<TAlgebra>::isMergingMaster(size_t level)
{
	return levels[level]->bHasBeenMerged && m_agglomerateLevel != level;
}

//////////////////////////////////////////////////////////////////////////////////////////////////////
template<typename TAlgebra>
bool AMGBase<TAlgebra>::isNotMerging(size_t level)
{
	return levels[level]->bHasBeenMerged && m_agglomerateLevel != level;
}
#endif

//////////////////////////////////////////////////////////////////////////////////////////////////////
template<typename TAlgebra>
bool AMGBase<TAlgebra>::add_correction_and_update_defect(vector_type &c, vector_type &d, size_t level)
{
	try{
	AMGLevel &L = *levels[level];
	matrix_operator_type &A = *L.pA;
	THROW_IF_NOT_EQUAL_3(c.size(), d.size(), A.num_rows());
#ifdef UG_PARALLEL
	if(levels[level]->bLevelHasMergers)
	{
		//matrix_operator_type &collA = *L.pAgglomeratedA;
		gather_vertical(d, L.collD, level, PST_ADDITIVE);

		// set 0 anyway
		//gather_vertical(c, L.collC, level, PST_CONSISTENT);

		bool b=true;
		if(!isMergingSlave(level))
		{
			L.collC.set(0.0);
			b = add_correction_and_update_defect2(L.collC, L.collD, *L.pAgglomeratedA, level);
		}

		broadcast_vertical(c, L.collC, level, PST_CONSISTENT);
		broadcast_vertical(d, L.collD, level, PST_ADDITIVE);

		return b;
	}
	else
		return add_correction_and_update_defect2(c, d, A, level);
	}UG_CATCH_THROW("AMGBase::add_correction_and_update_defect failed, level = " << level << " hasMergers = " << (levels[level]->bLevelHasMergers ? "true" : "false"));
#else
	return add_correction_and_update_defect2(c, d, A, level);
}UG_CATCH_THROW("AMGBase::add_correction_and_update_defect failed, level = " << level);
#endif


}

#ifdef UG_PARALLEL

template<typename TAlgebra>
bool AMGBase<TAlgebra>::agglomerate_from(const std::map<int, int> &merge, const std::vector<int> &mergeWith, size_t level, ParallelNodes &PN)
{
	try{
	AMGLevel &L = *levels[level];
	matrix_operator_type &A = *L.pA;

	GlobalLayout globalMasterLayout, globalSlaveLayout;
	CreateGlobalLayout(globalMasterLayout, A.layouts()->master(), PN);
	CreateGlobalLayout(globalSlaveLayout, A.layouts()->slave(), PN);

	L.bHasBeenMerged = true;
	IF_DEBUG(AMGAgglomerationDebug, 2)
	{
		UG_LOG("I am the merging father for pids ");
		for(size_t i=0; i<mergeWith.size(); i++) { UG_LOG(mergeWith[i] << " "); }
		UG_LOG("\n");
	}


	//PrintGlobalLayout(globalMasterLayout, "globalmasterLayout before");
	//PrintGlobalLayout(globalSlaveLayout, "globalSlaveLayout before");

	ReceiveMatrix(A, *L.collectedA, L.agglomerateMasterLayout, mergeWith, PN);
	ReceiveGlobalLayout(A.layouts()->comm(), mergeWith, globalMasterLayout, globalSlaveLayout);

	//PrintGlobalLayout(globalMasterLayout, "received globalmasterLayout");
	//PrintGlobalLayout(globalSlaveLayout, "received globalSlaveLayout");

	MergeGlobalLayout(globalMasterLayout, merge);
	MergeGlobalLayout(globalSlaveLayout, merge);

	//PrintGlobalLayout(globalMasterLayout, "AM globalmasterLayout");
	//PrintGlobalLayout(globalSlaveLayout, "AM globalSlaveLayout");

	CreateLayoutFromGlobalLayout(L.spAgglomeratedLayouts->master(), globalMasterLayout, PN);
	CreateLayoutFromGlobalLayout(L.spAgglomeratedLayouts->slave(), globalSlaveLayout, PN);

	L.spAgglomeratedLayouts->proc_comm() = A.layouts()->proc_comm().create_sub_communicator(true);

	L.collectedA->set_layouts(L.spAgglomeratedLayouts);
	L.pAgglomeratedA  = L.collectedA;

	L.collC.resize(L.collectedA->num_rows());
	SetParallelVectorAsMatrix(L.collC, *L.collectedA, PST_CONSISTENT);
	L.collD.resize(L.collectedA->num_rows());
	SetParallelVectorAsMatrix(L.collD, *L.collectedA, PST_ADDITIVE);

	/*PRINTPC(L.spAgglomeratedLayouts->proc_comm());
	UG_LOG("L.masterLayout2\n");
	PrintLayout(L.masterLayout2);
	UG_LOG("L.slaveLayout2\n");
	PrintLayout(L.slaveLayout2);
	UG_LOG("level = " << level << ". old size of matrix at this level: " <<L.collectedA.num_rows() << "\n");*/

	if(m_amghelper.has_positions())
		m_amghelper.receive_agglomerate_positions(level, A.layouts()->comm(), L.agglomerateMasterLayout, L.collectedA->num_rows());
	//merging_as_master(L);
	}UG_CATCH_THROW("AMGBase::agglomerate_from failed, level = " << level);
	return true;
}

template<typename TAlgebra>
bool AMGBase<TAlgebra>::agglomerate_to(const std::map<int, int> &merge, int mergeTo, size_t level, ParallelNodes &PN)
{
	try{
	AMGLevel &L = *levels[level];
	matrix_operator_type &A = *L.pA;

	GlobalLayout globalMasterLayout, globalSlaveLayout;
	CreateGlobalLayout(globalMasterLayout, A.layouts()->master(), PN);
	CreateGlobalLayout(globalSlaveLayout, A.layouts()->slave(), PN);

	L.bHasBeenMerged = true;
	UG_DLOG(AMGAgglomerationDebug, 2, "MERGING to " << mergeTo << "\n");

	SendMatrix(A, agglomerateSlaveLayout, mergeTo, PN);

	MergeGlobalLayout(globalMasterLayout, merge);
	MergeGlobalLayout(globalSlaveLayout, merge);

	// remove interfaces to mergeTo, so we do not unnecessary send them
	globalMasterLayout.erase(mergeTo);
	globalSlaveLayout.erase(mergeTo);

	//PrintGlobalLayout(globalMasterLayout, "send globalmasterLayout");
	//PrintGlobalLayout(globalSlaveLayout, "send globalSlaveLayout");
	SendGlobalLayout(A.layouts()->comm(), globalMasterLayout, globalSlaveLayout, mergeTo);

	L.spAgglomeratedLayouts->proc_comm() = A.layouts()->proc_comm().create_sub_communicator(false);
	/*PRINTPC(L.spAgglomeratedLayouts->proc_comm());
	UG_LOG("agglomerateSlaveLayout\n");
	PrintLayout(agglomerateSlaveLayout);
	UG_LOG("MERGING done \n");*/

	/*L.masterLayout.clear();
	L.slaveLayout.clear();
	A.set_layouts(L.masterLayout, L.slaveLayout);
	A.set_process_communicator(L.spAgglomeratedLayouts->proc_comm());*/
	m_agglomerateLevel = level;

	L.pAgglomeratedA  = SPNULL;
	if(m_amghelper.has_positions())
		m_amghelper.send_agglomerate_positions(level, A.layouts()->comm(), agglomerateSlaveLayout);
	//merging_as_slave(L);
	}UG_CATCH_THROW("AMGBase::agglomerate_to failed, mergeTo = " << mergeTo << ", level = " << level);
	return true;
}

template<typename TAlgebra>
bool AMGBase<TAlgebra>::agglomerate_adjust_interfaces(const std::map<int, int> &merge, size_t level, ParallelNodes &PN)
{
	try{
	AMGLevel &L = *levels[level];
	matrix_operator_type &A = *L.pA;

	GlobalLayout globalMasterLayout, globalSlaveLayout;
	CreateGlobalLayout(globalMasterLayout, A.layouts()->master(), PN);
	CreateGlobalLayout(globalSlaveLayout, A.layouts()->slave(), PN);

	L.bHasBeenMerged = false;
	UG_DLOG(AMGAgglomerationDebug, 2, "Not merging.\n");
	L.spAgglomeratedLayouts->proc_comm() = A.layouts()->proc_comm().create_sub_communicator(true);

	MergeGlobalLayout(globalMasterLayout, merge);
	MergeGlobalLayout(globalSlaveLayout, merge);
	CreateLayoutFromGlobalLayout(L.spAgglomeratedLayouts->master(), globalMasterLayout, PN);
	CreateLayoutFromGlobalLayout(L.spAgglomeratedLayouts->slave(), globalSlaveLayout, PN);

	//PrintGlobalLayout(globalMasterLayout, "AM globalmasterLayout");
	//PrintGlobalLayout(globalSlaveLayout, "AM globalSlaveLayout");

	//PRINTPC(L.spAgglomeratedLayouts->proc_comm());
	/*UG_LOG("others merged, might have to change Layouts:");
	UG_LOG("masterLayout2\n");
	PrintLayout(L.masterLayout2);
	UG_LOG("slaveLayout2\n");
	PrintLayout(L.masterLayout2);*/

	// todo: this is some waste because we're just changing layouts.
	// change vecs tooo....
	*L.collectedA = A;
	L.collectedA->set_layouts(L.spAgglomeratedLayouts);
	L.pAgglomeratedA  = L.collectedA;

	L.collC.resize(L.collectedA->num_rows());
	SetParallelVectorAsMatrix(L.collC, *L.collectedA, PST_CONSISTENT);
	L.collD.resize(L.collectedA->num_rows());
	SetParallelVectorAsMatrix(L.collD, *L.collectedA, PST_ADDITIVE);
	}UG_CATCH_THROW("AMGBase::agglomerate_adjust_interfaces failed, level = " << level);
	return true;
}


template<typename TAlgebra>
bool AMGBase<TAlgebra>::agglomerate(size_t level)
{
	try{
	//UG_ASSERT(0, "not working");
	AMGLevel &L = *levels[level];
	matrix_operator_type &A = *L.pA;

	if(m_levelInformation[level].get_nr_of_nodes_min() > m_minNodesOnOneProcessor
			|| A.layouts()->proc_comm().size() == 1)
	{
		L.bHasBeenMerged = false;
		L.pAgglomeratedA = L.pA;
		L.bLevelHasMergers = false;
		L.spAgglomeratedLayouts = L.spLayouts;
		return true;
	}

	UG_DLOG(AMGAgglomerationDebug, 1, "Minimal nodes over all processes is " << m_levelInformation[level].get_nr_of_nodes_min()
			<< " < " << m_minNodesOnOneProcessor << ": agglomerating.\n"
			<< pcl::ToString(A.layouts()->proc_comm()));


	L.bLevelHasMergers=true;

	// write unagglomerated matrix
	AMGWriteToFile(A, level, level,
			GetAMGFilenamePath(A, level, "A_unagglomerated", level, ".mat", true),
			m_amghelper);


	// 1. create global Algebra IDs

	pcl::InterfaceCommunicator<IndexLayout> &communicator = A.layouts()->comm();
	const pcl::ProcessCommunicator &pc = A.layouts()->proc_comm();

	ParallelNodes PN(A.layouts(), A.num_rows());

	// 1. ein prozessor erhaelt alle daten
	// 2. ein mass definiert, wie leicht man agglomerieren kann:
	// - reduktion des interfaces

	BinaryBuffer rstream;

	if(pcl::GetProcRank() == pc.get_proc_id(0))
	{
		UG_DLOG(AMGAgglomerationDebug, 2, "Receiving interface information.\n");

		// receive interface information and
		typedef std::map<int, BinaryBuffer> BufferMap;
		std::vector<BinaryBuffer> receivepack;
		receivepack.resize(pc.size());
		for(size_t i=1; i<pc.size(); i++)
			communicator.receive_raw(pc.get_proc_id(i), receivepack[i]);

		communicator.communicate();

		// the participating processors can be something like 5, 3, 22, 7
		// so we need to map 5 -> 0, 3 -> 1, 22 -> 2 and 7 -> 3
		// global -> "local"
		// so that agglomeration algorithms are easier to implement.

		std::map<int, int> globalToLocalPID;
		for(size_t i=0; i<pc.size(); i++)
		{
			globalToLocalPID[pc.get_proc_id(i)] = i;
			//UG_LOG("global " << pc.get_proc_id(i) << " local " << i << "\n");
		}

		std::vector<size_t> sizes;
		sizes.resize(pc.size());
		std::vector<std::map<int, size_t> > connections;
		connections.resize(pc.size());

		sizes[0] = A.num_rows();
		for(IndexLayout::const_iterator iter = A.layouts()->master().begin(); iter != A.layouts()->master().end(); ++iter)
		{
			size_t s = A.layouts()->master().interface(iter).size();
			int pid = A.layouts()->master().proc_id(iter);
			int localID = globalToLocalPID[pid];
			connections[0][localID] += s;
			connections[localID][0] += s;
		}

		for(size_t i=1; i<pc.size(); i++)
		{
			BinaryBuffer &stream = receivepack[i];

			size_t s; int pid;
			Deserialize(stream, s);
			sizes[i] = s;
			while(!stream.eof())
			{
				Deserialize(stream, pid);
				Deserialize(stream, s);
				int localID = globalToLocalPID[pid];
				connections[localID][i] += s;
				connections[i][localID] += s;
			}
		}

		IF_DEBUG(AMGAgglomerationDebug, 2)
		{
			for(size_t i=0; i<pc.size(); i++)
			{
				UG_DLOG(AMGAgglomerationDebug, 2, "Processor nr. " << i << " has PID " << pc.get_proc_id(i) << ", " << sizes[i] << " nodes and connections to ");
				for(std::map<int, size_t>::iterator it = connections[i].begin(); it != connections[i].end(); ++it){
					UG_DLOG(AMGAgglomerationDebug, 2, "processor nr. " << (*it).first << " PID " << pc.get_proc_id((*it).first) << " (" << (*it).second << " connections) ");
				}
				UG_DLOG(AMGAgglomerationDebug, 2, "\n");
			}
		}

		std::vector<std::vector<int> > mergeWith(pc.size());
		EasyAgglomeration(sizes, connections, mergeWith, m_minNodesOnOneProcessor, m_preferredNodesOnOneProcessor);


		for(size_t i=0; i<pc.size(); i++)
		{
			// send each processor
			// size
			// if size > 0
			// 		size = nr of processors which merge with you
			//		processors [1..n]
			// else
			//		master processor to merge to
			// endif
			// connected neighbors which merge with others

			BinaryBuffer stream;

			std::vector<int> mergeWithGlobal;
			for(size_t j=0; j<mergeWith[i].size(); j++)
				mergeWithGlobal.push_back(pc.get_proc_id(mergeWith[i][j]));
			Serialize(stream, mergeWithGlobal);

			std::map<int, int> merge; // maps old -> new

			for(std::map<int, size_t>::iterator it = connections[i].begin(); it != connections[i].end(); ++it)
			{
				int j = it->first;
				if(mergeWith[j][0] != j)
					merge[pc.get_proc_id(j)] = pc.get_proc_id(mergeWith[j][0]);
			}

			Serialize(stream, merge);

			if(i == 0)
				rstream = stream;
			else
				communicator.send_raw(pc.get_proc_id(i),
					stream.buffer(), stream.write_pos(), false);

		}
		communicator.communicate();
	}
	else
	{
		UG_DLOG(AMGAgglomerationDebug, 2, "Sending interface information.\n");
		BinaryBuffer sstream;
		Serialize(sstream, A.num_rows());
		for(IndexLayout::const_iterator iter = A.layouts()->master().begin(); iter != A.layouts()->master().end(); ++iter)
		{
			int pid = A.layouts()->master().proc_id(iter);
			size_t s = A.layouts()->master().interface(iter).size();

			Serialize(sstream, pid);
			Serialize(sstream, s);

			UG_DLOG(AMGAgglomerationDebug, 2, "Sending " << s << " bytes to processor " << pid << ".\n");
		}
		communicator.send_raw(pc.get_proc_id(0), sstream.buffer(), sstream.write_pos(), false);
		communicator.communicate();

		communicator.receive_raw(pc.get_proc_id(0), rstream);
		communicator.communicate();
	}

	// receive informations
	// 1. - should i merge and i am master? -> list of merging processors
	//    - should i merge and i am slave? -> master processor
	//    - i shouldnt merge
	// 2. which processors are active afterwards

	matrix_operator_type M;

	std::vector<int> mergeWith;
	Deserialize(rstream, mergeWith);

	std::map<int, int> merge; // maps old -> new
	Deserialize(rstream, merge);

	/*UG_LOG("Neighbor merge information: ");
	for(std::map<int, int>::iterator it = merge.begin(); it != merge.end(); ++it) UG_LOG(it->first << " -> " << it->second);
	UG_LOG("\n");*/


	if(mergeWith.size() > 1)
	{
		mergeWith.erase(mergeWith.begin());
		agglomerate_from(merge, mergeWith, level, PN);
	}
	else
	{
		int mergeTo = mergeWith[0];
		if(mergeTo == pcl::GetProcRank())
			agglomerate_adjust_interfaces(merge, level, PN);
		else
			agglomerate_to(merge, mergeTo, level, PN);
	}

	}UG_CATCH_THROW("AMGBase::agglomerate failed, level = " << level);
	return true;
}
#endif

// end group plugins_amg
/// \}

}
