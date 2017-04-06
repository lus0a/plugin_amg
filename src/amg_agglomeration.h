/*
 * Copyright (c) 2012-2015:  G-CSC, Goethe University Frankfurt
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

/**
 * Implementation of agglomeration of cores for AMG.
 */

#include "agglomeration.h"
#ifdef UG_PARALLEL
#include "lib_algebra/parallelization/global_layout.h"
#include "pcl/pcl_tostring.h"
#endif

namespace ug
{

/// \addtogroup plugins_amg
/// \{

#ifdef UG_PARALLEL
//////////////////////////////////////////////////////////////////////////////////////////////////////
template<typename TAlgebra>
bool AMGBase<TAlgebra>::gather_vertical(const vector_type &vec, vector_type &collectedVec, size_t level,
		ParallelStorageType type)
{
	AMGLevelType &L = *levels[level];
	matrix_operator_type &A = *L.pA;
	pcl::InterfaceCommunicator<IndexLayout> &com = A.layouts()->comm();

	if(!levels[level]->bHasBeenMerged)
	{
		UG_COND_THROW(!vec.has_storage_type(type), "storage type is " << vec.get_storage_type() << ", not " << type);
		if(&vec == &collectedVec) return true;
		else
		{
			for(size_t i=0; i<vec.size(); i++)
				collectedVec[i] = vec[i]; // don't use that
		}
		collectedVec.set_storage_type(type);
	}
	else
	{
		bool bRoot = m_agglomerateLevel != level;
		GatherVectorOnOne(L.agglomerateMasterLayout, agglomerateSlaveLayout, com,
				collectedVec, vec, type, bRoot);
		collectedVec.set_storage_type(type);
	}
	return true;
}

//////////////////////////////////////////////////////////////////////////////////////////////////////
template<typename TAlgebra>
bool AMGBase<TAlgebra>::broadcast_vertical(vector_type &vec, const vector_type &collectedVec, size_t level,
		ParallelStorageType type)
{
	AMGLevelType &L = *levels[level];
	matrix_operator_type &A = *L.pA;
	pcl::InterfaceCommunicator<IndexLayout> &com = A.layouts()->comm();

	if(!levels[level]->bHasBeenMerged)
	{
		UG_COND_THROW(!collectedVec.has_storage_type(type), "storage type is " << collectedVec.get_storage_type() << ", not " << type);
		//UG_LOG(collectedVec.get_storage_type() << " " << vec.get_storage_type() << " "<< type << "\n");
		if(&vec == &collectedVec) return true;
		else
		{
			for(size_t i=0; i<vec.size(); i++)
				vec[i] = collectedVec[i]; // don't use that
		}
		vec.set_storage_type(type);
	}
	else
	{
		/// /!\ does not work with type == PST_ADDITIVE (will throw)
		bool bRoot = m_agglomerateLevel != level;
		BroadcastVectorFromOne(L.agglomerateMasterLayout, agglomerateSlaveLayout, com,
				vec, collectedVec, type, bRoot);
		vec.set_storage_type(type);
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
	AMGLevelType &L = *levels[level];


#ifdef UG_PARALLEL
	if(levels[level]->bLevelHasMergers)
	{
		//matrix_operator_type &collA = *L.pAgglomeratedA;
		if(m_bDebugWriteAllVectors)
		{
			writevec(Stringify() << "AMG_" << iteration_glboal << "gd_L", d, level);
			writevec(Stringify() << "AMG_" << iteration_glboal << "gc_L", c, level);
		}
		gather_vertical(d, L.collD, level, PST_ADDITIVE);
		gather_vertical(c, L.collC, level, PST_CONSISTENT);

		if(m_bDebugWriteAllVectors)
		{
			if(!isMergingSlave(level))
			{
				writevec(Stringify() << "AMG_" << iteration_glboal << "coll_gd_L", L.collD, level);
				writevec(Stringify() << "AMG_" << iteration_glboal << "coll_gc_L", L.collC, level);
			}
		}

		// set 0 anyway
		//gather_vertical(c, L.collC, level, PST_CONSISTENT);

		bool b=true;
		if(!isMergingSlave(level))
			b = add_correction_and_update_defect2(L.collC, L.collD, *L.pAgglomeratedA, level);

		if(m_bDebugWriteAllVectors)
		{
			if(!isMergingSlave(level))
			{
				writevec(Stringify() << "AMG_" << iteration_glboal << "coll2_gd_L", L.collD, level);
				writevec(Stringify() << "AMG_" << iteration_glboal << "coll2_gc_L", L.collC, level);
			}
		}

		broadcast_vertical(c, L.collC, level, PST_CONSISTENT);

#if 1
		L.collD.change_storage_type(PST_CONSISTENT);
		broadcast_vertical(d, L.collD, level, PST_CONSISTENT);
		d.change_storage_type(PST_ADDITIVE);
		L.collD.change_storage_type(PST_ADDITIVE);
#else
		// this code is NOT working.
		// broadcast_vertical with PST_ADDITIVE is not possible
	//	broadcast_vertical(d, L.collD, level, PST_ADDITIVE);
#endif

		if(m_bDebugWriteAllVectors)
		{
			writevec(Stringify() << "AMG_" << iteration_glboal << "gd2_L", d, level);
			writevec(Stringify() << "AMG_" << iteration_glboal << "gc2_L", c, level);
			iteration_glboal++;
		}

		return b;
	}
	else
		return add_correction_and_update_defect2(c, d, *L.pA, level);
	}UG_CATCH_THROW("AMGBase::add_correction_and_update_defect failed, level = " << level << " hasMergers = " << (levels[level]->bLevelHasMergers ? "true" : "false"));
#else
	return add_correction_and_update_defect2(c, d, *L.pA, level);
}UG_CATCH_THROW("AMGBase::add_correction_and_update_defect failed, level = " << level);
#endif


}

#ifdef UG_PARALLEL

template<typename TAlgebra>
bool AMGBase<TAlgebra>::agglomerate_from(const std::map<int, int> &merge, const std::vector<int> &mergeWith, size_t level, ParallelNodes &PN)
{
	try{
	AMGLevelType &L = *levels[level];
	matrix_operator_type &A = *L.pA;

	GlobalLayout globalMasterLayout, globalSlaveLayout;
	CreateGlobalLayout(globalMasterLayout, A.layouts()->master(), PN);
	CreateGlobalLayout(globalSlaveLayout, A.layouts()->slave(), PN);

	L.bHasBeenMerged = true;
	IF_DEBUG(DBG_AMG_AGGLOMERATION, 2)
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
	AMGLevelType &L = *levels[level];
	matrix_operator_type &A = *L.pA;

	GlobalLayout globalMasterLayout, globalSlaveLayout;
	CreateGlobalLayout(globalMasterLayout, A.layouts()->master(), PN);
	CreateGlobalLayout(globalSlaveLayout, A.layouts()->slave(), PN);

	L.bHasBeenMerged = true;
	UG_DLOG(DBG_AMG_AGGLOMERATION, 2, "MERGING to " << mergeTo << "\n");

	agglomerateSlaveLayout.clear();
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
	AMGLevelType &L = *levels[level];
	matrix_operator_type &A = *L.pA;

	GlobalLayout globalMasterLayout, globalSlaveLayout;
	CreateGlobalLayout(globalMasterLayout, A.layouts()->master(), PN);
	CreateGlobalLayout(globalSlaveLayout, A.layouts()->slave(), PN);

	L.bHasBeenMerged = false;
	UG_DLOG(DBG_AMG_AGGLOMERATION, 2, "Not merging.\n");
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
	AMGLevelType &L = *levels[level];
	matrix_operator_type &A = *L.pA;

	if(m_levelInformation[level].get_nr_of_nodes_min() > m_minNodesOnOneProcessor
			|| A.layouts()->proc_comm().size() == 1)
	{
		L.bHasBeenMerged = false;
		L.pAgglomeratedA = L.pA;
		L.bLevelHasMergers = false;
		L.spAgglomeratedLayouts = L.spLayouts;

		UG_DLOG(DBG_AMG_AGGLOMERATION, 1, "No Agglomeration. Minimal nodes over all processes is " << m_levelInformation[level].get_nr_of_nodes_min()
					<< " < " << m_minNodesOnOneProcessor << ".\n");
		return true;
	}

	UG_DLOG(DBG_AMG_AGGLOMERATION, 1, "Minimal nodes over all processes is " << m_levelInformation[level].get_nr_of_nodes_min()
			<< " < " << m_minNodesOnOneProcessor << ": agglomerating.\n"
			<< pcl::ToString(A.layouts()->proc_comm()));


	L.bLevelHasMergers=true;

	// write unagglomerated matrix
	if(m_writeMatrices)
		AMGWriteToFile(A, level, level,
			GetAMGFilenamePath(A, level, "A_unagglomerated", level, ".mat", true), m_amghelper);


	// 1. create global Algebra IDs

	pcl::InterfaceCommunicator<IndexLayout> &communicator = A.layouts()->comm();
	const pcl::ProcessCommunicator &pc = A.layouts()->proc_comm();

	ParallelNodes PN(A.layouts(), A.num_rows());

	// 1. ein prozessor erhaelt alle daten
	// 2. ein mass definiert, wie leicht man agglomerieren kann:
	// - reduktion des interfaces

	BinaryBuffer rstream;

	if(pcl::ProcRank() == pc.get_proc_id(0))
	{
		UG_DLOG(DBG_AMG_AGGLOMERATION, 2, "Receiving interface information.\n");

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

		sizes[0] = m_levelInformation[level].get_nr_of_master_nodes();
		//sizes[0] = A.num_rows();
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

		IF_DEBUG(DBG_AMG_AGGLOMERATION, 2)
		{
			for(size_t i=0; i<pc.size(); i++)
			{
				UG_DLOG(DBG_AMG_AGGLOMERATION, 2, "Processor nr. " << i << " has PID " << pc.get_proc_id(i) << ", " << sizes[i] << " nodes and connections to ");
				for(std::map<int, size_t>::iterator it = connections[i].begin(); it != connections[i].end(); ++it){
					UG_DLOG(DBG_AMG_AGGLOMERATION, 2, "processor nr. " << (*it).first << " PID " << pc.get_proc_id((*it).first) << " (" << (*it).second << " connections) ");
				}
				UG_DLOG(DBG_AMG_AGGLOMERATION, 2, "\n");
			}
		}

		std::vector<std::vector<int> > mergeWith(pc.size());
		EasyAgglomeration(sizes, connections, mergeWith, m_minNodesOnOneProcessor, m_preferredNodesOnOneProcessor);

		//		for(size_t i=0; i<mergeWith.size(); i++) PRINT_VECTOR(mergeWith[i], "merge with " << i);

		// prevent agglomerating 0 (or the output process)
		// this prevents a lot of output problems
		if(globalToLocalPID.find(0) != globalToLocalPID.end())
		{
			int id0 = globalToLocalPID[0];
			if(mergeWith[id0].size() == 1 && mergeWith[id0][0] != id0 )
			{
				int mt = mergeWith[id0][0];
				// merge list of 0 should be the merge list of mt
				// replace 0 in this list with mt
				// (0 should be first)
				mergeWith[id0] = mergeWith[mt];
				for(size_t i=0; i<mergeWith[0].size(); i++)
				{
					if(mergeWith[id0][i] == id0)
						mergeWith[id0][i] = mt;
					else if(mergeWith[id0][i] == mt)
						mergeWith[id0][i] = id0;
				}

				// mt should merge to 0
				mergeWith[mt].clear();
				mergeWith[mt].push_back(id0);

				// tell all processes wanted to be merged to mt that they should now merge to 0
				for(size_t i=0; i<mergeWith.size(); i++)
					if(mergeWith[i].size() == 1 && mergeWith[i][0] == mt)
						mergeWith[i][0] = id0;
			}
		}
		else
			UG_LOG("ERROR: BUG: 0 not in processes anymore!!!\n");

//		for(size_t i=0; i<mergeWith.size(); i++) PRINT_VECTOR(mergeWith[i], "merge with " << i);

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
		UG_DLOG(DBG_AMG_AGGLOMERATION, 2, "Sending interface information.\n");
		BinaryBuffer sstream;
		size_t s = m_levelInformation[level].get_nr_of_master_nodes();
		//s = A.num_rows();
		Serialize(sstream, s);
		for(IndexLayout::const_iterator iter = A.layouts()->master().begin(); iter != A.layouts()->master().end(); ++iter)
		{
			int pid = A.layouts()->master().proc_id(iter);
			size_t s = A.layouts()->master().interface(iter).size();

			Serialize(sstream, pid);
			Serialize(sstream, s);

			UG_DLOG(DBG_AMG_AGGLOMERATION, 2, "Sending " << s << " bytes to processor " << pid << ".\n");
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
		if(mergeTo == pcl::ProcRank())
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
