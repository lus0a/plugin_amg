/*
 * amg_level_information.h
 *
 *  Created on: 07.05.2013
 *      Author: mrupp
 */

#ifndef AMG_LEVEL_INFORMATION_H_
#define AMG_LEVEL_INFORMATION_H_

#include <sstream>

namespace ug
{
class LevelInformation
{
	public:
	LevelInformation(int _level = 0)
	{
		level = _level;
		m_dCoarseningRate = 0.0;
		m_dCreationTimeMS = 0.0;
		m_iNrOfNodes = 0;
		m_iNrOfNodesMin = 0;
		m_iNrOfNodesMax = 0;
		m_iNrOfNodesSum = 0;
		m_iNNZs =0;
		m_iNNZsMin = 0;
		m_iNNZsMax = 0;
		m_iNNZsSum = 0;
		m_iInterfaceElements = 0;
		m_connectionsMax = 0;
		m_iParticipating = 0;
	}
	double get_creation_time_ms() const
	{
		return m_dCreationTimeMS;
	}

	size_t get_nr_of_nodes_min() const
	{
		return m_iNrOfNodesMin;
	}
	size_t get_nr_of_nodes_max() const
	{
		return m_iNrOfNodesMax;
	}
	size_t get_nr_of_nodes() const
	{
		return m_iNrOfNodesSum;
	}
	size_t get_nnz() const
	{
		return m_iNNZsSum;
	}
	size_t get_nnz_min() const
	{
		return m_iNNZsMin;
	}
	size_t get_nnz_max() const
	{
		return m_iNNZsMax;
	}

	// nr of elements in master/slave interfaces (including multiplicities)
	size_t get_nr_of_interface_elements() const
	{
		return m_iInterfaceElements;
	}

	double get_fill_in() const
	{
		return ((double) m_iNNZsSum)
				/ (((double) m_iNrOfNodesSum) * ((double) m_iNrOfNodesSum));
	}
	double get_avg_nnz_per_row() const
	{
		return m_iNNZsSum / (double) m_iNrOfNodesSum;
	}
	size_t get_max_connections() const
	{
		return m_connectionsMax;
	}
	bool is_valid() const
	{
		return this != NULL;
	}

public:
	double get_coarsening_rate() const
	{
		return m_dCoarseningRate;
	}

	int get_participating()
	{
		return m_iParticipating;
	}

#ifdef UG_PARALLEL
	void calc_parallel_data(const pcl::ProcessCommunicator &pc)
	{
		m_dCreationTimeMS = pc.allreduce(m_dCreationTimeMS, PCL_RO_MAX);
		m_iNrOfNodesMin = pc.allreduce(m_iNrOfNodes, PCL_RO_MIN);
		m_iNrOfNodesMax = pc.allreduce(m_iNrOfNodes, PCL_RO_MAX);
		m_iNrOfNodesSum = pc.allreduce(m_iNrOfNodes, PCL_RO_SUM);

		m_iNNZsMin = pc.allreduce(m_iNNZs, PCL_RO_MIN);
		m_iNNZsMax = pc.allreduce(m_iNNZs, PCL_RO_MAX);
		m_iNNZsSum = pc.allreduce(m_iNNZs, PCL_RO_SUM);

		m_iInterfaceElements = pc.allreduce(m_iInterfaceElements, PCL_RO_MAX);
		m_connectionsMax = pc.allreduce(m_connectionsMax, PCL_RO_MAX);
	}

	static void parallel_update(const pcl::ProcessCommunicator &pc,
			std::vector<LevelInformation> &levelInformation, size_t totalUsedLevels)
	{
		for(size_t i=0; i<totalUsedLevels; i++)
		{
			if(levelInformation.size() < i+1)
				levelInformation.resize(i+1);

			LevelInformation &li = levelInformation[i];
			li.m_iParticipating=pc.allreduce(li.m_iParticipating, PCL_RO_MAX);
			li.level = i;
			li.calc_parallel_data(pc);
			if(i>0)
				li.m_dCoarseningRate = (double)li.get_nr_of_nodes()/
						((double)levelInformation[i-1].get_nr_of_nodes());
			else
				li.m_dCoarseningRate = 0.0;
		}
	}
#endif

	template<typename matrix_type>
	void calculate(matrix_type &A, int _level, double dCreationTime, LevelInformation *previousLevel)
	{
		level = _level;
		m_dCreationTimeMS = dCreationTime;
		m_iNNZs = GetNNZs(A); // A.total_num_connections();
		m_connectionsMax = GetMaxConnections(A);

#ifdef UG_PARALLEL
		const pcl::ProcessCommunicator &pc = A.layouts()->proc_comm(); //!
		m_iNrOfNodes = A.num_rows() - A.layouts()->slave().num_interface_elements();
		size_t localInterfaceElements = A.layouts()->master().num_interface_elements();
		m_iInterfaceElements = pc.allreduce(localInterfaceElements, PCL_RO_SUM);
		m_iParticipating=pc.size();
		calc_parallel_data(pc);
#else
		m_iNrOfNodes = A.num_rows();
		m_iInterfaceElements = 0;
		m_iParticipating=1;
		m_iNrOfNodesMin = m_iNrOfNodesMax = m_iNrOfNodesSum = m_iNrOfNodes;
		m_iNNZsMin = m_iNNZsMax = m_iNNZsSum = m_iNNZs;
#endif

		if(previousLevel != NULL)
		{
			m_dCoarseningRate = m_iNrOfNodesSum/((double)previousLevel->m_iNrOfNodesSum);
		}
		else
		{
			m_dCoarseningRate = 0.0;
		}

	}

	std::string tostring() const
	{
		std::stringstream ss;
		ss << "Level " << level << ": creation time: "
				<< get_creation_time_ms() << " ms. number of nodes: "
				<< get_nr_of_nodes() << ". fill in " << get_fill_in() * 100
				<< "%. ";
		if (level != 0)
			ss << "coarsening rate: " << get_coarsening_rate() * 100.0
					<< "%. ";
		ss << "nr of interface elements: " << get_nr_of_interface_elements()
				<< " ("
				<< (double) get_nr_of_interface_elements()
						/ (double) get_nr_of_nodes() * 100.0 << "%) "
				<< "nnzs: " << get_nnz() << " avgNNZs/row: "
				<< get_avg_nnz_per_row() << " maxCon: "
				<< get_max_connections();
	#ifdef UG_PARALLEL
		ss << " participating: " << m_iParticipating;
	#endif
		return ss.str();
	}



public:
	int level;
	double m_dCoarseningRate;
	double m_dCreationTimeMS;
	size_t m_iNrOfNodes;
	size_t m_iNNZs;

	size_t m_iNrOfNodesMin;
	size_t m_iNrOfNodesMax;
	size_t m_iNrOfNodesSum;

	size_t m_iNNZsMin;
	size_t m_iNNZsMax;
	size_t m_iNNZsSum;

	size_t m_iInterfaceElements;
	size_t m_connectionsMax;
	size_t m_iParticipating;
};

}
#endif /* AMG_LEVEL_INFORMATION_H_ */
