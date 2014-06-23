/*
 * amg_level.h
 *
 *  Created on: 12.06.2014
 *      Author: mrupp
 */

#ifndef AMG_LEVEL_H_
#define AMG_LEVEL_H_
namespace ug{

template<typename vector_type, typename prolongation_matrix_type, typename matrix_type>
class AMGLevel
{
	typedef MatrixOperator<matrix_type, vector_type> matrix_operator_type;
public:
	~AMGLevel()
	{
	}
	AMGLevel(int level) :
			pA(NULL), presmoother(NULL), postsmoother(NULL)
#ifdef UG_PARALLEL
					, spLayouts(new AlgebraLayouts), spAgglomeratedLayouts(new AlgebraLayouts)
					, collectedA(new matrix_operator_type)
#endif
	{
#ifdef UG_PARALLEL
		//slaveLayout.clear();
		//masterLayout.clear();
		agglomerateMasterLayout.clear();
		bHasBeenMerged = false;
#endif
	}

	vector_type corr;					///< temporary Vector for storing the correction made on this level
	vector_type cH;						///< temporary Vector for storing rH
	vector_type dH; 					///< temporary Vector for storing eH

	stdvector<bool> is_fine;			///< fine marks for f-smoothing

	prolongation_matrix_type R; 		///< R Restriction Matrices
	prolongation_matrix_type P; 		///< P Prolongation Matrices
	SmartPtr<matrix_operator_type> pA;	///< A Matrices

	SmartPtr<ILinearIterator<vector_type> > presmoother;///< presmoothers for each level
	SmartPtr<ILinearIterator<vector_type> > postsmoother;///< postsmoothers for each level

	/// inverse for diagonal, for f-smoothing
	stdvector<typename block_traits<typename matrix_type::value_type>::inverse_type> m_diagInv;

#ifdef UG_PARALLEL
	SmartPtr<AlgebraLayouts> spLayouts;
	//pcl::InterfaceCommunicator<IndexLayout> com; 						///< the communicator object on this level
	//IndexLayout slaveLayout, masterLayout;							///< layouts

	SmartPtr<AlgebraLayouts> spAgglomeratedLayouts;
	//IndexLayout slaveLayoutAfterAgglomeration, masterLayoutAfterAgglomeration;	//< since layouts can change due to neighbors merging, these are the new layouts
	//pcl::ProcessCommunicator agglomeratedPC;			///< process Communicator if level contains mergers

	// agglomeration
	bool bHasBeenMerged;///< if true, this core is either father or child of a collected group
	// level 0 - m_agglomerateLevel
	pcl::ProcessCommunicator *pProcessCommunicator;	///< process Communicator on this level

	// level 0 - m_agglomerateLevel-1
	IndexLayout agglomerateMasterLayout;

	vector_type collC, collD;
	SmartPtr<matrix_operator_type> collectedA;
	bool bLevelHasMergers;
#else
	int *spLayouts;
#endif
	SmartPtr<IConvergenceCheck<vector_type> > spYConvCheck;
	// in contrast to collC, collD, and collectedA, this can also point to other vectors.
	SmartPtr<matrix_operator_type> pAgglomeratedA;
};

struct AMGLevelConvergenceRates
{
	AMGLevelConvergenceRates()
	{
		totalSteps = 0;
		calls = 0;
		avgReduction = 1.0;
	}
	size_t totalSteps;
	size_t calls;
	double avgReduction;
};
}
#endif /* AMG_LEVEL_H_ */
