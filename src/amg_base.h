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
 * Goethe-Center for Scientific Computing 2010.
 */

#ifndef __H__UG__LIB_ALGEBRA__AMG_SOLVER__AMG_BASE_H__
#define __H__UG__LIB_ALGEBRA__AMG_SOLVER__AMG_BASE_H__

//#include "lib_algebra/lib_algebra.h"
#include <sys/stat.h>
#include "lib_algebra/operator/interface/operator_inverse.h"
#include "lib_algebra/operator/interface/preconditioner.h"
#include "lib_algebra/operator/interface/matrix_operator.h"
#include "lib_algebra/operator/interface/linear_operator_inverse.h"
#include "lib_algebra/operator/vector_writer.h"
#include "common/util/string_util.h"
#include <vector>
#include <iostream>
#include <sstream>
#include <string>
#include "amg_debug_helper.h"
#include "common/util/string_table_stream.h"
#include "amg_level_information.h"
#include "lib_algebra/common/stl_debug.h" // stdvector
#ifdef UG_PARALLEL
#include "pcl/pcl.h"
#include "lib_algebra/parallelization/parallel_nodes.h"
#endif

#include "lib_algebra/operator/fixed_convergence_check.h"
#include "lib_algebra/operator/convergence_check.h"
#include "common/debug_print.h"

#include "common/util/string_util.h"
#include "amg_level.h"

#include "famg/testvector_calculator.h"


#define PRINTLAYOUT(pc, com, Layout1, Layout2) MyPrintLayout(pc, com, (Layout1), (Layout2), #Layout1, #Layout2)

namespace ug {

/**
 * \defgroup plugins_amg Algebraic Multigrid
 * \ingroup plugins_experimental
 * \{
 */

#ifdef UG_PARALLEL
typedef AlgebraLayouts AMGLayouts;
#else
typedef int AMGLayouts;
#endif


//extern DebugID DebugRSAMGProlongation;
#define DebugRSAMGProlongation LIB_ALG_AMG

template<typename TAlgebra>
class AMGBase: public IPreconditioner<TAlgebra> {
public:
	using DebugWritingObject<TAlgebra>::debug_writer;

///	Algebra type
	typedef TAlgebra algebra_type;

///	Vector type
	typedef typename TAlgebra::vector_type vector_type;

///	Matrix type
	typedef typename TAlgebra::matrix_type matrix_type;
///	Prolongation Matrix type

//typedef typename CPUAlgebra::matrix_type prolongation_matrix_type;
	typedef typename TAlgebra::matrix_type prolongation_matrix_type;

///	Matrix Operator type
	typedef MatrixOperator<matrix_type, vector_type> matrix_operator_type;

	typedef typename matrix_type::value_type value_type;
	typedef IPreconditioner<TAlgebra> base_type;
	typedef AMGLevel<vector_type, prolongation_matrix_type, matrix_type> AMGLevelType;

	struct checkResult {
		double preSmoothing;
		double preFSmoothing;
		size_t iInnerIterations;
		double lastCoarseReduction;
		double coarseDefect;
		double coarseCorrection;
		double reduction;
		double postFSmoothing;
		double postSmoothing;
	};

//  functions
	AMGBase();
	/// copy constructor
	AMGBase(const AMGBase<TAlgebra> &parent);
	/// helper function to unify constructor and copy constructor
	void set_or_copy(const AMGBase<TAlgebra> *parent);

	//	Name of preconditioner
	virtual ~AMGBase();
	void cleanup();

/////// pure virtual functions //////

	virtual void create_AMG_level(matrix_type &AH, prolongation_matrix_type &R,
			const matrix_type &A, prolongation_matrix_type &P, size_t level,
			AMGLayouts &nextLevelLayouts) = 0;
	virtual void precalc_level(size_t level) = 0;

	virtual SmartPtr<ILinearIterator<vector_type> > clone() = 0;

	/// overwrite this method if your AMG method needs positions somehow
	virtual bool needs_positions() const { return m_writeMatrices; }

	virtual const char* name() const = 0;

public:
/////// IPreconditioner interface ///////

	///	returns information about configuration parameters
	virtual std::string config_string() const;

	/**
	 * creates MG Hierachy for with matrix_operator_type A and temporary vectors for higher levels
	 * @param A matrix A.
	 * @return true
	 */
	virtual bool preprocess(SmartPtr<MatrixOperator<matrix_type, vector_type> > pOp);

	///	nothing to do here
	virtual bool postprocess() {
		return true;
	}

//	Stepping routine
	virtual bool step(SmartPtr<MatrixOperator<matrix_type, vector_type> > pOp,
			vector_type& c, const vector_type& d)
	{
		UG_COND_THROW(!m_bInited, "not inited?");
#ifdef UG_PARALLEL
		if(m_bEmptyProcess)
		{
			c.set_storage_type(PST_CONSISTENT);
			return true;
		}
#endif
		return get_correction(c, d);
	}
	virtual bool supports_parallel() const { return true; }

protected:

/////// base ///////
	/**
	 * create direct solver on a specified level
	 * @param level
	 */
	void create_direct_solver(size_t level);

	/**
	 * solves A[level]c = d with the direct solver,
	 * then calculates d -= A[level]c
	 * @param c output correction
	 * @param d defect. will be updated
	 * @param level
	 * @return true
	 */
	bool solve_on_base(vector_type &c, vector_type &d, size_t level);


//////////// AMG CONSTRUCTION ////////////
	bool preprocess_create_level(size_t level, double &createAMGlevelTiming, int oldOutputProcess);
	void preprocess_init_positions(matrix_operator_type& mat);
	bool preprocess_check_requirements();

	/**
	 * every coarse node c gets a unique new index newIndex[c]
	 * calculate PnewIndices(r, newIndex[c]) = PoldIndex(r, c)
	 *
	 * @param PoldIndices			prolongation matrix with old indices
	 * @param PnewIndices			prolongation matrix with new indices
	 * @param N						overlap 0 size (size of master+slave+inner nodes)
	 * @param amgnodes				amgnodes used for coarse/fine
	 * @param newIndex				new indices are stored here
	 * @param dEpsilonTruncation	connections in P with P(i,j) < dEpsilonTruction * max_k P(i,k) are dropped.
	 */
	template<typename TAMGNodes>
	void create_new_indices(prolongation_matrix_type &PoldIndices,
			prolongation_matrix_type &PnewIndices, size_t N,
			TAMGNodes &amgnodes, stdvector<int> &newIndex,
			double dEpsilonTruncation);
#ifndef UG_PARALLEL
	/**
	 * this function processes the prolongation in old indices and returns a prolongation with new indices,
	 * computing the new indices and and everything associated with the new indices (fine marks, parent indices etc.)
	 * in detail, this is
	 * - create the new indices by counting used coarse nodes and creation of PnewIndices @sa create_new_indices
	 * - create parent index (for debugging) @sa create_parent_index
	 * - create the coarser level of positions (not only for debugging) @sa make_coarse_level
	 * - create fine marks for f-smoothing @sa create_fine_marks
	 *
	 * @note this function is the serial version of @sa parallel_process_prolongation
	 *
	 * @param 	PoldIndices				used Prolongation matrix with old (fine grid) indices. in/out.
	 * @param	PnewIndices				Prolongation matrix with new indices.
	 * @param	dEpsilonTruncation		used to trucation the prolongation /sa create_new_indices
	 * @param	level					the amg level
	 * @param	amgnodes				coarse/fine ratings of the nodes
	 * @param	nextLevelMasterLayout	the master layout on the next coarser level
	 * @param	nextLevelSLaveLayout	the slave layout on the next coarser level
	 */
	template<typename TAMGNodes>
	void serial_process_prolongation(prolongation_matrix_type &PoldIndices, prolongation_matrix_type &PnewIndices,
			double dEpsilonTruncation, int level, TAMGNodes &amgnodes, stdvector<int> &newIndex);
#else

	/**
	 * this function processes the prolongation in old indices and returns a prolongation with new indices,
	 * computing the new indices and everything associated with the new indices (fine marks, parent indices etc.)
	 * in detail, this is
	 * - since we calculate the prolongation only in the master/inner nodes, we need to send this communication to the slave
	 *  nodes (so we get additive AH=RAP again)
	 * 	 @sa communicate_prolongation
	 * - post set nodes which are interpolated from because of newly received prolongation rows @sa postset_coarse
	 * - create a minimal layout for the next level: we start with the total overlap layout we have in ParallelNodes PN, which
	 *  also includes some node which might have
	 *   been added because of the sending of the prolongation rows, and then send a 'delete' to the master proc
	 *   @sa create_minimal_layout_for_prolongation
	 * - create the new indices by counting used coarse nodes and creation of PnewIndices @sa create_new_indices
	 * - Replace Indices in the layout so we have our next level layouts with new indices @sa ReplaceIndicesInLayout
	 * - create parent index (for debugging) @sa create_parent_index
	 * - create the coarser level of positions (for debugging) @sa make_coarse_level
	 * - create fine marks for f-smoothing @sa create_fine_marks
	 *
	 * \note this function is the parallel version of @sa serial_process_prolongation
	 *
	 * @param 	PoldIndices				used Prolongation matrix with old (fine grid) indices. in/out.
	 * @param	PnewIndices				Prolongation matrix with new indices.
	 * @param	dEpsilonTruncation		used to trucation the prolongation /sa create_new_indices
	 * @param	level					the amg level
	 * @param	amgnodes				coarse/fine ratings of the nodes
	 * @param	nextLevelMasterLayout	the master layout on the next coarser level
	 * @param	nextLevelSLaveLayout	the slave layout on the next coarser level
	 */
	template<typename TAMGNodes>
	void parallel_process_prolongation(prolongation_matrix_type &PoldIndices,
			prolongation_matrix_type &PnewIndices, double dEpsilonTruncation,
			int level, TAMGNodes &amgnodes, ParallelNodes &PN,
			bool bCreateNewNodes, AlgebraLayouts &nextLevelLayout, stdvector<int> &newIndex);

	/** since we calculate the prolongation only in the master/inner nodes, we need to send this communication to the slave nodes
	 *  (so we get additive AH=RAP again)
	 * This function sends the matrix rows of PoldIndices from the master to the slave
	 * and receives them, adds the new nodes which need to be created,
	 * and then sets the newly generated rows
	 * @param	PN				ParallelNodes. here we use master/slaveLayouts
	 * @param 	PoldIndices		used Prolongation matrix with old (fine grid) indices. in/out.
	 * @note this function can create new nodes in ParallelNodes. This is because for example a Slave node may be interpolated by
	 * a node which is in Overlap1.
	 */
	void communicate_prolongation_additive_to_consistent(ParallelNodes &PN,
			prolongation_matrix_type &PoldIndices, bool bCreateNewNodes);

	/**
	 * post set nodes which are interpolated from because of newly received prolongation rows
	 * sets all nodes in amgnodes coarse which are connected to by some node in PoldIndices.
	 * These can only be the neighbors of slave nodes, since these rows come from other processes.
	 *
	 * for all slave layouts S
	 *  	for all i in S
	 *  		for all j with PoldIndices[i, j] != 0
	 *  			set node j coarse
	 *
	 * @param PN			Parallel nodes structure
	 * @param PoldIndices	prolongation with original indices
	  * @param	amgnodes		struct to set fine/coarse.
	 */
	template<typename TAMGNodes>
	void postset_coarse(ParallelNodes &PN,
			prolongation_matrix_type &PoldIndices, TAMGNodes &nodes);

	/**
	 * P is a matrix which has nodes from ParallelNodes PN. We want to get a minimal layout for P.
	 * For this we first mark all slaves j for which there exist i with P(i, j) != 0.
	 * Since the P is P : Coarse -> Fine, this means we mark all slave elements we use on this processor.
	 * Then we create a slave layout out of it an communicate it to the master. This is necessary
	 * since the master cannot know which nodes we use.
	 *
	 * @param PN			Parallel Nodes structure
	 * @param P				used prolongation matrix
	 * @param newLayouts	new layouts for P, condensed from PN.get_total_master_layout() and PN.get_total_slave_layout();
	 */
	void create_minimal_layout_for_prolongation(ParallelNodes &PN,
			prolongation_matrix_type &P, AlgebraLayouts &newLayouts);

	/**
	 * since we can turn off cores, we need to collect the data from them
	 * here we collect the data from vec into collectedVec.
	 * @param vec			uncollected vec
	 * @param collectedVec	collectedVec (result)
	 * @param level			level we operate on
	 * @param type			the type we want the result to be in (PST_ADDITIVE or PST_CONSISTENT)
	 */
	bool gather_vertical(const vector_type &vec, vector_type &collectedVec,
			size_t level, ParallelStorageType type);

	/**
	 * in this function, we broadcast the data to cores we have shut off
	 * i.e. we write the data to vec
	 * @param vec			uncollected vec (result)
	 * @param collectedVec	collectedVec
	 * @param level			level we operate on
	 * @param type			the type we want the result to be in (PST_ADDITIVE or PST_CONSISTENT)
	 */
	bool broadcast_vertical(vector_type &vec, const vector_type &collectedVec,
			size_t level, ParallelStorageType type);

#endif


/////// public set/get functions ///////
public:
	void set_num_smooth(size_t n)
	{
		set_num_presmooth(n);
		set_num_postsmooth(n);
	}

	void set_num_presmooth(size_t new_presmooth) 				{ m_numPreSmooth = new_presmooth; }
	size_t get_num_presmooth() const 							{ return m_numPreSmooth; }

	void set_num_postsmooth(size_t new_postsmooth) 				{ m_numPostSmooth = new_postsmooth; }
	size_t get_num_postsmooth() const 							{ return m_numPostSmooth; }

	void set_cycle_type(size_t new_cycletype)
	{
		UG_COND_THROW(new_cycletype <= 0, "cannot set cycle type " << new_cycletype << ", has to be > 0");
		set_Y_check(make_sp(new FixedConvergenceCheck<vector_type>(new_cycletype)));
		if(new_cycletype == 1)
			m_cycleType = "V-Cycle";
		else if(new_cycletype == 2)
			m_cycleType = "W-Cycle";
		else
			m_cycleType = "Gamma=" + ToString(new_cycletype) + "-Cycle";
	}
	void set_Y_cycle(int maxIterations, double dYminAbsDefect, double dYrelReduction)
	{
		set_Y_check(make_sp(new StdConvCheck<vector_type>(maxIterations, dYminAbsDefect, dYrelReduction, false, false)));
	}
	void set_Y_check(SmartPtr<IConvergenceCheck<vector_type> > spCC)
	{
		spYConvCheck = spCC;
		for(size_t i=0; i<levels.size(); i++)
			levels[i]->spYConvCheck = spCC->clone();
		m_cycleType = std::string("YCycle (ConvCheck = ") + spCC->config_string() + std::string(")");
	}


	std::string get_cycle_type() const {
		return m_cycleType;
	}

	void set_max_levels(size_t new_max_levels) {
		UG_COND_THROW(new_max_levels == 0, "cannot set max levels " << new_max_levels << ", has to be > 0");
		m_maxLevels = new_max_levels;
	}
	size_t get_max_levels() const 								{ return m_maxLevels; }

	void set_fsmoothing(bool enable) 							{ m_bFSmoothing = enable; }
	bool get_fsmoothing() const 								{ return m_bFSmoothing;	}

	void set_max_nodes_for_base(size_t newMaxNodesForBase) 		{ m_maxNodesForBase = newMaxNodesForBase; }
	size_t get_max_nodes_for_base() const 						{ return m_maxNodesForBase; }

	void set_max_fill_before_base(double newMaxFillBeforeBase)  { m_dMaxFillBeforeBase = newMaxFillBeforeBase;	}
	double get_max_fill_before_base() const 					{ return m_dMaxFillBeforeBase; }

	void set_min_nodes_on_one_processor(size_t newMinNodes) 	{ m_minNodesOnOneProcessor = newMinNodes; }
	size_t get_min_nodes_on_one_processor() const 				{ return m_minNodesOnOneProcessor; }

	void set_preferred_nodes_on_one_processor(size_t i) 		{ m_preferredNodesOnOneProcessor = i; }
	size_t get_preferred_nodes_on_one_processor() const 		{ return m_preferredNodesOnOneProcessor; }

	void set_checkLevel_post_iterations(size_t i) 				{ m_checkLevelPostIterations = i; }
	size_t get_checkLevel_post_iterations() const 				{ return m_checkLevelPostIterations; }

	void set_verbose(int verbose) 								{ m_verbose = verbose; }

	void set_presmoother(SmartPtr<ILinearIterator<vector_type> > presmoother) {
		m_presmoother = presmoother;
	}
	void set_postsmoother(
			SmartPtr<ILinearIterator<vector_type> > postsmoother) {
		m_postsmoother = postsmoother;
	}
	void set_base_solver(
		SmartPtr<ILinearOperatorInverse<vector_type> > basesolver) {
		m_basesolver = basesolver;
	}
	void set_position_provider(SmartPtr<IPositionProvider<2> > ppp2d) {
		m_pPositionProvider2d = ppp2d;
	}

	void set_position_provider(SmartPtr<IPositionProvider<3> > ppp3d) {
		m_pPositionProvider3d = ppp3d;
	}

	void set_positions(const MathVector<2> *pos, size_t size) {
		m_dbgPositions.resize(size);
		for (size_t i = 0; i < size; ++i) {
			m_dbgPositions[i].x() = pos[i].x();
			m_dbgPositions[i].y() = pos[i].y();
			m_dbgPositions[i].z() = 0.0;
		}
		m_dbgDimension = 2;
	}
	void set_positions(const MathVector<3> *pos, size_t size) {
		m_dbgPositions.resize(size);
		for (size_t i = 0; i < size; ++i)
			m_dbgPositions[i] = pos[i];
		m_dbgDimension = 3;
	}

	void set_matrix_write_path(const char *path) {
		m_writeMatrixPath = path;
		m_writeMatrices = true;
	}

	void set_debug_write_all_vectors(bool b) 				{ m_bDebugWriteAllVectors=b; 	}
	void set_nr_of_preiterations_at_check(size_t i) 		{ m_iNrOfPreiterationsCheck = i; }
	void set_preiterations_mimum_defect_at_check(double d) 	{ m_dPreiterationsMimumDefect = d; }
	void	set_galerkin_truncation(double d)				{ m_dGalerkinTruncation = d; }
	double	get_galerkin_truncation() const					{ return m_dGalerkinTruncation; }

	/// if true, do not re-create AMG we preprocess is called
	void set_one_init(bool b) 								{ m_bOneInit = b; }


////// results, output ////////
// public results and output functions
public:
	/// print level informations
	std::string level_information_str() const
	{
		if(m_bInited==false) return "no level information: not initialized.";
		if(m_bEmptyProcess) return "this is an empty process.";
		std::stringstream ss;
		ss.setf(std::ios::fixed);
		ss << "Operator Complexity: " << get_operator_complexity() << ", grid complexity: "
				<< get_grid_complexity() << ".\n";
		ss << "Whole setup took " << get_timing_whole_setup_ms() << " ms, coarse solver setup took "
				<< get_timing_coarse_solver_setup_ms() << " ms.\n";

		/*for(size_t i = 0; i<get_used_levels(); i++)
			ss << get_level_information(i)->tostring() << "\n";*/
		StringTableStream sts;
		AMGLevelInformation::get_header(sts);
		for(size_t i = 0; i<get_used_levels(); i++)
			get_level_information(i)->add_data(sts);
		ss << sts;

		return ss.str();
	}

	std::string level_conv_rates()
	{
		if(is_fixed_conv_check())
			return "not available for fixed conv check.";
		else
		{
			std::stringstream ss;
//			double totalCalls = m_vLvlConvRates[0].calls;
			ss << "Level 0 (cycle)\n";
			for(size_t i = 1; i<get_used_levels()-1; i++)
			{
//				double wSteps = pow(2, i);
				ss << "Level " << i <<
						" avg gamma: " << m_vLvlConvRates[i].totalSteps/((double)m_vLvlConvRates[i].calls) <<
						" avg reduction: " << m_vLvlConvRates[i].avgReduction << " total Steps " << m_vLvlConvRates[i].totalSteps << " calls: " << m_vLvlConvRates[i].calls << "\n";
			}
			ss << "Level " << get_used_levels()-1 << " (basesolver)\n";
			return ss.str();
		}

	}
	void clear_level_conv_rates()
	{
		size_t s = m_vLvlConvRates.size();
		m_vLvlConvRates.clear();
		m_vLvlConvRates.resize(s);
	}

	const AMGLevelInformation *get_level_information(size_t level) const {
		if (level < m_levelInformation.size())
			return &m_levelInformation[level];
		else
			return NULL;
	}

	size_t my_multigrid_levels()
	{
#ifdef UG_PARALLEL
		if(m_agglomerateLevel == (size_t)-1)
			return m_usedLevels-1;
		else
			return m_agglomerateLevel;
#else
		return m_usedLevels-1;
#endif
	}

	/// @return c_A = total nnz of all matrices divided by nnz of matrix A
	double get_operator_complexity() const 					{ return m_dOperatorComplexity; }
	/// @return c_G = total number of nodes of all levels divided by number of nodes on level 0
	double get_grid_complexity() const 						{ return m_dGridComplexity;	}
	/// @return the time spent on the whole setup in ms
	double get_timing_whole_setup_ms() const 				{ return m_dTimingWholeSetupMS;	}
	/// @return the time spent in the coarse solver setup in ms
	double get_timing_coarse_solver_setup_ms() const 		{ return m_dTimingCoarseSolverSetupMS;	}

	/*
	 size_t get_nr_of_coarse(size_t level)
	 {
	 assert(level+1 < m_usedLevels);
	 return A[level+1]->length;
	 }
	 */
	/// returns the number of AMG levels
	size_t get_used_levels() const {
		if(m_bEmptyProcess) return 0;
		return m_totalUsedLevels;
	}

protected:
/////// debug ///////
	/**
	 * reads in the positions from the position provider
	 * @sa set_position_provider, set_positions
	 */
	void update_positions();

	std::string GetAMGFilenamePath(const matrix_type &m, int levelDir, std::string name, int level, std::string extension, bool bWriteHeader=false);
	std::string GetAMGFilename(const matrix_type &m, std::string name, int level, std::string extension, bool bWriteHeader=false);

	void write_debug_matrices(matrix_type &AH, prolongation_matrix_type &R, const matrix_type &A, prolongation_matrix_type &P, size_t level);

	template<typename TMatrix>
	void write_debug_matrix(TMatrix &mat, size_t fromlevel, size_t tolevel, const char *name, bool bAddMarkers=false);

	template<typename TNodeType>
	void write_debug_matrix_markers(size_t level, const TNodeType &nodes);
	/**
	 * writes the vector d in a connection-viewer-vec format
	 *
	 * @param 	filename
	 * @param	d
	 * @param	level
	 */
	template<typename TVector>
	bool writevec(std::string filename, const TVector &d, size_t level,
			const TVector *solution = NULL);
	std::string get_output_directory(int level);

public:
	/// debug output of interfaces written in ConnectionViewer-format
	void write_interfaces();

/////// checks ///////
public:
	bool check2(const vector_type &const_c, const vector_type &const_d);
	bool check_fsmoothing();
	bool check_level(vector_type &c, vector_type &d, matrix_type &A,
			size_t level, checkResult &res, const vector_type *solution = NULL);
//	bool check(IMatrixOperator const vector_type &const_c, const vector_type &const_d);

	bool check(const vector_type &const_c, const vector_type &const_d);
	bool check();
	bool check(vector_type &c, vector_type &d, bool bRandom);

	double check_interpolation(vector_type &v, const char *filename);
	void use_checks(bool b)
	{
		m_bChecks = b;
	}


/////// smoothing ///////
public:
	bool init_smoothers();
protected:
	void init_smoother(size_t level);

	/// init f-smoothing (get consistent matrix)
	void init_fsmoothing();
	/// f-smoothing
	bool f_smoothing(vector_type &corr, vector_type &d, size_t level);

	/**
	 * create fine marks for F-Smoothing
	 * @param level		current AMG level
	 * @param amgnodes 	for getting amgnodes[i].is_fine()
	 * @param N			number of elements in amgnodes with is_fine() == true.
	 */
	template<typename TAMGNodes>
	void create_fine_marks(int level, TAMGNodes &amgnodes, size_t N);

	/**
	 * calculate so that parentIndex[level+1][newIndex[i]] = i
	 * @param level		 current AMG level
	 * @param newIndex	 new indices of the nodes. -1 means node is fine.
	 * @param nrOfCoarse nr of Coarse nodes.
	 */
	void create_parent_index(int level, stdvector<int> newIndex,
			size_t nrOfCoarse);

/////// cycle ///////
protected:
	/**
	 * handles agglomeration, and then calls
	 * add_correction_and_update_defect2 with collected vectors
	 * @param c			result: computed correction is added
	 * @param d			defect, result: updated defect.
	 * @param level		AMG level we operate on
	 * @return			true
	 */
	bool add_correction_and_update_defect(vector_type &c, vector_type &d,
			size_t level = 0);
	/**
	 * adds the correction dc to c and updates the defect d with d -= A dc.
	 * uses the normal recursive AMG cycle for the coarse problem
	 * actual AMG cycle.
	 * @param c			result: computed correction is added
	 * @param d			defect, result: updated defect.
	 * @param A			matrix A to use
	 * @return			true
	 */
	bool add_correction_and_update_defect2(vector_type &c, vector_type &d,
			matrix_operator_type &A, size_t level = 0);

	/**
	 * for the purposes: like a normal AMG cycle, except that it will try to solve
	 * the coarse problem at exactLevel exactly by linear iteration
	 * @param c				result: computed correction is added
	 * @param d				defect, result: updated defect.
	 * @param level			current level
	 * @param exactLevel	level to solve problem exactly
	 * @return				true
	 */
	bool add_correction_and_update_defect(vector_type &c, vector_type &d,
			size_t level, size_t exactLevel);

	bool get_correction(vector_type &c, const vector_type &d);


	/// calculates vH[j] = v[m_parentIndex[level+1][j]];
	bool injection(vector_type &vH, const vector_type &v, size_t level);

	bool is_fixed_conv_check() const
	{
		return spYConvCheck. template cast_dynamic<FixedConvergenceCheck<vector_type> > () .valid();
	}

/////// other ///////
	double Dnorm(vector_type &t, size_t level);
	double Anorm(vector_type &t, size_t level);

	double AnormConst(const vector_type &t, size_t level)
	{
		SmartPtr<vector_type> tmp = t.clone();
		return Anorm(*tmp, level);
	}


private:
/////// agglomeration ///////
#ifdef UG_PARALLEL
	/// starts the agglomeration process
	bool agglomerate(size_t level);

	/**
	 * agglomerate the vectors and matrix from all processors in mergeWith
	 * this processor becomes father for this level
	 * @param merge		this map maps neighbors which got agglomerated to their new cores
	 * @param mergeWith the processors to agglomerate with
	 * @param level
	 */
	bool agglomerate_from(const std::map<int, int> &merge,
			const std::vector<int> &mergeWith, size_t level, ParallelNodes &PN);

	/**
	 * agglomerates this processor on the level to mergeTo
	 * this processor becomes child for this level
	 * @param merge		this map maps neighbors which got agglomerated to their new cores
	 * @param mergeTo 	the processor to merge to
	 * @param level		the AMG level0
	 */
	bool agglomerate_to(const std::map<int, int> &merge, int mergeTo,
			size_t level, ParallelNodes &PN);
	/**
	 * even though this processor does not get agglomerated and does not agglomerate,
	 * he might need to change his interfaces because neighbors could have gotten
	 * agglomerated
	 * @param merge		this map maps neighbors which got agglomerated to their new cores
	 * @param level		the AMG level0
	 */
	bool agglomerate_adjust_interfaces(const std::map<int, int> &merge,
			size_t level, ParallelNodes &PN);
#endif

/////// DATA ///////
protected:

	// other options
	double m_dGalerkinTruncation;

/////// agglomeration ///////

	size_t m_maxNodesForBase;///< max nr of coarse nodes before Base solver is used
	double m_dMaxFillBeforeBase;///< max fill rate before Base solver is used
	size_t m_minNodesOnOneProcessor;///< min nr of nodes on one processor (for agglomeration)
	size_t m_preferredNodesOnOneProcessor;///< preferred nr of nodes on one processor (for agglomeration)
	bool m_bUseCollectedSolver;



/////// checks ///////

	size_t m_iNrOfPreiterationsCheck;///< nr of mg cycles performed before checking
	double m_dPreiterationsMimumDefect;	///< minimum defect for preiterations before checking
	size_t m_checkLevelPostIterations;


	vector_type *m_vec4;	///< temporary Vector for defect (in get_correction)
	int iteration_glboal;


/////// smoothing ///////

	size_t m_numPreSmooth;				///< nu_1 : nr. of pre-smoothing steps
	size_t m_numPostSmooth;				///< nu_2: nr. of post-smoothing steps
	SmartPtr<ILinearIterator<vector_type> > m_presmoother;///< presmoother template
	SmartPtr<ILinearIterator<vector_type> > m_postsmoother;	///< postsmoother template \note: may be pre=post, is optimized away.
	bool m_bFSmoothing;						///< if true, use F-Smoothing

/////// base ///////

	size_t m_maxLevels;						///< max. nr of levels used for FAMG
	SmartPtr<ILinearOperatorInverse<vector_type> > m_basesolver; ///< the base solver


/////// cycles ///////

	std::string m_cycleType;			///< gamma: cycle type (1 = V-Cycle, 2 = W-Cycle)
	SmartPtr<IConvergenceCheck<vector_type> > spYConvCheck;

	void calc_level_conv_rates(SmartPtr<IConvergenceCheck<vector_type> > pCheck, size_t level);
	std::vector<AMGLevelConvergenceRates> m_vLvlConvRates;



	bool m_bInited;	///< true if inited. needed since preprocess doesnt give us a InterfaceCommunicator atm.

/////// results ///////
	double m_dOperatorComplexity;
	double m_dGridComplexity;
	double m_dTimingWholeSetupMS;
	double m_dTimingCoarseSolverSetupMS;

	size_t m_usedLevels;						///< nr of FAMG levels used
	size_t m_totalUsedLevels;

	stdvector<SmartPtr<AMGLevelType> > levels;

/////// debug ///////
	bool m_writeMatrices;
	std::string m_writeMatrixPath;
	SmartPtr<IPositionProvider<2> > m_pPositionProvider2d;
	SmartPtr<IPositionProvider<3> > m_pPositionProvider3d;
	cAMG_helper m_amghelper;///< helper struct for viewing matrices (optional)
	stdvector<MathVector<3> > m_dbgPositions;///< positions of geometric grid (optional)
	int m_dbgDimension;				///< dimension of geometric grid (optional)
	int m_verbose;							///< verbose level. 0 = no output.
	stdvector<stdvector<int> > m_parentIndex;///< parentIndex[i] is the index of i on the finer level


/////// misc ///////
	stdvector<AMGLevelInformation> m_levelInformation;
	void calculate_level_information(size_t level, double createAMGlevelTiming);


	bool m_bOneInit;///< if true, do not re-create AMG we preprocess is called

	bool m_bDebugWriteAllVectors;
	bool m_bChecks;

	/// always false in serial
	bool m_bEmptyProcess;

/////// agglomeration ///////

#ifdef UG_PARALLEL
	// on level m_agglomerateLevel
	IndexLayout agglomerateSlaveLayout;
	size_t m_agglomerateLevel;
	size_t m_baseLevel;

	// stuff from old level 0 agglomeration
	pcl::ProcessCommunicator m_emptyPC;
	IndexLayout m_emptyLayout;

	/// return true if we have been merged to another core on this level
	bool isMergingSlave(size_t level);
	/// return true if we are the master for agglomeration on this level
	bool isMergingMaster(size_t level);
	/// return true if we are neither horizontal master nor slave on this level
	bool isNotMerging(size_t level);

	bool isMasterBaseLevel(size_t level)
	{
		return level == m_usedLevels-1 && isMergingSlave(level) == false;
	}
#endif

public:
	void check_vector(SmartPtr<ITestvectorCalculator<TAlgebra> > testvectorCalculator, std::string desc);

protected:

};

// end group plugins_amg
/// \}

}// namespace ug

#include "amg_base_impl.h"

#endif // __H__UG__LIB_ALGEBRA__AMG_SOLVER__AMG_BASE_H__
