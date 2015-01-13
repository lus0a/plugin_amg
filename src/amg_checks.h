/**
 * \file amg_checks.h
 * \author Martin Rupp
 *
 * \date 03.08.11
 *
 * Goethe-Center for Scientific Computing 2011
 *
 *
 * in this file are functions for checking the convergence behaviour on the AMG levels (for bugfixing etc.)
 */
#include <sstream>
#include <sys/stat.h>
#include "common/util/file_util.h"
#include "common/util/stringify.h"
#include "lib_algebra/operator/energy_convergence_check.h"
namespace ug
{

/// \addtogroup plugins_amg
/// \{

#ifdef UG_PARALLEL
template<typename TAlgebra>
void AMGBase<TAlgebra>::write_interfaces()
{
	try{
	AMG_PROFILE_FUNC();
	if(!m_writeMatrices)
	{
		UG_LOG("no position data given.\n");
		return;
	}
	cAMG_helper &h = m_amghelper;
	for(size_t level=0; level<m_usedLevels; level++)
	{
		std::string filename =
				Stringify() << get_output_directory(level) << "AMG_interface_" << pcl::ProcRank() << ".mat";
		std::fstream file(filename.c_str(), std::ios::out);

		ConnectionViewer::WriteGridHeader(file, h.positions[level], h.dimension);

		const IndexLayout &masterLayout = levels[level]->pA->layouts()->master();
		const IndexLayout &slaveLayout = levels[level]->pA->layouts()->slave();
		std::vector<std::string> v;
		v.resize(h.positions[level].size());


		//UG_LOG("==============================\n LEVEL "<< level << "\n============================\n")
		//PRINTLAYOUT(levels[level]->com, levels[level]->masterLayout, levels[level]->slaveLayout);
		//PRINTLAYOUT(levels[level]->com, masterLayout, slaveLayout);

		for(IndexLayout::const_iterator iter = masterLayout.begin(); iter != masterLayout.end(); ++iter)
		{
			const IndexLayout::Interface &interface = masterLayout.interface(iter);
			int pid = masterLayout.proc_id(iter);
			for(IndexLayout::Interface::const_iterator iter = interface.begin(); iter != interface.end(); ++iter)
			{
				size_t i = interface.get_element(iter);
				v[i].append(Stringify() << "M" << pid << " ");
			}
		}

		for(IndexLayout::const_iterator iter = slaveLayout.begin(); iter != slaveLayout.end(); ++iter)
		{
			const IndexLayout::Interface &interface = slaveLayout.interface(iter);
			int pid = slaveLayout.proc_id(iter);
			for(IndexLayout::Interface::const_iterator iter = interface.begin(); iter != interface.end(); ++iter)
			{
				size_t i = interface.get_element(iter);
				v[i].append(Stringify() << "S" << pid << " ");
			}
		}

		for(size_t i=0;i<v.size(); i++)
		{
			std::string &s = v[i];
			if(s.size() > 0)
				file << i << " " << i << " " << s << std::endl;
		}
	}
	}UG_CATCH_THROW_FUNC();
}
#else
template<typename TAlgebra>
void AMGBase<TAlgebra>::write_interfaces()
{
	return;
}
#endif


template<typename TAlgebra>
bool AMGBase<TAlgebra>::check()
{
	vector_type c, d;
//	UG_LOG_ALL_PROCS("check\n");

	CloneVector(c, levels[0]->corr);
	CloneVector(d, levels[0]->corr);
	return check(c, d, true);

}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// check
//--------------------------------
/**
 * this function checks convergence level for level by using \sa check_level
 *
 * \param 	const_c				c
 * \param	const_d				d
 */
template<typename TAlgebra>
bool AMGBase<TAlgebra>::check(vector_type &c, vector_type &d, bool bRandom)
{
	try{
	int oldOutputProcess = GetLogAssistant().get_output_process();
	//UG_LOG("\nLEVEL 0\n\n");
#ifdef UG_PARALLEL
	if(m_agglomerateLevel == (size_t)-1)
		{UG_LOG("Number of levels is " << m_usedLevels << " (no agglomeration)\n");}
	else
		{UG_LOG("agglomerate Level is " << m_agglomerateLevel << ", used levels is " << m_usedLevels << "\n");}
#endif

	std::vector<checkResult> checkRes(m_usedLevels-1);
	for(size_t level=0; level<m_usedLevels-1; level++)
	{
//		UG_LOG_ALL_PROCS((Stringify() << "level " << level << "\n").str());

		vector_type *pC, *pD;

		if(level==0)
		{
			pC = &c;
			pD = &d;
		}
		else
		{
			pC = &levels[level-1]->cH;
			pD = &levels[level-1]->dH;
			pC->set(0.0); // think about that
		}

		SmartPtr<matrix_operator_type> pA = levels[level]->pA;

#ifdef UG_PARALLEL
		if(m_agglomerateLevel==(unsigned int)(-1) && level == m_usedLevels-1)
		{
//			UG_LOG_ALL_PROCS("breaking2 " << level << "\n");
			break;
		}
		if(levels[level]->bLevelHasMergers)
		{
			UG_LOG_ALL_PROCS("level " << level << " has mergers\n");
			if(!bRandom)
			{
				gather_vertical(*pD, levels[level]->collD, level, PST_ADDITIVE);
				gather_vertical(*pC, levels[level]->collC, level, PST_CONSISTENT);
			}
			if(isMergingSlave(level))
			{
//				UG_LOG_ALL_PROCS("breaking, because i'm slave on level " << level << "\n");
				break;
			}
			pD = &levels[level]->collD;
			pC = &levels[level]->collC;
			pA = levels[level]->collectedA;

			const pcl::ProcessCommunicator &pc = levels[level]->pAgglomeratedA->layouts()->proc_comm();
			if(oldOutputProcess == 0 && GetLogAssistant().get_output_process() != pcl::ProcRank()
				&& pc.get_proc_id(0) == pcl::ProcRank())
			{
				GetLogAssistant().set_output_process(pcl::ProcRank());
//				UG_LOG_ALL_PROCS("Processor " << pcl::ProcRank() << " is now output process.\n");
			}
		}
#endif

		vector_type c2, d2;
		CloneVector(d2, *pD);
		CloneVector(c2, *pC);
		d2 = *pD;
		c2 = *pC;

		UG_LOG("--[[ CHECKING LEVEL" << level << " ]]--\n");
		if(bRandom)
		{
			// Ac = 0
			d2.set(0.0);
			c2.set_random(-1.0, 1.0);
			// d = 0 - Ac
			pA->matmul_minus(d2, c2);
#ifdef UG_PARALLEL
			c2.change_storage_type(PST_CONSISTENT);
			d2.change_storage_type(PST_ADDITIVE);
#endif
		}
		check_level(c2, d2, *pA, level, checkRes[level]);


		if(!bRandom)
			levels[level]->R.apply(levels[level]->dH, d2);
	}

	UG_LOG("out\n");
#ifdef UG_PARALLEL
//	const_c.layouts()->proc_comm().barrier();
#endif
	GetLogAssistant().set_output_process(oldOutputProcess);
	}UG_CATCH_THROW_FUNC();
	return true;
}



///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// check
//--------------------------------
/**
 * this function checks convergence level for level by using \sa check_level
 *
 * \param 	const_c				c
 * \param	const_d				d
 */
template<typename TAlgebra>
bool AMGBase<TAlgebra>::check(const vector_type &const_c, const vector_type &const_d)
{
	AMG_PROFILE_FUNC();
	UG_LOG("\n");
	UG_LOG("         check     \n");
	UG_LOG("========================\n");

	/*if(m_usedLevels <= 1)
	{
		UG_LOG("No Multigrid hierachy.\n");
		return true;
	}*/
	vector_type c, d;
	CloneVector(c, const_c);
	CloneVector(d, const_d);
	c = const_c;
	d = const_d;
	return check(c, d, false);
}




template<typename TVector>
double ConstTwoNorm(const TVector &const_d)
{
	TVector d;
	CloneVector(d, const_d);
	d = const_d;
	return d.norm();

}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// check_level:
//--------------------------------
/**
 * this function checks convergence on a specified level.
 * it will log all convergence rates:
 * - presmoother
 * - fsmoother
 * - coarse correction
 * - f-postsmoother
 * - postsmoother
 *
 * as coarse correction we use the normal multigrid on a level below as often as necessary
 * to reduce the error on the coarse grid by 0.01 to distinguish between bad convergence on this level
 * and bad convergence on other levels
 *
 * \param 	c
 * \param	d
 * \param	level
 */
template<typename TAlgebra>
bool AMGBase<TAlgebra>::check_level(vector_type &c, vector_type &d, matrix_type &A,
		size_t level, checkResult &res, const vector_type *solution)
{
	size_t line;
	try{
		d.layouts()->proc_comm().barrier();
	AMG_PROFILE_FUNC();
	AMGLevelType &L = *levels[level];

	UG_ASSERT(c.size() == d.size() && c.size() == A.num_rows(),
			"c.size = " << c.size() << ", d.size = " << d.size() << ", A.size = " << A.num_rows() << ": not matching");

#ifdef UG_PARALLEL
	if(!d.has_storage_type(PST_ADDITIVE) || !c.has_storage_type(PST_CONSISTENT))
	{
		UG_LOG("ERROR: In 'amg::check':Inadequate storage format of Vectors.\n");
		return false;
	}
#endif

	line=__LINE__;
	std::stringstream reason;
	//UG_LOG("preprenorm: " << d.norm() << std::endl);
	vector_type &corr = levels[level]->corr;

	double firstnorm = ConstTwoNorm(d), n2=firstnorm, n1=n2;
	//UG_LOG("firstnorm " <<  firstnorm << "\n");
	double prevrate=0;
	if(m_writeMatrices) writevec("cc0", c, level, solution);
	if(m_writeMatrices) writevec("c0", c, level);
	if(m_writeMatrices) writevec("d0", d, level);

	UG_LOG("\tType    " << "  Iter#     Defect         Rate\n");
	UG_LOG("\t                 " << std::scientific << firstnorm << "\n");


	// V/W-check
	{
		vector_type c2;
		CloneVector(c2, c);
		c2 = c;
		vector_type d2;
		CloneVector(d2, d);
		d2 = d;

		line=__LINE__;
		for(size_t i=0; i<m_iNrOfPreiterationsCheck; i++)
		{
			line=__LINE__;
			add_correction_and_update_defect2(c2, d2, *L.pAgglomeratedA, level);
			line=__LINE__;
			n1 = n2;
			n2 = ConstTwoNorm(d2);

			const double deltaRate=0.005;
			double rate = n2/n1;

			UG_LOG("\tpreIter " << std::setw(4) << i+1 << ":    " << std::scientific << n2 <<
								"    " << std::scientific << rate);
			if((prevrate/rate < (1+deltaRate) && prevrate/rate > (1-deltaRate)) )
			{
				UG_LOG(" (stable)\n");
				break;
			}

			UG_LOG("\n");
			prevrate = rate;


		}
	}

	prevrate = 0;
	n2=firstnorm;
	n1=n2;

	SmartPtr<IConvergenceCheck<vector_type> > pCheckOld = levels[level]->spYConvCheck;
	levels[level]->spYConvCheck = make_sp(new StdConvCheck<vector_type>(m_iNrOfPreiterationsCheck, m_dPreiterationsMimumDefect, 1e-10, false, true));

	line=__LINE__;
	// 2-grid check
	for(size_t i=0; i<m_iNrOfPreiterationsCheck; i++)
	{
		line=__LINE__;
		add_correction_and_update_defect2(c, d, *L.pAgglomeratedA, level);
		line=__LINE__;

		if(m_writeMatrices) writevec(Stringify() << "cc" << i+1, c, level, solution);
		if(m_writeMatrices) writevec(Stringify() << "c" << i+1, c, level);
		if(m_writeMatrices) writevec(Stringify() << "d" << i+1, d, level);

		n1 = n2;
		n2 = ConstTwoNorm(d);

		const double deltaRate=0.005;
		double rate = n2/n1;
		d *= 1/n2;
		c *= 1/n2;
		n2 = 1.0;

		UG_LOG("\t2-preIter " << std::setw(4) << i+1 << ":    " << std::scientific << "          " <<
				"    " << std::scientific << rate);

		if((prevrate/rate < (1+deltaRate) && prevrate/rate > (1-deltaRate)) )
		{
			UG_LOG(" (stable)\n");
			break;
		}
		UG_LOG("\n");

		prevrate = rate;

/*#ifdef UG_PARALLEL
		if(level+1 < levels.size() && !levels[level]->bLevelHasMergers)
#else
			if(level+1 < levels.size())
#endif*/
		if(pcl::NumProcs() == 1)
		{
			vector_type vH;
			CloneVector(vH, L.cH);
			vector_type v2;
			CloneVector(v2, c);

			vector_type *pC = &c;
			L.R.apply(vH, c);
			if(m_writeMatrices) writevec(std::string("Rc")+ToString(i+1), vH, level+1);
			L.P.apply(v2, vH);
			if(m_writeMatrices) writevec(std::string("PRc")+ToString(i+1), v2, level);
			if(m_writeMatrices) writevec(std::string("(1-PR)c")+ToString(i+1), v2, level, pC);

			injection(vH, c, level);
			if(m_writeMatrices) writevec(std::string("Rinj_c")+ToString(i+1), vH, level+1);
			L.P.apply(v2, vH);
			if(m_writeMatrices) writevec(std::string("PRinj_c")+ToString(i+1), v2, level);
			if(m_writeMatrices) writevec(std::string("(1-PRinj)c")+ToString(i+1), v2, level, pC);
		}
	}


	line=__LINE__;
	corr.set(0.0);

	double n1A = AnormConst(d, level);

	for(size_t iterations=0; iterations<1; iterations++)
	{
		UG_LOG("\n");

			double prenorm = ConstTwoNorm(d);
				UG_LOG("\tdefect below is in relation to " << std::scientific << prenorm << ".\n");
				UG_LOG("\tType    " << "  Iter#     Defect         Rate  \n");
				//UG_LOG("Prenorm = " << prenorm << "\n");

				n1 = ConstTwoNorm(d);
				n2=n1;

				if(n1 < 1e-13) reason << "- error too small\n";

		//c.set(0.0);
	#ifdef UG_PARALLEL
		corr.set_storage_type(PST_CONSISTENT);
	#endif

		if(m_writeMatrices) writevec(ToString(iterations)+"S0_d_before", d, level);
		if(m_writeMatrices) writevec(ToString(iterations)+"S0_c_before", c, level, solution);

		// presmooth
		for(size_t i=0; i < m_numPreSmooth; i++)
		{
			L.presmoother->apply_update_defect(corr, d);
			c += corr;
			n2 = ConstTwoNorm(d);
			UG_LOG("\tpreSm   " << std::setw(4) << i+1 << ":    " << std::scientific << n2 <<
							"    " << std::scientific << n2/n1 << "\n");
			n1 = n2;
		}
		line=__LINE__;
		res.preSmoothing = n2/prenorm;

		if(m_writeMatrices) writevec(ToString(iterations)+"S1_d_presmoothed", d, level);
		if(m_writeMatrices) writevec(ToString(iterations)+"S1_c_presmoothed", c, level, solution);

		if(m_bFSmoothing)
		{
			f_smoothing(corr, d, level);
			c+=corr;
			n2 = ConstTwoNorm(d);
			UG_LOG("\tpreFSm  " << std::setw(4) << 0 << ":    " << std::scientific << n2 <<
									"    " << std::scientific << n2/n1 << "\n");
			n1 = n2;
		}
		line=__LINE__;
		res.preFSmoothing = n2/prenorm;

		if(m_writeMatrices) writevec(ToString(iterations)+"S1_d_fpresmoothed", d, level);
		if(m_writeMatrices) writevec(ToString(iterations)+"S1_c_fpresmoothed", c, level, solution);

		vector_type &cH = levels[level]->cH;
		vector_type &dH = levels[level]->dH;
		cH.set(0.0);
		dH.set(0.0);
	#ifdef UG_PARALLEL
		cH.set_storage_type(PST_CONSISTENT);
	#endif
		// restrict defect
		// dH = m_R[level]*d;
		L.R.apply(dH, d);

		if(m_writeMatrices) writevec(ToString(iterations)+"S2_dH", dH, level+1);
		double nH1 = ConstTwoNorm(dH);


		//double preHnorm=nH1;
		size_t i=0;
		line=__LINE__;

		{
			//UG_LOG_ALL_PROCS(__LINE__);
			cH.set(0.0);
			double pre_nH=nH1;
			double nH2=nH1;
			SmartPtr<IConvergenceCheck<vector_type> > pCheckOld = levels[level+1]->spYConvCheck;
			levels[level+1]->spYConvCheck = make_sp(new StdConvCheck<vector_type>(100, 1e-20, 1e-10, false, false));
			//levels[level+1]->spYConvCheck = make_sp(new StdConvCheck<vector_type>(100, 1e-20, 1e-10, true, true));

#if 0
			SmartPtr<IConvergenceCheck<vector_type> > pCheck = make_sp(new StdConvCheck<vector_type>(100, 1e-20, 1e-10, true, true));
			pCheck->set_offset(level);
			pCheck->set_name((Stringify() << "Cycle " << level).str());
		  //SmartPtr<IConvergenceCheck<vector_type> > pCheck(new FixedConvergenceCheck<vector_type>(1));
			//FixedConvergenceCheck<vector_type> pCheck(m_iYCycle);
			pCheck->start(dH);

			//UG_LOG("\n" <<repeat(' ', level+1) << "." <<"\n");
			do
			{
				add_correction_and_update_defect(cH, dH, level+1);
				pCheck->update(dH);
			}
			while(!pCheck->iteration_ended());


#else
			for(i=0; i<100; i++)
			{
				add_correction_and_update_defect(cH, dH, level+1);
				nH1 = nH2;
				nH2 = ConstTwoNorm(dH);
				//UG_LOG("coarse correction (on coarse) " << i+1 << ": " << nH2/nH1 << " " << nH2 << "\n");
				/*if(i == 0)
				{	UG_LOG("coarse correction (on coarse) " << i+1 << ": " << nH2/nH1 << " " << nH2/preHnorm << "\n"); }
				else UG_LOG(".");*/
				if(m_writeMatrices) writevec(ToString(iterations)+(std::string("AMG_V")+ToString(i)+std::string("_cH_L")).c_str(), cH, level+1);
				if(m_writeMatrices) writevec(ToString(iterations)+(std::string("AMG_V")+ToString(i)+std::string("_dH_L")).c_str(), dH, level+1);
				if(nH2/pre_nH < 1e-12 || nH2 < 1e-12)
					break;
			}
#endif

			levels[level+1]->spYConvCheck = pCheckOld;

			res.iInnerIterations = i+1;
			res.coarseDefect = nH2;
			res.lastCoarseReduction = nH2/nH1;

			UG_LOG("\tCoarse Correction (on coarse):\n");
			UG_LOG("\t  Norm of dH before = " << pre_nH << ", dH now = " << res.coarseDefect << "\n");
			UG_LOG("\t  Took " << res.iInnerIterations << " cycles to reduce by " << res.coarseDefect/pre_nH << "\n");
			UG_LOG("\t  last reduction was " << res.lastCoarseReduction <<  "\n");


			/*cH.set(0.0);
			for(i=0; i<100; i++)
			{
				for(int j=0; j< m_cycleType; j++)
					add_correction_and_update_defect(cH, dH, level+1);

				double nH2 = ConstTwoNorm(dH);
				if(i < 6)
				{	UG_LOG("\tcoarse correction (on coarse) " << i+1 << ": " << nH2/nH1 << "\n"); nH1 = nH2; }
				if(nH2/preHnorm < 0.01) { UG_LOG("coarse solver reduced by 0.01 in iteration " << i+1 << std::endl); break; }
			}*/
		}

		line=__LINE__;
		// interpolate correction
		// corr = m_P[level]*cH
		if(m_writeMatrices) writevec(ToString(iterations)+"S3_dH_solved", dH, level+1);
		line=__LINE__;
				if(m_writeMatrices) writevec(ToString(iterations)+"S3_cH_solved", cH, level+1);
				line=__LINE__;

				L.P.apply(corr, cH);
				line=__LINE__;

		if(m_writeMatrices) writevec(ToString(iterations)+"S3_corr_on_fine", corr, level);

	#ifdef UG_PARALLEL
		cH.set_storage_type(PST_CONSISTENT);
	#endif

		// add coarse grid correction to level correction
		line=__LINE__;
		c += corr;

		//update defect
		// d = d - Ah*corr
		A.matmul_minus(d, corr);
		line=__LINE__;
		n2 = ConstTwoNorm(d);
		UG_LOG("\tcoarse correct : " << n2/prenorm << "    " <<n2/n1 << "\n");	n1 = n2;
		res.coarseCorrection = n2/prenorm;

		if(m_writeMatrices) writevec(ToString(iterations)+"S4_d_coarseCorrected", d, level);
		if(m_writeMatrices) writevec(ToString(iterations)+"S4_c_coarseCorrected", c, level, solution);


		// post f-smoothing
		if(m_bFSmoothing)
		{
			f_smoothing(corr, d, level);
			c+=corr;
		}
		n2 = ConstTwoNorm(d);
		line=__LINE__;
		UG_LOG("\tpostFSm " << std::setw(4) << 0 << ":    " << std::scientific << n2 <<
										"    " << std::scientific << n2/n1 << "\n");
		res.postFSmoothing = n2/prenorm;
		line=__LINE__;
		if(m_writeMatrices) writevec(ToString(iterations)+"S5_d_fpost", d, level);
		if(m_writeMatrices) writevec(ToString(iterations)+"S5_c_fpost", c, level, solution);

		// postsmooth
		for(size_t i=0; i < m_numPostSmooth; i++)
		{
			L.postsmoother->apply_update_defect(corr, d);
			c += corr;
			n1 = n2;
			n2 = ConstTwoNorm(d);
			UG_LOG("\tpostSm  " << std::setw(4) << i+1 << ":    " << std::scientific << n2 <<
									"    " << std::scientific << n2/n1 << "\n");
		}

		res.postSmoothing = n2/prenorm;

		UG_LOG("Total Correction in A-Norm: " << AnormConst(d, level)/n1A << "\n");
		line=__LINE__;
		if(res.postSmoothing > 0.1)
		{
			UG_LOG("\n ERROR: Level " << level << " not working: total correction is " << res.postSmoothing << "\n");

			if(res.preSmoothing > 0.99)
			{	UG_LOG(" - presmoother not working (reduction rate " << res.preSmoothing << ")\n"); }
			if(res.lastCoarseReduction > 0.2 || res.coarseDefect > 1e-10)
			{
				UG_LOG(" - lower level not working (done " << res.iInnerIterations << " iterations):\n");
				if(res.lastCoarseReduction > 0.2) { UG_LOG("  - bad convergence rate (" << res.lastCoarseReduction << ")\n"); }
				if(res.coarseDefect > 1e-8) { UG_LOG("  - coarse reduction not sufficient (" << res.coarseDefect << ")\n"); }
			}
		}


		double postnorm = ConstTwoNorm(d);
		UG_LOG(" Level " << level << " reduction: " << postnorm/prenorm << std::endl);
		line=__LINE__;
		if(m_writeMatrices) writevec(ToString(iterations)+"S5_d_postsmoothed", d, level);
		if(m_writeMatrices) writevec(ToString(iterations)+"S5_c_postsmoothed", c, level, solution);

		UG_LOG(" Postnorm = " << postnorm << "\n\n");
		res.reduction = postnorm/prenorm;
	}

	levels[level]->spYConvCheck = pCheckOld;

	}UG_CATCH_THROW(PRETTY_FUNCTION << ", level = " << level << ", failed at line " << line );
	//return postnorm/prenorm;
	return 0;
}


///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// check2
//--------------------------------
/**
 * for all i= 1 .. usedLevels-2
 * 	we check the convergence behaviour of the multigrid level if, on level i, we use an exact solver
 * \param 	const_c				c
 * \param	const_d				d
 */
template<typename TAlgebra>
bool AMGBase<TAlgebra>::check2(const vector_type &const_c, const vector_type &const_d)
{
	try{
	AMG_PROFILE_FUNC();
	UG_LOG("\n");
	UG_LOG("            check2\n");
	UG_LOG("==================================\n");
	UG_LOG("we check the convergence behaviour of the multigrid level if, on level i, we use an exact solver");

	if(m_usedLevels <= 1)
	{
		UG_LOG("No Multigrid hierachy.\n");
		return true;
	}

	vector_type c, d;
	CloneVector(c, const_c);
	CloneVector(d, const_d);

	for(size_t i=1; i<m_usedLevels-1; i++)
	{
		c = const_c;
		d = const_d;
		UG_LOG("Testing Level 0 to i=" << i << ":\n");
		double prenorm = ConstTwoNorm(d);
		add_correction_and_update_defect(c, d, 0, i);
		double postnorm = ConstTwoNorm(d);
		UG_LOG("Reduction is " << postnorm/prenorm << "\n\n");
	}

	UG_LOG("finished check2.\n");

	}UG_CATCH_THROW_FUNC();
	return true;
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// add_correction_and_update_defect
//-----------------------------------
/**
 * \sa add_correction_and_update_defect
 * for all i= 1 .. usedLevels-2
 * 	we check the convergence behaviour of the multigrid level if, on level i, we use an exact solver
 * \param c
 * \param d
 * \param level
 * \param exactLevel	level on which we will use an exact solver, that is: doing MG on that level until norm(d) < 1e-12.
 */
template<typename TAlgebra>
bool AMGBase<TAlgebra>::add_correction_and_update_defect(vector_type &c, vector_type &d, size_t level, size_t exactLevel)
{
	try{
	AMG_PROFILE_FUNC();
	AMGLevelType &L = *levels[level];
	const matrix_operator_type &A = *L.pA;
	UG_ASSERT(c.size() == d.size() && c.size() == A.num_rows(),
			"c.size = " << c.size() << ", d.size = " << d.size() << ", A.size = " << A.num_rows() << ": not matching");

#ifdef UG_PARALLEL
	if(!d.has_storage_type(PST_ADDITIVE) || !c.has_storage_type(PST_CONSISTENT))
	{
		UG_LOG("ERROR: In 'amg::check':Inadequate storage format of Vectors.\n");
		return false;
	}
#endif


	vector_type &corr = levels[level]->corr;
#ifdef UG_PARALLEL
	corr.set_storage_type(PST_CONSISTENT);
#endif

	// presmooth
	for(size_t i=0; i < m_numPreSmooth; i++)
	{
		L.presmoother->apply_update_defect(corr, d);
		c += corr;
	}

	// pre f-smoothing
	if(m_bFSmoothing)
	{
		f_smoothing(corr, d, level);
		c+=corr;
	}

	vector_type &cH = levels[level]->cH;
	vector_type &dH = levels[level]->dH;
	cH.set(0.0);
	dH.set(0.0);
#ifdef UG_PARALLEL
	cH.set_storage_type(PST_CONSISTENT);
#endif

	L.R.apply(dH, d);

	if(level+1 == m_usedLevels-1)
		solve_on_base(cH, dH, level+1);
	else
	{
		cH.set(0.0);
		if(level+1 < exactLevel)
		{
			SmartPtr<IConvergenceCheck<vector_type> > pCheck = levels[level]->spYConvCheck;
			pCheck->start(dH);
			while(!pCheck->iteration_ended())
			{
				add_correction_and_update_defect(cH, dH, level+1, exactLevel);
				pCheck->update(dH);
			}
		}
		else
		{
			UG_LOG("using normal MG cycle from level " << level+1 << "on: ");
			int k;

			for(k=0; k<50; k++)
			{
				add_correction_and_update_defect(cH, dH, level+1);
				UG_LOG(".");
				if(ConstTwoNorm(dH) < 1e-12)
					break;
			}
			if(ConstTwoNorm(dH) > 1e-12) { UG_LOG("not converged after " << k << " steps\n"); }
			else { UG_LOG("took " << k << " steps\n") }
		}
	}

	L.P.apply(corr, cH);
	c += corr;
	A.matmul_minus(d, corr);

	// post f-smoothing
	if(m_bFSmoothing)
	{
		f_smoothing(corr, d, level);
		c+=corr;
	}

	// postsmooth
	for(size_t i=0; i < m_numPostSmooth; i++)
	{
		L.postsmoother->apply_update_defect(corr, d);
		c += corr;
	}
	}UG_CATCH_THROW_FUNC();
	return true;
}


template<typename TAlgebra>
double AMGBase<TAlgebra>::check_interpolation(vector_type &v, const char *filename)
{
	try{
	double anorm = Anorm(v, 0);
	double dnorm = Dnorm(v, 0);
	AMGLevelType &L = *levels[0];
	vector_type &vH = L.cH;
	vector_type &v2 = L.corr;
	//L.R.apply(vH, v);
	injection(vH, v, 0);
	L.P.apply(v2, vH);
	if(m_writeMatrices)
	{
		writevec(std::string(filename)+std::string("_v.vec"), v, 0);
		writevec(std::string(filename)+std::string("_Rv.vec"), vH, 1);
		writevec(std::string(filename)+std::string("_PRv.vec"), v2, 0);
		v -= v2;
		writevec(std::string(filename)+std::string("_(1-PR)v.vec"), v, 0);
	}
	else v -= v2;
	double dnorm2 = Dnorm(v, 0);
	UG_LOG("Anorm(v) = " << anorm << " Dnorm(v) = " << dnorm << " Dnorm( (1-PR)v ) = " << dnorm2 << "\n");
	return dnorm2;
	}UG_CATCH_THROW_FUNC();
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// check_fsmoothing
//--------------------------------
/**
 */
template<typename TAlgebra>
bool AMGBase<TAlgebra>::check_fsmoothing()
{
	try{
	AMG_PROFILE_FUNC();
	UG_LOG("\n");
	UG_LOG("            check_fsmoothing\n");
	UG_LOG("==========================================\n")
	UG_LOG("check the convergence behaviour of a pure f-smoother")

	if(m_usedLevels == 0) return false;

	for(size_t level=0; level<m_usedLevels-1; level++)
	{
		UG_LOG("LEVEL " << level << "\n");
		vector_type c, d;
		CloneVector(c, levels[level]->corr);
		CloneVector(d, levels[level]->corr);

		for(size_t j=0; j<c.size(); j++)
		{
			if(levels[level]->is_fine[j])
				c[j] = urand<double>(-1.0, 1.0);
			else
				c[j] = 0.0;
		}

#ifdef UG_PARALLEL
		c.set_storage_type(PST_ADDITIVE);
		c.change_storage_type(PST_CONSISTENT);
#endif
		// d = A*c
		levels[level]->pA->apply(d, c);
		c.set(0.0);
		double n=ConstTwoNorm(d);

		for(size_t j=0; j<10; j++)
		{
			f_smoothing(c, d, level);
			double n2 = ConstTwoNorm(d);
			UG_LOG(j << ": " << n2/n << "\n");
			n = n2;
			if(m_writeMatrices) writevec(std::string("fsmoothing_c_") + ToString(j), c, level);
			if(m_writeMatrices) writevec(std::string("fsmoothing_d_") + ToString(j), d, level);
		}
	}

	UG_LOG("finished check_fsmoothing.\n");
	}UG_CATCH_THROW_FUNC();
	return true;
}

// end group plugins_amg
/// \}

}
