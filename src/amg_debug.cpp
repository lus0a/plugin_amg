#include "common/log.h"
#ifdef UG_PARALLEL
#include "pcl/pcl_layout_tests.h"
#include "lib_algebra/parallelization/parallel_index_layout.h"
#endif

namespace ug
{

DebugID RSAMGCoarseningDebug("RSAMG.Coarsening");
DebugID RSAMGProlongationDebug("RSAMG.Prolongation");
DebugID AMGAgglomerationDebug("AMG.Agglomeration");

#ifdef UG_PARALLEL
void MyPrintLayout(const pcl::ProcessCommunicator &pc,
				   pcl::InterfaceCommunicator<IndexLayout> &communicator,
				   const IndexLayout &layout1, const IndexLayout &layout2,
				   const char *name1, const char *name2)
{
	UG_LOG("\n========================================\n");
	UG_LOG(name1<<"\n");
	PrintLayout(layout1);
	UG_LOG("\n"<<name2<<"\n");
	PrintLayout(layout2);
	UG_LOG("\n"<<name1 << "-" << name2 <<"\n");
	PrintLayout(pc, communicator, layout1, layout2);
	UG_LOG("\n\n");
}
#endif


}
