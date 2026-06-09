/* nmesh_evolve.h */
/* (c) Wolfgang Tichy 2/2019 */

#include "nmesh.h"
#include "evolve.h"


int nmesh_evolve(tMesh *mesh)
{
  printf("Adding evolve\n");

  /* functions */
  AddFun(POST_PARAMETERS, evolve_set_EvolveGlobals);
  AddFun(EVOLVE, evolve_output_timers);
  AddFun(EVOLVE, loadtimer_reset_mesh);
  AddFun(EVOLVE, EVOLVE_timer_start);
  AddFun(EVOLVE, evolve_myln);
  AddFun(FILTER, EVOLVE_timer_stop);
  AddFun(FILTER, evolve_filter_evosys_mesh);
  AddFun(FINALIZEMESH, evolve_free_evosys);
  AddFun(POST_PARAMETERS, evolve_free_evosys);

  /* variables */

  /* parameters */
  AddPar("evolve_method", "RK4", "[Euler,RK4,sspRK3]");
  AddPar("evolve_u_surfaces", "yes", "MPI-surface exchange for u [yes,no]"
         "(see comment before pVLList_set_surfacezones in evolve_init_evosys)");
  AddPar("evolve_loadtime", "3", "3: time RHS & limiter, 7: also time TROUBLE "
         "[3,7]"); //could extend: 1:time RHS, 2:time limiter, 4:time trouble
  AddPar("evolve_redo_troubled", "step", "what we redo in troubled elms "
         "[step,substep]");
  AddPar("evolve_trouble_BackInterpScheme", "1", "scheme we use to interp "
         "back from UNIFORM to LGL grids [1,2]. 1: Lagrange, 2: matrix R");
  AddPar("evolve_filter_all", "no", "whether we filter all evo. vars [no,yes]");
  AddPar("evolve_filter_varlist", "", "1st comps of other vars we filter");
  AddPar("evolve_filter_alp", "36", "alp in e^{-alp (i/(n0-1-dn))^s}");
  AddPar("evolve_filter_s",   "32",   "s in e^{-alp (i/(n0-1-dn))^s}");
  AddPar("evolve_filter_dn",  "0",   "dn in e^{-alp (i/(n0-1-dn))^s}");
  AddPar("evolve_output_timers", "no", "output EVOLVE-time and loadtime [no,yes]");
  AddPar("evolve_compute_change", "no", "compute change over 1 evo step [no,yes]");

  /* just a test, not needed for anything else */
  if(Getv(Par("evolve_method"), "evolve_test"))
  {
    AddFun(INITIALDATA, evolve_test_init);
    AddFun(ANALYZE, evolve_test_analyze);
    AddEvoVar("evolve_test_u", "", "test var1");
    AddEvoVar("evolve_test_v", "", "test var2");
    AddEvoVar("evolve_test_s", "", "source in v eqn");
    AddAuxVar("evolve_test_u_err", "", "error in var1");
    AddAuxVar("evolve_test_v_err", "", "error in var2");
    AddPar("evolve_method_order", "0", "expected order of convergence");
  }

  /* Old parameters that are now banned */
  if(ParLax("evolve_filter")>=0 && Getb(Par("evolve_filter")))
    BanPar("evolve_filter", "use evolve_filter_all instead");

  return 0;
}
