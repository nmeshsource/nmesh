/* nmesh_evolve.h */
/* (c) Wolfgang Tichy 2/2019 */

#include "nmesh.h"
#include "evolve.h"


int nmesh_evolve(tMesh *mesh)
{
  printf("Adding evolve\n");

  /* functions */
  AddFun(EVOLVE, evolve_myln);
  AddFun(FILTER, evolve_filter_evosys_mesh);
  AddFun(FINALIZEMESH, evolve_free_evosys);
  AddFun(POST_PARAMETERS, evolve_free_evosys);

  /* variables */

  /* parameters */
  AddPar("evolve_method", "RK4", "[Euler,RK4,sspRK3]");
  AddPar("evolve_filter", "no", "whether we filter all evo. vars [no,yes]");
  AddPar("evolve_filter_varlist", "", "1st comps of other vars we filter");
  AddPar("evolve_filter_alp", "36", "alp in e^{-alp (i/(n0-1))^s}");
  AddPar("evolve_filter_s",   "32",   "s in e^{-alp (i/(n0-1))^s}");
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

  return 0;
}
