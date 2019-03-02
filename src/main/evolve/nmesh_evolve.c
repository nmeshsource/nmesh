/* nmesh_evolve.h */
/* (c) Wolfgang Tichy 2/2019 */

#include "nmesh.h"
#include "evolve.h"


int nmesh_evolve(tMesh *mesh)
{
  printf("Adding evolve\n");

  /* functions */
  AddFun(EVOLVE, evolve_myln);
  AddFun(FINALIZEMESH, evolve_free_evosys);
  AddFun(POST_PARAMETERS, evolve_free_evosys);

  /* variables */

  /* parameters */
  AddPar("evolve_method", "RK", "[Euler,RK]");

  /* just a test, not needed for anything else */
  if(Getv(Par("physics"), "evolve_test"))
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
