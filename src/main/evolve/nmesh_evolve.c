/* nmesh_evolve.h */
/* (c) Wolfgang Tichy 2/2019 */

#include "nmesh.h"
#include "evolve.h"


int nmesh_evolve(tMesh *mesh)
{
  printf("Adding evolve\n");

  /* functions */
  AddFun(EVOLVE, evolve_myln);
  AddFun(FINALIZEMESH, evolve_finalize);

  /* variables */

  /* parameters */
  //AddPar("evolve_beans", "no", "[no,yes]");

  /* just a test, not needed for anything else */
  if(Getv(Par("physics"), "evolve_test"))
  {
    AddFun(INITIALDATA, evolve_test_init);
    AddFun(ANALYZE, evolve_test_analyze);
    AddEvoVar("evolve_u", "", "test var1");
    AddEvoVar("evolve_v", "", "test var2");
    AddEvoVar("evolve_s", "", "source in v eqn");
    AddAuxVar("evolve_u_err", "", "error in var1");
    AddAuxVar("evolve_v_err", "", "error in var2");
  }

  return 0;
}
