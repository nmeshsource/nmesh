/* nmesh_misc.h */
/* (c) Wolfgang Tichy 2/2019 */

#include "nmesh.h"
#include "misc.h"


int nmesh_misc(tMesh *mesh)
{
  if(!Getv(Par("physics"), "misc")) return 0;

  printf("Adding misc\n");

  /* functions */
  AddFun(INITIALDATA, misc_test);

  /* variables */
  AddEvoVar("misc_u", "", "some var");
  AddEvoVar("misc_v", "", "another var");
  //AddVarDim("misc_us", "", "that needs surface exchange",7,0,6);

  /* parameters */
  AddPar("misc_ajsurf_v_init", "test_func(x,y,z)", "how we init v test_ajsurf"
         "[test_func(x,y,z),test_func(lam,A,B)]");

  return 0;
}
