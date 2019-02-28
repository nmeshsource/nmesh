/* nmesh_evolve.h */
/* (c) Wolfgang Tichy 2/2019 */

#include "nmesh.h"
#include "evolve.h"


int nmesh_evolve(tMesh *mesh)
{
  printf("Adding evolve\n");

  /* functions */
  //AddFun(INITIALDATA, evolve_test);

  /* variables */
  AddEvoVar("evolve_u", "", "some var");
  AddEvoVar("evolve_v", "", "another var");
  //AddVarDim("evolve_us", "", "that needs surface exchange",7,0,6);

  /* parameters */
  //AddPar("evolve_lowlatency", "no", 
  //       "send many small rather than few large messages");

  return 0;
}
