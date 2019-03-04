/* nmesh_advection1.h */
/* (c) Wolfgang Tichy 2/2019 */

#include "nmesh.h"
#include "advection1.h"


int nmesh_advection1(tMesh *mesh)
{
  printf("Adding advection1\n");

  /* functions */
  AddFun(INITIALDATA, advection1_init);
  AddFun(ANALYZE, advection1_analyze);

  /* variables */
  AddEvoVar("advection1_u", "",  "field we advect");
  AddAuxVar("advection1_u", "i", "deriv of u");
  AddAuxVar("advection1_u_err", "", "error in u");

  /* parameters */
  AddPar("advection1_direction", "1 0 0", "propagation direction");

  return 0;
}
