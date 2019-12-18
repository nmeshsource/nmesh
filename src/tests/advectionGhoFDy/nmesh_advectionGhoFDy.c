/* nmesh_advectionGhoFDy.h */
/* (c) Wolfgang Tichy 12/2019 */

#include "nmesh.h"
#include "advectionGhoFDy.h"


int nmesh_advectionGhoFDy(tMesh *mesh)
{
  if(!Getv(Par("physics"), "advectionGhoFDy")) return 0;

  printf("Adding advectionGhoFDy\n");

  /* functions */
  AddFun(INITIALDATA, advectionGhoFDy_init);
  AddFun(ANALYZE, advectionGhoFDy_analyze);

  /* variables */
  AddEvoVar("advectionGhoFDy_u", "",     "field we advect");
  AddAuxVar("advectionGhoFDy_u_err", "", "error in u");

  AddAuxVar("advectionGhoFDy_x", "", "x-coord");
  AddAuxVar("advectionGhoFDy_y", "", "y-coord");
  AddAuxVar("advectionGhoFDy_z", "", "z-coord");

  /* parameters */
  AddPar("advectionGhoFDy_profile", "sin", "initial profile [sin,square]");

  return 0;
}
