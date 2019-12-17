/* nmesh_advectionFDy.h */
/* (c) Wolfgang Tichy 12/2019 */

#include "nmesh.h"
#include "advectionFDy.h"


int nmesh_advectionFDy(tMesh *mesh)
{
  if(!Getv(Par("physics"), "advectionFDy")) return 0;

  printf("Adding advectionFDy\n");

  /* functions */
  AddFun(INITIALDATA, advectionFDy_init);
  AddFun(ANALYZE, advectionFDy_analyze);

  /* variables */
  AddEvoVar("advectionFDy_u", "",     "field we advect");
  AddAuxVar("advectionFDy_u_err", "", "error in u");

  AddAuxVar("advectionFDy_x", "", "x-coord");
  AddAuxVar("advectionFDy_y", "", "y-coord");
  AddAuxVar("advectionFDy_z", "", "z-coord");

  /* parameters */
  AddPar("advectionFDy_profile", "sin", "initial profile [sin,square]");

  return 0;
}
