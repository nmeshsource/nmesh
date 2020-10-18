/* nmesh_basis.h */
/* (c) Wolfgang Tichy 2/2019 */

#include "nmesh.h"
#include "basis.h"


int nmesh_basis(tMesh *mesh)
{
  printf("Adding basis\n");

  /* functions */
  AddFun(FIRST, init_gridpoints);
  AddFun(FINALIZEMESH, free_gridpoints);
  AddFun(POST_PARAMETERS, basis_init_globals);

  /* variables */
  //AddVar("basis_temp1", "", "temporary variable");

  /* parameters */
  //AddPar("basis0", "Legendre", "basis functions we use");
  AddPar("basis_expfilter_JacobianPower", "1",
         "multiply var with Jacobian to this power before filtering it");

  return 0;
}
