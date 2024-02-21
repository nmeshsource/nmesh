/* nmesh_basis.h */
/* (c) Wolfgang Tichy 2/2019 */

#include "nmesh.h"
#include "basis.h"


int nmesh_basis(tMesh *mesh)
{
  printf("Adding basis\n");

  /* functions */
  AddFun(POST_INITLIBS, init_gridpoints);
  AddFun(POST_FINALIZEMESH, free_gridpoints);
  AddFun(POST_PARAMETERS, basis_init_globals);

  /* variables */
  //AddVar("basis_temp1", "", "temporary variable");

  /* parameters */
  //AddPar("basis0", "Legendre", "basis functions we use");
  AddPar("basis_expfilter_JacobianPower", "1",
         "multiply var with Jacobian to this power before filtering it");
  AddPar("basis_filter_fv", "yes", "filter in fv nodes [yes,no]");
  AddPar("basis_Lagrange_interp_np", "n", "interp point number [n,6,8,10]");
  AddPar("basis_WENO_interp_np", "6", "interp point number [4,6]");

  return 0;
}
