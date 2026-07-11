/* nmesh_units.c */
/* Wolfgang Tichy 6/2019 */

#include "nmesh.h"
#include "units.h"



int nmesh_units(tMesh *mesh)
{
  printf("Adding units\n");

  /* functions */
  AddFun(POST_INITLIBS, units_set_SItoGMc1_mesh);
  AddFun(POST_INITIALDATA, print_unit_conversion_factors);

  /* variables */
  //AddAuxVar("units_temp1", "", "temporary variable(for vol. integrals)");

  /* parameters */
  //AddPar("0dunits", "", "variables set units for");

  return 0;
}
