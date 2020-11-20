/* nmesh_dg.c */
/* Wolfgang Tichy, 4/2019 */

#include "nmesh.h"
#include "dg.h"


int nmesh_dg(tMesh *mesh)
{
  printf("Adding dg\n");

  /* functions */
//AddFun(FIRST, pr_weight_ratios);
  //AddFun(INITIALDATA, dg_startup);
  //AddFun(ANALYZE, dg_analyze);

  /* variables */
  //AddAuxVar("dg_u",      "",    "test function");
   
  /* parameters */
  AddPar("dg_outerBC_flux_fac", "1 1 1", "factors in all 3 directions that "
         "can be used to switch off all surface flux terms (e.g. use "
         "'dg_outerBC_flux_fac = 1 0 0' for 1d meshes)");
  AddPar("dg_numerical_flux", "LLF", "numerical flux [LLF]");

  return 0;
}
