/* nmesh_dg.c */
/* Wolfgang Tichy, 4/2019 */

#include "nmesh.h"
#include "dg.h"


int nmesh_dg(tMesh *mesh)
{
  printf("Adding dg\n");

  /* functions */
  //AddFun(FIRST, pr_weight_ratios);
  AddFun(POST_PARAMETERS, WENOweights_init_globals);
  AddFun(FINALIZE,        WENOweights_free_globals);
  AddFun(POST_PARAMETERS, dg_set_DGglobals);
  AddFun(PRE_EVOLVE,      dg_set_DGglobals); // in case we update fv_rec

  /* variables */
  //AddAuxVar("dg_u",      "",    "test function");
   
  /* parameters */
  AddPar("dg_outerBC_flux_fac", "1 1 1", "factors in all 3 directions that "
         "can be used to switch off all surface flux terms (e.g. use "
         "'dg_outerBC_flux_fac = 1 0 0' for 1d meshes)");
  AddPar("dg_numerical_flux", "LLF", "numerical flux [LLF]");
  AddPar("fv_rec", "1", "how we reconstruct with fv "
         "[1,WENOm3_2,WENO3if1away_1,WENO3if2away_1,WENO3_2,WENO3_2g]");

  return 0;
}
