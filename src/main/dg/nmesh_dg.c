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
  AddFun(POST_PARAMETERS, dg_print_DGglobals);
  AddFun(PRE_EVOLVE,      dg_set_DGglobals); // in case we update fv_rec

  /* variables */
  //AddAuxVar("dg_u",      "",    "test function");
   
  /* parameters */
  /* Discontinous Galerkin (dg) pars */
  AddPar("dg_outerBC_flux_fac", "1 1 1", "factors in all 3 directions that "
         "can be used to switch off all surface flux terms (e.g. use "
         "'dg_outerBC_flux_fac = 1 0 0' for 1d meshes)");
  AddPar("dg_numerical_flux", "LLF", "numerical flux [LLF]");
  /* Finite Volume (fv) pars */
  AddPar("fv_rec", "WENOm3_2", "how we reconstruct with fv [1,WENOm3_2,"
         "WENOm5_2,WENOmZ_2,WENOm3_1,WENOm5_1,WENOmZ_1,"
         "WENO3if1away_1,WENO3if2away_1,WENO3_2,WENO3_2g]");
  AddPar("fv_WENOm3_opt_weightratio", "2", "2=standard, 3=mod value");
  AddPar("fv_flux", "fnum_minus_fi", "[fnum_minus_fi, fnum]");
  AddPar("fv_divf_extrap", "dnfn_extrap1", "how we extrap div(f) [no,"
         "divf_extrap1,dnfn_extrap1]");
  AddPar("fv_divf_extrap_s1", "0", "use 0th order extrapolation if "
         "|u1-u0|/|u2-u1| < 0.75 s1, [0, 0.5, ...]");
  AddPar("fv_divf_extrap_s2", "DBL_MAX", "use 0th order extrapolation if "
         "|u1-u0|/|u2-u1| > 0.75 s2, [1e50, 1e5, ...]");
  AddPar("fv_divf_use_only_right_flux", "yes", "if yes we only compute fnumR "
         "and assume fnumL=-fnumR, otherwise we also compute fnumL [yes,no]");
  AddPar("fv_surface_interp", "linear", "how we interpolate on node surfaces "
         "[linear,parabolic]");
  /* Finite Differences (fd) pars */
  /* fv_surface_interp is also used for fd */
  AddPar("fd_stencilsize", "3", "Fin.Diff. accuracy = fd_stencilsize-1");
  AddPar("fd_lopsidesize", "1", "magnitude of shift in lop-sided stencils");

  return 0;
}
