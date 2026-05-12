/* nmesh_dg.c */
/* Wolfgang Tichy, 4/2019 */

#include "nmesh.h"
#include "dg.h"


int nmesh_dg(tMesh *mesh)
{
  printf("Adding dg\n");

  /* functions */
  //AddFun(POST_INITLIBS, pr_weight_ratios);
  AddFun(POST_PARAMETERS, WENOweights_init_globals);
  AddFun(FINALIZE,        WENOweights_free_globals);
  AddFun(POST_PARAMETERS, dg_set_DGglobals);
  AddFun(POST_PARAMETERS, dg_print_DGglobals);
  AddFun(PRE_EVOLVE,      dg_set_DGglobals); // in case we update fv_rec
  AddFun(FINALIZE,        dg_free_DGglobals);

  /* variables */
  //AddAuxVar("dg_u",      "",    "test function");
   
  /* parameters */
  /* Discontinous Galerkin (dg) pars */
  AddPar("dg_outerBC_flux_fac", "1 1 1", "factors in all 3 directions that "
         "can be used to switch off all surface flux terms (e.g. use "
         "'dg_outerBC_flux_fac = 1 0 0' for 1d meshes)");
  /* Finite Volume (fv) pars */
  AddPar("fv_rec", "WENOm3_2", "how we reconstruct with fv [1,WENOm3_2,"
         "WENOm5_2,WENOmT_2,WENOmZ_2,WENOm3_1,WENOm5_1,WENOmZ_1,"
         "WENOm3,WENOm5,WENOmZ,"
         "WENO3if1away_1,WENO3if2away_1,WENO3_2,WENO3_2g]");
  AddPar("fv_WENOm3_opt_weightratio", "2", "2=standard, 3=mod value");
  AddPar("fv_rec_WENOm_s1", "0", "use 0th order extrapolation if "
         "|u1-u0|/|u2-u1| < s1, [0, 0.333, ...]");
  AddPar("fv_rec_WENOm_s2", "DBL_MAX", "use 0th order extrapolation if "
         "|u1-u0|/|u2-u1| > s2, [1e99, 3, ...]");
  AddPar("fv_rec_WENOm_opt", "0", "possible options to fall back to 0th "
         "order extrapolation [0,-1,1,2] (see rec1d_u_in1_weightfac)");
  AddPar("fv_flux", "fnum", "[fnum_minus_fi, fnum]");
  AddPar("fv_divf_extrap", "dnfn_extrap1", "how we extrap div(f) [no,"
         "divf_extrap1,dnfn_extrap1,dnfn_extrap3,dnfn_extrapWENO3]");
  AddPar("fv_divf_extrap_s1", "0", "use 0th order extrapolation if "
         "|u1-u0|/|u2-u1| < 0.75 s1, [0, 0.5, ...]");
  AddPar("fv_divf_extrap_s2", "DBL_MAX", "use 0th order extrapolation if "
         "|u1-u0|/|u2-u1| > 0.75 s2, [1e50, 1e5, ...]");
  AddPar("fv_divf_extrap_opt", "0", "possible options to fall back to 0th "
         "order extrapolation [0,-1,1,2] (see rec1d_compute_1s1_u)");
  AddPar("fv_divf_use_only_right_flux", "yes", "if yes we only compute "
         "fnumR_{i} for gridpoint i and assume fnumL_{i} = -fnumR_{i-1}. "
         "Otherwise we also compute fnumL_{i} [yes,no]");
  AddPar("fv_surface_interp", "linear", "how we interpolate on node surfaces "
         "[linear,parabolic,WENO]");
  /* Pars for interpolation between fv and dg */
  AddPar("fv2dg_interp_use_extrap1", "no", "whether we extrap fv vars to "
         "faces before we interpolate from fv to dg [no,yes]");
  AddPar("fv2dg_interp_use_extrap1_vars", "", "list of vars to which the par "
         "fv2dg_interp_use_extrap1 applies");
  /* Finite Differences (fd) pars */
  /* fv_surface_interp is also used for fd */
  AddPar("fd_stencilsize", "3", "Fin.Diff. accuracy = fd_stencilsize-1");
  AddPar("fd_lopsidesize", "1", "magnitude of shift in lop-sided stencils");

  /* do some tests */
  AddPar("fv_tests", "no", "[no,yes]");
  if(Getb(Par("fv_tests")))
  {
    AddFun(EVOLVE, fv_tests); //put in EVOLVE since fv_divf needs nb. surf
    AddEvoVar("fv_test_divf", "i", "d_k f^k(u_i), where we simply set "
              "f^k(u_x) = (1,0,0), f^k(u_y) = (0,1,0), f^k(u_z) = (0,0,1) so "
              "d_k f^k(u_i) = (1/J) d_{kb} (J sqrt(g^{kb,kb} n^{kb}_i) "
              "n^{kb}_i = i-component of normal in dir Xb^{k}");
  }

  return 0;
}
