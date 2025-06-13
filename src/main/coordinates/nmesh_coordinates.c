/* nmesh_coordinates.h */
/* (c) Wolfgang Tichy 2/2019 */

#include "nmesh.h"
#include "coordinates.h"


int nmesh_coordinates(tMesh *mesh)
{
  printf("Adding coordinates\n");

  /* functions */
  AddFun(POST_PARAMETERS, coordinates_set_globals);
  AddFun(COORDINATES, coordinates_init); //is also called in make_child_elm!

  /* variables */
  AddAuxVar("X", "",    "coord0 in each patch e.g. lam"); // don't put any-
  AddAuxVar("Y", "",    "coord1 in each patch e.g. A");   // thing in between
  AddAuxVar("Z", "",    "coord2 in each patch e.g. B");   // these 3 lines
  AddAuxVar("x", "",    "Cartesian x coordinate"); // don't put any-
  AddAuxVar("y", "",    "Cartesian y coordinate"); // thing in between
  AddAuxVar("z", "",    "Cartesian z coordinate"); // these 3 lines
  AddAuxVar("dXd", "i", "coord derivative dX/dx^i"); // don't put any-
  AddAuxVar("dYd", "i", "coord derivative dY/dx^i"); // thing in between
  AddAuxVar("dZd", "i", "coord derivative dZ/dx^i"); // these 3 lines
  AddAuxVar("det_dXbdx", "", "determinant of dXb/dx");
  AddAuxVar("sqrtdet2g_o_det3gamma","@", "sqrt(det(2g)/det(3gamma_ij)), here "
            "2g is 2-metric induced on surface_f of node in Xb-coords and "
            "3gamma_ij is 3-metric in x-coords");
  AddAuxVar("sqrtgdiag", "I", "sqrt of diagonal components of upper index "
            "3-metric in Xb-coords");
  AddAuxVar("coordinates_tmp1", "", "temp. space"); // don't put any-
  AddAuxVar("coordinates_tmp2", "", "temp. space"); // thing in between
  AddAuxVar("coordinates_tmp3", "", "temp. space"); // the tmp var lines
  //AddAuxVar("oC0_", "@", "coord0 on faces"); // don't put anything
  //AddAuxVar("oC1_", "@", "coord1 on faces"); // between these 2 lines

  /* Vars that contain cub. sph. sigma_{0/1} and their derivs at grid points.
     We should not set these six explicitly, rather we use the funcs
     CI->FSurf[0],CI->dFSurfdC[0], CI->FSurf[1],CI->dFSurfdC[1] to set them. */
  AddAuxVarDim("CubedSphere_sigma0",     "", "sigma_{0}", 1,-1,-1);
  AddAuxVarDim("CubedSphere_dsigma0_dA", "", "d/dA sigma_{0}", 1,-1,-1);
  AddAuxVarDim("CubedSphere_dsigma0_dB", "", "d/dB sigma_{0}", 1,-1,-1);
  AddAuxVarDim("CubedSphere_sigma1",     "", "sigma_{1}", 1,-1,-1);
  AddAuxVarDim("CubedSphere_dsigma1_dA", "", "d/dA sigma_{1}", 1,-1,-1);
  AddAuxVarDim("CubedSphere_dsigma1_dB", "", "d/dB sigma_{1}", 1,-1,-1);

  /* parameters */
  AddPar("coordinates_verbose", "yes", "verbose [yes,no]");
  AddPar("coordinates_3metric", "", "[flat,any var name]");
  AddPar("coordinates_surface_metric", "sqrtgdiag", "metric for faces, set "
         "e.g. in coordinates_init_node [sqrtgdiag,sqrtdet2g_o_det3gamma]");

  /* special CubedSphere pars and vars */
  AddPar("CubedSphere_sigma01_lmax", "8", "lmax for Ylm's used in the "
         "func FSurf_CubSph_sigma01");
  AddPar("CubedSphere_sigma01_def", "no", "switch on sigma01_def [no,yes]");
  AddPar("CubedSphere_sigma01_test", "no", "run tests [no,yes]");
  if(Getb(Par("CubedSphere_sigma01_test")))
  {
    Sets(Par("CubedSphere_sigma01_def"), "yes");
    AddFun(POST_COORDINATES, CubSphTest_CI_Fcoef0_for_deformed_sigma);
  }
  if(Getb(Par("CubedSphere_sigma01_def")))
  {
    /* Var to we use to compute CI->FSurf[0]: */
    AddAuxVar("CubedSphere_sigma0_def",    "", "var we use to define and set "
              "the sigma_0 in CubedSphere_sigma0*");
    /* The strategy is to put the Ylm coeffs of CubedSphere_sigma0_def into
       CI->Fcoef[0] for the 1st patch p0. CI->Fcoef[1] of the enclosed patch
       p0-6 should be the same as CI->Fcoef[0] from patch p0. We thus use
       CI->Fcoef[0] from patch p0 instead of CI->Fcoef[1] for patch p0-6.
       Fcoef[0] is then used in the funcs FSurf and dFSurfdC.
       Note that CubedSphere_sigma0_def needs to have storage only in the
       elements that touch face0 of patch p. */
    /* As in sgrid we use more than 1 point in the radial direction to have
       room to store sigma0 from previous iterations. To checkpoint
       CubedSphere_sigma0_def we make it a DATAVAR: */
    MeshVarSetType(mesh, Ind("CubedSphere_sigma0_def"), DATAVAR);

    /* For some unknown reason there is a large error in the coefs of
       Ylm expansiom if we go beyond lmax=3. */
    if(Geti(Par("CubedSphere_sigma01_lmax"))>=4)
      errorexit("CubedSphere_sigma01_lmax is too large! For some unknown "
                "reason the Fcoef for l=4,m=0 has a large error...");
  }

  /* Variables to store some things in the middle between two regular grid
     points of a node */
  AddPar("coordinates_midpoint_data", "no", "[no,yes,yes all]");
  if(Getb(Par("coordinates_midpoint_data")))
  {
    AddAuxVar("Xm_dXd", "i", "dX/dx^i on X mid points"); // don't put any-
    AddAuxVar("Xm_dYd", "i", "dY/dx^i on X mid points"); // thing in between
    AddAuxVar("Xm_dZd", "i", "dZ/dx^i on X mid points"); // these 3 lines
    AddAuxVar("Ym_dXd", "i", "dX/dx^i on Y mid points"); // don't put any-
    AddAuxVar("Ym_dYd", "i", "dY/dx^i on Y mid points"); // thing in between
    AddAuxVar("Ym_dZd", "i", "dZ/dx^i on Y mid points"); // these 3 lines
    AddAuxVar("Zm_dXd", "i", "dX/dx^i on Z mid points"); // don't put any-
    AddAuxVar("Zm_dYd", "i", "dY/dx^i on Z mid points"); // thing in between
    AddAuxVar("Zm_dZd", "i", "dZ/dx^i on Z mid points"); // these 3 lines
    AddAuxVar("Xm_det_dXbdx", "", "det(dXb/dx) on X mid points");
    AddAuxVar("Ym_det_dXbdx", "", "det(dXb/dx) on Y mid points");
    AddAuxVar("Zm_det_dXbdx", "", "det(dXb/dx) on Z mid points");
    AddAuxVar("Xm_sqrtgdiag", "I", "sqrt of diag comps of upper index 3-"
              "metric in Xb-coords on mid points in X-dir");
    AddAuxVar("Ym_sqrtgdiag", "I", "sqrt of diag comps of upper index 3-"
              "metric in Xb-coords on mid points in Y-dir");
    AddAuxVar("Zm_sqrtgdiag", "I", "sqrt of diag comps of upper index 3-"
              "metric in Xb-coords on mid points in Z-dir");
    /* Maybe we only need Xm_dXdx^i, Ym_dYdx^i, Zm_dZd^i and
       Xm_sqrtgdiagx, Ym_sqrtgdiagy, Zm_sqrtgdiagz ??? Then we could
       drastically reduce the number of variables! */
  }

  /* do some tests */
  AddPar("coordinates_tests", "no", "[no,yes]");
  if(Getb(Par("coordinates_tests")))
  {
    AddFun(POST_COORDINATES, coordinates_tests);
    AddAuxVar("divb_J_sqrtgdiag_n","i", "d_{ib}(J sqrtgdiag^{ib} n^{ib}_i)");
    AddAuxVar("ooJ_Db_J_sqrtgdiag_n","i",
              "(1/J) D_{ib} (J sqrtgdiag^{ib} n^{ib}_i)");
  }

  return 0;
}
