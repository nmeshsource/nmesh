/* nmesh_coordinates.h */
/* (c) Wolfgang Tichy 2/2019 */

#include "nmesh.h"
#include "coordinates.h"


int nmesh_coordinates(tMesh *mesh)
{
  printf("Adding coordinates\n");

  /* functions */
  AddFun(POST_PARAMETERS, coordinates_set_globals);
  AddFun(COORDINATES, coordinates_init); //is also called in make_child_node!

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

  /* create vars that contain cub. sph. sigma_{0/1} and their derivs */
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
  AddPar("CubedSphere_sigma01_lmax", "8", "lmax for Ylm's "
         "used in FSurf_CubSph_sigma01_func");

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

  return 0;
}
