/* nmesh_coordinates.h */
/* (c) Wolfgang Tichy 2/2019 */

#include "nmesh.h"
#include "coordinates.h"


int nmesh_coordinates(tMesh *mesh)
{
  printf("Adding coordinates\n");

  /* functions */
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
  AddAuxVar("sqrtdet2gamma","@", "sqrt(det(2gamma)), here 2gamma is 2-metric "
            "induced on surface_f of node in Xb-coords");
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
  AddPar("coordinates_surface_metric", "sqrtdet2gamma", "surface metric info "
         "we set in coordinates_init_node [sqrtdet2gamma,sqrtgdiag]");
  AddPar("CubedSphere_sigma01_lmax", "8", "lmax for Ylm's "
         "used in FSurf_CubSph_sigma01_func");

  return 0;
}
