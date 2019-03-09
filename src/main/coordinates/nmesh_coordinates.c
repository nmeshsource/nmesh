/* nmesh_coordinates.h */
/* (c) Wolfgang Tichy 2/2019 */

#include "nmesh.h"
#include "coordinates.h"


int nmesh_coordinates(tMesh *mesh)
{
  printf("Adding coordinates\n");

  /* functions */
  AddFun(COORDINATES, coordinates_init);

  /* variables */
  AddAuxVar("X", "",    "coord0 in each patch e.g. lam"); // don't put any-
  AddAuxVar("Y", "",    "coord1 in each patch e.g. A");   // thing in between
  AddAuxVar("Z", "",    "coord2 in each patch e.g. B");   // these 3 lines
  AddAuxVar("dXd", "i", "coord derivative dX/dx^i"); // don't put any-
  AddAuxVar("dYd", "i", "coord derivative dY/dx^i"); // thing in between
  AddAuxVar("dZd", "i", "coord derivative dZ/dx^i"); // these 3 lines
  AddAuxVar("det_dXbdx", "", "determinant of dXb/dx");
  AddAuxVar("x", "", "Cartesian x coordinate"); // don't put any-
  AddAuxVar("y", "", "Cartesian y coordinate"); // thing in between
  AddAuxVar("z", "", "Cartesian z coordinate"); // these 3 lines

  /* create vars that contain cub. sph. sigma_{0/1} and their derivs */
  AddAuxVarDim("CubedSphere_sigma01",     "", "sigma_{0/1}", 2,-1,-1);
  AddAuxVarDim("CubedSphere_dsigma01_dA", "", "d/dA sigma_{0/1}", 2,-1,-1);
  AddAuxVarDim("CubedSphere_dsigma01_dB", "", "d/dB sigma_{0/1}", 2,-1,-1);

  /* parameters */
  AddPar("Coordinates_verbose", "yes", "verbose [yes,no]");
  AddPar("CubedSphere_sigma01_lmax", "8", "lmax for Ylm's "
         "used in FSurf_CubSph_sigma01_func");

  return 0;
}
