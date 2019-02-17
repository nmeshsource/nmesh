/* nmesh_coordinates.h */
/* (c) Wolfgang Tichy 2/2019 */

#include "nmesh.h"
#include "coordinates.h"


int nmesh_coordinates(tMesh *mesh)
{
  printf("Adding coordinates\n");

  /* functions */
  AddFun(COORDINATES, coordinates_init);
  AddFun(REINIT, coordinates_reinit);

  /* variables */
  AddAuxVar("X", "",    "coordinate 0 in each patch e.g. lambda");
  AddAuxVar("Y", "",    "coordinate 1 in each patch e.g. A");
  AddAuxVar("Z", "",    "coordinate 2 in each patch e.g. B");
  AddAuxVar("dXd", "i", "coordinate derivative dX/dx^i"); // don't put any-
  AddAuxVar("dYd", "i", "coordinate derivative dY/dx^i"); // thing in between
  AddAuxVar("dZd", "i", "coordinate derivative dZ/dx^i"); // these 3 lines
  AddAuxVar("x", "", "Cartesian x coordinate");
  AddAuxVar("y", "", "Cartesian y coordinate");
  AddAuxVar("z", "", "Cartesian z coordinate");
 
  /* parameters */
  AddPar("Coordinates_verbose", "yes", "verbose [yes,no]");

  return 0;
}
