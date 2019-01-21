/* nmesh_amr.c */
/* Wolfgang Tichy, 1/2019 */

#include "nmesh.h"
#include "amr.h"



int nmesh_amr(tMesh *mesh)
{
  printf("Adding amr\n");

  /* variables */
  AddVar("X", "",
  "coordinate 1 used for output and in spectral expansion e.g. rho");
  AddVar("Y", "",
  "coordinate 2 used for output and in spectral expansion e.g. theta");
  AddVar("Z", "",
  "coordinate 3 used for output and in spectral expansion e.g. z");

  /* parameters */
  AddPar("storage_verbose", "no", 
	 "verbose mode for memory allocation [no,yes]");
  return 0;
}
