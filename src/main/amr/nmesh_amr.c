/* nmesh_amr.c */
/* Wolfgang Tichy, 1/2019 */

#include "nmesh.h"
#include "amr.h"



int nmesh_amr(tMesh *mesh)
{
  printf("Adding amr\n");

  /* functions */
/* for testing */
AddFun(POST_PARAMETERS, setup_test_mesh);

  /* variables */
  AddAuxVar("X", "",
  "coordinate 1 used for output and in spectral expansion e.g. rho");
  AddAuxVar("Y", "",
  "coordinate 2 used for output and in spectral expansion e.g. theta");
  AddAuxVar("Z", "",
  "coordinate 3 used for output and in spectral expansion e.g. z");

  /* parameters */
  AddPar("storage_verbose", "no", 
	 "verbose mode for memory allocation [no,yes]");
  return 0;
}
