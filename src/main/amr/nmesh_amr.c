/* nmesh_amr.c */
/* Wolfgang Tichy, 1/2019 */

#include "nmesh.h"
#include "amr.h"



int nmesh_amr(tMesh *mesh)
{
  printf("Adding amr\n");

  /* functions */
  AddFun(POST_PARAMETERS, setup_mesh);

  /* variables */
  //AddAuxVar("Xb", "", "coord0 inside each node"); // we may not need these???
  //AddAuxVar("Yb", "", "coord1 inside each node");
  //AddAuxVar("Zb", "", "coord2 inside each node");

  /* parameters */
  AddPar("storage_verbose", "no",
         "verbose mode for memory allocation [no,yes]");
  AddPar("amr_mesh_type", "test_mesh",
	 "mesh we start with [test_mesh,l2_mesh]");
  return 0;
}
