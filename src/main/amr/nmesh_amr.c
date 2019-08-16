/* nmesh_amr.c */
/* Wolfgang Tichy, 1/2019 */

#include "nmesh.h"
#include "amr.h"



int nmesh_amr(tMesh *mesh)
{
  printf("Adding amr\n");

  /* functions */
  //AddFun(POST_PARAMETERS, setup_mesh);

  /* variables */
  //AddAuxVar("Xb", "", "coord0 inside each node"); // we may not need these???
  //AddAuxVar("Yb", "", "coord1 inside each node");
  //AddAuxVar("Zb", "", "coord2 inside each node");

  /* parameters */
  AddPar("amr_mesh_type", "test_mesh",
	 "mesh we start with [test_mesh,l2_mesh]");
  AddPar("amr_n0", "5", "number of points in dir. 0 in one node");
  AddPar("amr_n1", "5", "number of points in dir. 1 in one node");
  AddPar("amr_n2", "5", "number of points in dir. 2 in one node");
  AddPar("amr_nmax", "55", "max number of points in all 3 dir. in one node");
  AddPar("amr_luni", "0",  "level up to which each patch is refined initially");
  AddPar("amr_refine_p", "-1", "patch that we refine one level further");
  AddPar("amr_BoxMesh_xc", "0",   "center for boxes");
  AddPar("amr_BoxMesh_dout", "1", "box radius");

  AddPar("bface_options", "face2_order3", "how we set some bface flags "
         "[none,face2_order0,face2_order1,face2_order2,face2_order3]");

  return 0;
}
