/* nmesh_center.c */
/* Wolfgang Tichy 1/2024 */

#include "nmesh.h"
#include "center.h"



int nmesh_center(tMesh *mesh)
{
  printf("Adding center\n");

  /* functions */
  AddFun(POST_PARAMETERS, center_init_globals);
  AddFun(ANALYZE, center_update);

  /* variables */
  //AddAuxVar("center_temp1", "", "temporary variable(for vol. integrals)");

  /* parameters */
  AddPar("center0_x", "0", "x-coord of center0");
  AddPar("center0_y", "0", "y-coord of center0");
  AddPar("center0_z", "0", "z-coord of center0");
  AddPar("center1_x", "0", "x-coord of center1");
  AddPar("center1_y", "0", "y-coord of center1");
  AddPar("center1_z", "0", "z-coord of center1");
  AddPar("center2_x", "0", "x-coord of center2");
  AddPar("center2_y", "0", "y-coord of center2");
  AddPar("center2_z", "0", "z-coord of center2");
  AddPar("center0_track_method", "0", "0: do not track, 3: track CM");
  AddPar("center1_track_method", "0", "0: do not track, 1/2: track max/min");
  AddPar("center2_track_method", "0", "0: do not track, 1/2: track max/min");
  AddPar("center1_track_var", "", "variable we track [ADM_alpha,GRHD_rho0]");
  AddPar("center2_track_var", "", "variable we track [ADM_alpha,GRHD_rho0]");

  return 0;
}
