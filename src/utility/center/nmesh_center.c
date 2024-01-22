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
  AddPar("center0_x", "", "x-coord of center0"); //
  AddPar("center0_y", "", "y-coord of center0"); // KEEP these
  AddPar("center0_z", "", "z-coord of center0"); //
  AddPar("center1_x", "", "x-coord of center1"); // center pars in
  AddPar("center1_y", "", "y-coord of center1"); //
  AddPar("center1_z", "", "z-coord of center1"); // this particular
  AddPar("center2_x", "", "x-coord of center2"); //
  AddPar("center2_y", "", "y-coord of center2"); // order !!!
  AddPar("center2_z", "", "z-coord of center2"); //
  AddPar("center0_track", "no", "[no,CM]");
  AddPar("center1_track", "no", "[no,min,max]");
  AddPar("center2_track", "no", "[no,min,max]");
  AddPar("center1_track_var", "", "variable we track [ADM_alpha,GRHD_rho0]");
  AddPar("center2_track_var", "", "variable we track [ADM_alpha,GRHD_rho0]");
  AddPar("center_track_minmove", "0.01" , "don't move more than value*h");
  AddPar("center1_mass", "1", "mass1 used for CM calculation");
  AddPar("center2_mass", "1", "mass2 used for CM calculation");

  return 0;
}
