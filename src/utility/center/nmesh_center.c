/* nmesh_center.c */
/* Wolfgang Tichy 1/2024 */

#include "nmesh.h"
#include "center.h"



int nmesh_center(tMesh *mesh)
{
  printf("Adding center\n");

  /* functions */
  AddFun(POST_PARAMETERS, center_init_globals);
  AddFun(POST_EVOLVE, center_update);

  /* variables */
  //AddAuxVar("center_temp1", "", "temporary variable(for vol. integrals)");

  /* parameters */
  AddPar("center_verbose", "no", "print center pars [no,yes]");
  AddPar("center0_x", "", "x-coord of center0");
  AddPar("center0_y", "", "y-coord of center0");
  AddPar("center0_z", "", "z-coord of center0");
  AddPar("center1_x", "", "x-coord of center1");
  AddPar("center1_y", "", "y-coord of center1");
  AddPar("center1_z", "", "z-coord of center1");
  AddPar("center2_x", "", "x-coord of center2");
  AddPar("center2_y", "", "y-coord of center2");
  AddPar("center2_z", "", "z-coord of center2");
  AddPar("center0_track", "no", "[no,CM]");
  AddPar("center1_track", "no", "[no,min,max]");
  AddPar("center2_track", "no", "[no,min,max]");
  AddPar("center1_track_var", "", "variable we track [ADM_alpha,GRHD_rho0]");
  AddPar("center2_track_var", "", "variable we track [ADM_alpha,GRHD_rho0]");
  AddPar("center_track_minmove", "0.01" , "skip moves of less than value*h");
  AddPar("center1_mass", "1", "mass1 used for CM calculation");
  AddPar("center2_mass", "1", "mass2 used for CM calculation");

  return 0;
}
