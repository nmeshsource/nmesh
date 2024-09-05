/* nmesh_center.c */
/* Wolfgang Tichy 1/2024 */

#include "nmesh.h"
#include "center.h"



int nmesh_center(tMesh *mesh)
{
  printf("Adding center\n");

  /* functions */
  AddFun(INITMESH, center_amr);
  AddFun(POST_PARAMETERS, center_init_globals);
  AddFun(POST_EVOLVE, center_update);
  AddFun(POST_EVOLVE, center_1_2_distance_output);
  AddFun(AMR, center_amr);

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
  /* for AMR following centers 1 & 2 */
  AddPar("center_amr_time", "-1", "when to call center_amr (-1=never)");
  AddPar("center1_amr_lmax", "0", "finest level at center1");
  AddPar("center1_amr_radius", "8", "radius of finest level around center1");
  AddPar("center2_amr_lmax", "center1_amr_lmax", "finest level at center2");
  AddPar("center2_amr_radius", "8", "radius of finest level around center2");
  /* for computing distances */
  AddPar("center_distance_output_time", "-1", "time for distance (-1=never)");
  AddPar("center_distance_radius1", "0", "radius we exclude from distance");
  AddPar("center_distance_radius2", "0", "radius we exclude from distance");
  AddPar("center_distance_metric", "ADM_gxx", "metric for proper distance");

  return 0;
}
