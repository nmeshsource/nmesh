/* center.c */
/* Wolfgang Tichy 1/2024 */

#include "nmesh.h"
#include "center.h"

/* global pars for center */
tcenter center[1];


/* init global center vars */
int center_init_globals(tMesh *mesh)
{
  PRFs(":\n");
  center->center0_x = Par("center0_x");
  center->center0_y = Par("center0_y");
  center->center0_z = Par("center0_z");
  center->center1_x = Par("center1_x");
  center->center1_y = Par("center1_y");
  center->center1_z = Par("center1_z");
  center->center2_x = Par("center2_x");
  center->center2_y = Par("center2_y");
  center->center2_z = Par("center2_z");

  /* make sure some pars are saved in checkpoint */
  if(Geti(Par("center0_track_method"))
    checkpoint_save_pars_append(mesh, "center0_x center0_y center0_z");
  if(Geti(Par("center1_track_method"))
    checkpoint_save_pars_append(mesh, "center1_x center1_y center1_z");
  if(Geti(Par("center2_track_method"))
    checkpoint_save_pars_append(mesh, "center2_x center2_y center2_z");

  return 0;
}

/* update positions of centers */
int center_update(tMesh *mesh)
{
  PRFs(":\n");
  return 0;
}
