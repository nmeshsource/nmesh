/* nmesh_checkpoint.c */
/* Wolfgang Tichy, 8/2019 */

#include "nmesh.h"
#include "checkpoint.h"


int nmesh_checkpoint(tMesh *mesh)
{
  printf("Adding checkpoint\n");

  /* functions */

  /* variables */

  /* parameters */
  AddPar("checkpoint_hours", "-1", "hours after which we checkpoint");
  AddPar("checkpoint_hours_quit", "-1", "when to quit after checkpoint");
  AddPar("checkpoint_save_pars", "", "pars in this list will be saved, all "
         "others come from parfile or earlier call to initialize_libraries "
         "when we restart");

  return 0;
}
