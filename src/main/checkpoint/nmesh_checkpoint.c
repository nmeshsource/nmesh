/* nmesh_checkpoint.c */
/* Wolfgang Tichy, 8/2019 */

#include "nmesh.h"
#include "checkpoint.h"


int nmesh_checkpoint(tMesh *mesh)
{
  printf("Adding checkpoint\n");

  /* functions */
  /* NOTE: checkpoint is not scheduled here because it is directly called
           by function evolve_mesh from main.c. In evolve_mesh we have:
     ...
     RunFun(ANALYZE);
     ...
     checkpoint_save_if_needed(mesh, 0);
     ...  */

  /* variables */

  /* parameters */
  AddPar("checkpoint", "no", "whether we save checkpoints [no,yes,load_mesh]");
  AddPar("checkpoint_hours", "-1", "hours after which we checkpoint");
  AddPar("checkpoint_hours_quit", "-1", "when to quit after checkpoint");
  AddPar("checkpoint_save_pars", "", "pars in this list will be saved, all "
         "others come from parfile or earlier call to initialize_libraries "
         "when we restart");

  return 0;
}
