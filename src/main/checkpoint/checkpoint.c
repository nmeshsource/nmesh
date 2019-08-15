/* checkpoint.c */
/* Wolfgang Tichy, 8/2019 */


#include "nmesh.h"
#include "checkpoint.h"


/* save a checkpoint if the time is right for it */
int checkpoint_save_if_needed(tMesh *mesh)
{
  static double last_checkpoint_time = 0.;
  double hours      = Getd("checkpoint_hours");
  double hours_quit = Getd("checkpoint_hours_quit");
  double time       = getTimeIn_s()/3600.;
  double time_since_checkpoint;
  int do_checkpoint = 0;

  /* test if it is time */
  if(Rank0
  {
    /* test based on walltime */
    time_since_checkpoint = time - last_checkpoint_time;
    if((hours      > 0. && hours      <= time_since_checkpoint) ||
       (hours_quit > 0. && hours_quit <= time))
    {
      do_checkpoint = 1;
      last_checkpoint_time = time;
      printf("checkpoint by walltime, dt = %g hours\n", time_since_checkpoint)
    }
  }

  /* broadcast do_checkpoint from rank0 to all others */
  //

  if(do_checkpoint)
    checkpoint_save(mesh);

  return 0;
}
