/* nmesh_main.c */
/* Wolfgang Tichy, 1/2019 */

#include "nmesh.h"
#include "main.h"



int nmesh_main(tMesh *mesh)
{
  printf("Adding main\n");

  /* functions */
  AddFun(OUTPUT, write_all_timers);
  AddFun(FINALIZEMESH, free_all_timers);
  AddFun(PRE_COORDINATES, nmesh_update_parameters);
  AddFun(PRE_EVOLVE, nmesh_update_parameters);

  /* variables */

  /* parameters */
  AddPar("physics", "", "what problem to solve");
  AddPar("dt", "auto", "time step dt [#,auto], for auto: dt = dtfac * hmin");
  AddPar("dtfac", "0.25", "Courant factor in: dt = dtfac * hmin");

  AddPar("iterations", "0", "number of mesh iterations");
  AddPar("finaltime", "0", "iterate until mesh reaches this time");
  AddPar("iterate_parameters", "no", "whether to iterate certain parameters");

  AddPar("errorexit", "exit", "how we exit in case of error [exit,abort]");
  AddPar("logfile_creation", "yes","create logfile [yes,no]");

  AddPar("timer_on", "no", "whether we use timers [yes,no]");
  AddPar("timer_MPI_barrier", "no", "barrier in timer [yes,no]");

  return 0;
}
