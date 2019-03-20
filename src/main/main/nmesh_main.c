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

  /* variables */

  /* parameters */
  AddPar("physics", "", "what problem to solve");
  AddPar("dt"  , "0.125", "time step dt");
  
  AddPar("iterations", "0", "number of mesh iterations");
  AddPar("finaltime", "0", "iterate until mesh reaches this time");
  AddPar("iterate_parameters", "no", "whether to iterate certain parameters");

  AddPar("errorexit", "exit", "how we exit in case of error [exit,abort]");
  AddPar("logfile_creation", "append","how to create logfile [no,yes,append]");

  AddPar("timer_on", "no", "whether we use timers [yes,no]");
  AddPar("timer_MPI_barrier", "no", "barrier in timer [yes,no]");

  return 0;
}
