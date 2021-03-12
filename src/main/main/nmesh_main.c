/* nmesh_main.c */
/* Wolfgang Tichy, 1/2019 */

#include "nmesh.h"
#include "main.h"



int nmesh_main(tMesh *mesh)
{
  printf("Adding main\n");

  /* functions */
  AddFun(FIRST, print_endian_info);
  AddFun(OUTPUT, write_all_timers);
  AddFun(POST_OUTPUT, nan_checker);
  AddFun(POST_OUTPUT, sysmon);
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

  AddPar("nan_check", "", "variables we check for NANs [varnames]");
  AddPar("sysmon_hours", "0.5", "hours after which we write sysmon.log");
  AddPar("sysmon_output_per_rank", "no","output data for each rank [no,yes]");
  AddPar("update_parameters_hours", "0.2",
         "hours after which we read nmesh_update_parameters.par");

  AddPar("fs_bufsize", "1048576", "use setvbuf with this size for IO");
  AddPar("fs_sync", "no", "whether fs_sync is active [yes,no]");

  return 0;
}
