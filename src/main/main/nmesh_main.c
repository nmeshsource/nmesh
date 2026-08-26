/* nmesh_main.c */
/* Wolfgang Tichy, 1/2019 */

#include "nmesh.h"
#include "main.h"



int nmesh_main(tMesh *mesh)
{
  printf("Adding main (rev %s)\n", RevStr(MAINREV));

  /* functions */
  AddFun(POST_INITLIBS, print_endian_info);
  AddFun(POST_INITLIBS, check_compiledphysics);
  AddFun(POST_INITLIBS, CheckForBannedPars);
  AddFun(OUTPUT, write_all_timers);
  AddFun(POST_OUTPUT, nan_checker);
  AddFun(POST_OUTPUT, sysmon);
  AddFun(FINALIZEMESH, free_all_timers);
  AddFun(POST_FINALIZEMESH, FreeBannedParList);
  AddFun(PRE_COORDINATES, nmesh_update_parameters);
  AddFun(PRE_EVOLVE, nmesh_update_parameters);

  /* variables */

  /* parameters */
  AddPar("physics", "", "what problem to solve");
  AddPar("dt", "auto", "time step dt [#,auto,auto2]. auto computes dt. "
         "auto: dt=dtfac*hmin, auto2 (from CFL limit): dt=dtfac*dtlim");
  AddPar("dtfac", "0.25", "Courant factor in: dt = dtfac * hmin");
  AddPar("uniform_dtfac", "0.125", "Courant factor on uniform grids where: "
         "dt = uniform_dtfac * hmin");
  AddPar("iterations", "0", "number of mesh iterations");
  AddPar("finaltime", "0", "iterate until mesh reaches this time");
  AddPar("iterate_parameters", "no", "whether to iterate certain parameters");

  AddPar("errorexit", "exit", "how we stop on errors [exit,abort;errno]");
  AddPar("logfile_creation", "no","create logfile [yes,no,/dev/null,none]");
  AddPar("logfile_reproducible", "no","printf only stuff that is the same "
         "for each run (e.g. no timing) [yes,no]");

  AddPar("timer_on", "no", "whether we use timers [yes,no]");
  AddPar("timer_MPI_barrier", "no", "barrier in timer [yes,no]");

  AddPar("nan_check", "", "variables we check for NANs [varnames]");
  AddPar("sysmon_hours", "0.5", "hours after which we write sysmon.log");
  AddPar("sysmon_output_per_rank", "no","output data for each rank [no,yes]");
  AddPar("update_parameters_hours", "0.2",
         "hours after which we read nmesh_update_parameters.par");

  AddPar("fwrite_bufsize","1048576", "use setvbuf with this size for write");
  AddPar("fread_bufsize", "0",       "setvbuf size for read, 0 means we "
         "don't call setvbuf so that we use C's default buffering");
  AddPar("fs_sync", "no", "whether fs_sync is active [yes,no]");
  AddPar("file_sync", "no", "how we sync a file [fdatasync,fsync,sync]");

  return 0;
}
