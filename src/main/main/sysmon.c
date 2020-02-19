/* sysmon.c */
/* Wolfgang Tichy, Feb. 2020 */

#include "nmesh.h"

#include <sys/time.h>      /* for getrusage */
#include <sys/resource.h>  /* for getrusage */



/* write into sysmon.log if the time is right for it */
int sysmon(tMesh *mesh)
{
  static int firstcall = 1;
  static double last_sysmon_time = 0.;
  static double last_mesh_time = 0.;
  double hours = Getd(Par("sysmon_hours"));
  double time  = getTimeIn_s()/3600.;
  double time_since_sysmon;
  int do_sysmon = 0;

  ///* is sysmon on? */
  //if(!Getb(Par("sysmon"))) return 0;

  /* reset firstcall & last_mesh_time if we are again at interation 0 */
  if(mesh->iteration==0 && firstcall==0)
  {
    firstcall = 1;
    last_mesh_time = mesh->time;
  }

  /* test if it is time */
  time_since_sysmon = time - last_sysmon_time;
  if(Rank0)
  {
    /* test based on walltime */
    if((hours >= 0. && hours <= time_since_sysmon) || firstcall)
      do_sysmon = 1; /* yes, we want to write sysmon data */
  }

  /* broadcast do_sysmon from rank0 to all others */
  nMPI_Bcast(&do_sysmon, 1, nMPI_INT, 0);

  /* now do it if needed */
  if(do_sysmon)
  {
    int i, n;
    int count = 1;
    double dat[count];
    double datall[count];
    struct rusage usage[1];
    long ru_maxrss = 0;

    /* get info */
    getrusage(RUSAGE_SELF, usage);
#ifndef NO_RU_MAXRSS
    ru_maxrss = usage->ru_maxrss;  /* struct rusage has ru_maxrss */
#endif
    n = 0;
    dat[n++] = ru_maxrss;  /* maximum resident set size */
    for(i=0; i<n; i++) datall[i] = dat[i]; /* in case MPI is not there */
    nMPI_Allreduce(dat, datall, n, nMPI_DOUBLE, nMPI_SUM);

    if(Rank0)
    {
      char fname[8192];
      FILE *fp;

      /* open destination file */
      snprintf(fname,8192, "%s/%s", Gets(Par("outdir")), "sysmon.log");
      fp = fopen(fname, "a");
      if(fp)
      {
        double dpt = mesh->time - last_mesh_time;
        double dt  = time       - last_sysmon_time;

        if(firstcall)
          fprintf(fp, "#    PhysTime        WallTime"
                      "    dtPhys/dtWall         max(RSS)\n");

        fprintf(fp, "%13g  %13gh  %13g/h  %13gGi\n",
                mesh->time, time, dpt/dt, datall[0]/1048576);
        fclose(fp);
      }
    }

    /* update times */
    last_sysmon_time = time;
    last_mesh_time = mesh->time;
    firstcall = 0;
  }
  return 0;
}
