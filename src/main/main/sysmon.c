/* sysmon.c */
/* Wolfgang Tichy, Feb. 2020 */

#include "nmesh.h"

#include <sys/time.h>      /* for getrusage */
#include <sys/resource.h>  /* for getrusage */


/* functions needed */
void write_sysmon(tMesh *mesh, double last_mesh_time, const char *name,
                  double time, double last_sysmon_time, double *data,
                  int addheader);



/* write into sysmon.log if the time is right for it */
int sysmon(tMesh *mesh)
{
  static int firstcall = 1;
  static double last_sysmon_time = 0.;
  static double last_mesh_time = 0.;
  double time  = getTimeIn_s()/3600.;
  double time_since_sysmon;
  double hours = Getd(Par("sysmon_hours"));
  int do_sysmon = 0;
  int output_per_rank = Getb(Par("sysmon_output_per_rank"));

  ///* is sysmon on? */
  //if(!Getb(Par("sysmon"))) return 0;

  /* Reset firstcall if mesh->time seems to have gone backwards.
     This can happen if checkpoint=yes and iterate_parameters=yes causes
     the loading of more than one checkpoint for one of the later pars in
     iterate_parameter1. Then we can get a negative dtPhys/dtWall. */
  if( (mesh->time < last_mesh_time) || mesh->time==0. )
    firstcall = 1;

  /* test if it is time */
  time_since_sysmon = time - last_sysmon_time;
  if(Rank0)
  {
    /* test based on walltime */
    if(hours >= 0. && (hours <= time_since_sysmon || firstcall))
      do_sysmon = 1; /* yes, we want to write sysmon data */
  }

  /* broadcast do_sysmon from rank0 to all others */
  MCK( nMPI_Bcast(&do_sysmon, 1, nMPI_INT, 0) );

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
    MCK( nMPI_Allreduce(dat, datall, n, nMPI_DOUBLE, nMPI_SUM) );

    /* output datall on rank 0 */
    if(Rank0)
    {
      write_sysmon(mesh, last_mesh_time, "sysmon.log",
                   time, last_sysmon_time, datall, firstcall);
    }

    /* output dat if needed */
    if(output_per_rank)
    {
      char form[100];
      char name[100];
      /* filename format and file name for files */
      snprintf(form,99, "sysmon.%%0%dd", (int) log10(nMPI_size())+1);
      snprintf(name,99, form, nMPI_rank());
      write_sysmon(mesh, last_mesh_time, name,
                   time, last_sysmon_time, dat, firstcall);
    }

    /* update times */
    last_sysmon_time = time;
    last_mesh_time = mesh->time;
    firstcall = 0;
  }
  return 0;
}

/* write sysmon results (in data) into file outdir/name */
void write_sysmon(tMesh *mesh, double last_mesh_time, const char *name,
                  double time, double last_sysmon_time, double *data,
                  int addheader)
{
  char fname[8192];
  FILE *fp;

  /* open destination file */
  snprintf(fname,8192, "%s/%s", Gets(Par("outdir")), name);
  fp = fopen(fname, "a");
  if(fp)
  {
    double dpt = mesh->time - last_mesh_time;
    double dt  = time       - last_sysmon_time;

    if(addheader)
      fprintf(fp, "#    PhysTime        WallTime"
                  "    dtPhys/dtWall         max(RSS)\n");

    fprintf(fp, "%13g  %13gh  %13g/h  %13gGi\n",
            mesh->time, time, dpt/dt, data[0]/1048576);
    fclose(fp);
  }
}


/* ignore the par sysmon_hours and write into sysmon.log now */
int sysmon_now(tMesh *mesh)
{
  char *sysmon_hours_bak = strdup(Gets(Par("sysmon_hours")));
  double time;

  /* signal that we want sysmon right now */
  Sets(Par("sysmon_hours"), "0");
  time = getTimeIn_s()/3600.;
  sysmon(mesh);
  PRF;printf(": at WallTime = %13gh\n", time);

  /* restore par sysmon_hours */
  Sets(Par("sysmon_hours"), sysmon_hours_bak);
  free(sysmon_hours_bak);

  return 0;
}
