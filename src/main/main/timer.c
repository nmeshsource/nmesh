/* timer.c */
/* Wolfgang Tichy 3/2019 */


#include "nmesh.h"


extern tMesh *main_mesh;


/* timer data base */
typedef struct tTIMER {
  char *name;    /* name of function we time */
  double start;  /* time when timer was started */
  double time;   /* time in s spent in this func */
  long n;        /* number of times func was timed */
} tTimer;


/* this is a case where global variables are a good idea */
static tTimer **tdb = NULL; /* list of timers */
int ntimers = 0;

static int timer_on = -1;          /* -1 marks that no timer was used yet */
static int timer_MPI_barrier = 0;




/* get or make a new timer for the function named "funcname",
   returns pointer to entry for name */
tTimer *timer_get(const char *funcname)
{
  int i;
  int tid = TID;
  int namelen = strlen(funcname) + 64;
  char *name = malloc(sizeof(name[0])*namelen);

  /* append thread ID to name, if we do not call from thread0 */
  if(tid) snprintf(name,namelen, "%s %d", funcname, TID);
  else    snprintf(name,namelen, "%s", funcname);

  GEN_Pragma(omp critical (timer_get))
  {
    /* look for timer with correct name */
    for(i = 0; i < ntimers; i++)
      if(strcmp(tdb[i]->name, name) == 0)
        break;

    /* if we cannot find timer */
    if(i >= ntimers)
    {
      /* make a new timer */
      i = ntimers;
      ntimers++;
      /* realloc timer list */
      tdb = realloc(tdb, ntimers * sizeof(tdb[0]));
      /* alloc and set new timer */
      tdb[i] = calloc(1, sizeof(*tdb[0]));
      tdb[i]->name  = name;
      tdb[i]->start = -1.;
      tdb[i]->time  = 0.;
      tdb[i]->n     = 0;
    }
    else /* timer has index i and now we do not need name anymore */
    {
      free(name);
    }
  }

  return tdb[i];
}


/* start a timer */
tTimer *timer_start(const char *name)
{
  tTimer *t = 0;

  if(timer_on)
  {
    /* is this the first time we start a timer? */
    if(timer_on == -1)
    {
      tMesh *mesh = main_mesh;
      timer_on = Getv(Par("timer_on"), "yes");
      timer_MPI_barrier = Getv(Par("timer_MPI_barrier"), "yes");
    }

    /* get timer and save current time in it */
    t = timer_get(name);
    if(timer_MPI_barrier) nMPI_barrier();
    t->start = getTimeIn_s();
    t->n += 1; /* func now has been called */
  }
  return t;
}

/* stop timer */
tTimer *timer_stop(const char *name)
{
  tTimer *t = 0;

  if(timer_on)
  {
    t = timer_get(name);

    /* t->start = -1 marks that timer is stopped already */
    if(t->start < 0.) return t;

    if(timer_MPI_barrier) nMPI_barrier();

    /* save timing info */
    t->time += getTimeIn_s() - t->start;

    /* t->start = -1 marks that timer is stopped now */
    t->start = -1.;
  }
  return t;
}

/* update time field of timer without stopping it */
tTimer *timer_update(const char *name)
{
  tTimer *t = 0;

  if(timer_on)
  {
    timer_stop(name);
    t = timer_start(name);
    t->n -= 1;
  }
  return t;
}

/* read out timer without stopping it */
double timer_read(const char *name)
{
  if(timer_on)
  {
    tTimer *t = timer_get(name);

    /* t->start = -1 marks that timer is stopped already */
    if(t->start < 0) return -1.;

    return getTimeIn_s() - t->start + t->time;
  }
  return -1.;
}


/* free all timers */
int free_all_timers(tMesh *mesh)
{
  int i;
  for(i = 0; i < ntimers; i++)
  {
    free(tdb[i]->name);
    free(tdb[i]);
  }
  free(tdb);
  tdb = NULL;
  return 0;
}


/* compare times for qsort */
int timer_compar(const void *x1, const void *x2)
{
  tTimer *const *t1 = x1;
  tTimer *const *t2 = x2;

  if(t1[0]->time < t2[0]->time) return 1;
  if(t1[0]->time > t2[0]->time) return -1;
  return 0;
}


/* write timer info into files */
int write_all_timers(tMesh *mesh)
{
  char *outdir = Gets(Par("outdir"));
  char f[100], s[1000];
  FILE *fp;
  tTimer *t;
  double pct, total;
  int i;

  if (!timer_on) return 0;

  /* read timer for main() */
  t = timer_update("main");
  total = t->time;

  /* open file */
  snprintf(f, 100, "%%s/timer.%%0%dd", (int) log10(nMPI_size())+1);
  snprintf(s, 1000, f, outdir, nMPI_rank());
  fp = fopen(s, "a");
  if(!fp) errorexits("could not open %s", s);
  fprintf(fp, "---------------------------------------------------------"
          "---------------------\n");
  fprintf(fp, "iteration %d, time %g\n", mesh->iteration, mesh->time);
  fprintf(fp, "Functions with timer                                %%"
          "      time/s       calls\n");
  fprintf(fp, "---------------------------------------------------------"
          "---------------------\n");

  /* sort according to time spent */
  qsort(tdb, ntimers, sizeof(tdb[0]), timer_compar);

  /* output */
  for(i = 0; i < ntimers; i++)
  {
    t   = tdb[i];
    pct = 100. * t->time/total;
    fprintf(fp, "%-49s %6.2f %9.2e %10ld\n", t->name, pct, t->time, t->n);
  }
  fprintf(fp, "\n\n");
  fclose(fp);
  return 0;
}
