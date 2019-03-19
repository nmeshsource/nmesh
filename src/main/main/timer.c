/* timer.c */
/* Wolfgang Tichy 3/2019 */


#include "nmesh.h"

#define PR 0

extern tMesh *main_mesh;


/* timer data base */
typedef struct {
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




/* get or make a new timer for the function named "name",
   returns pointer to entry for name */
tTimer *timer_get(char *name)
{
  int i;

  if(tdb != NULL)
  {
    /* return timer with correct name */
    for(i = 0; i < ntimers; i++)
      if(strcmp(tdb[i]->name, name) == 0)
        return tdb[i];
  }

  /* make a new timer */
  i = ntimers;
  ntimers++;
  /* realloc timer list */
  tdb = realloc(tdb, ntimers * sizeof(tdb[0]));
  /* alloc and set new timer */
  tdb[i] = calloc(1, sizeof(*tdb[0]));
  tdb[i]->name  = strdup(name);
  tdb[i]->start = -1;
  tdb[i]->time  = 0;
  tdb[i]->n     = 0;

  return tdb[i];
}


/* start a timer */
void timer_start(char *name)
{
  if(timer_on)
  {
    tTimer *t;

    /* is this the first time we start a timer? */
    if(timer_on == -1)
    {
      tMesh *mesh = main_mesh;
      timer_on = Getv(Par("timer_on"), "yes");
      timer_MPI_barrier = Getv(Par("timer_MPI_barrier"), "yes");
    }

    t = timer_get(name);
    if(timer_MPI_barrier) nMPI_barrier();
    t->start = getTimeIn_s();
  }
}

/* stop timer */
void timer_stop(char *name)
{
  if(timer_on)
  {
    tTimer *t = timer_get(name);

    /* t->start = -1 marks that timer is stopped already */
    if(t->start < 0.) return;

    if(timer_MPI_barrier) nMPI_barrier();

    /* save timing info */
    t->time += getTimeIn_s() - t->start;
    t->n    += 1;

    /* t->start = -1 marks that timer is stopped now */
    t->start = -1.;
  }
}

/* read out timer without stopping it */
double timer_read(char *name)
{
  if(timer_on)
  {
    tTimer *t = timer_get(name);

    /* t->start = -1 marks that timer is stopped already */
    if(t->start < 0) return -1.;

    return getTimeIn_s() - t->start;
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


/* compare times for qsort_r */
int timer_compar(const void *x1, const void *x2)
{
  tTimer *const *t1 = x1;
  tTimer *const *t2 = x2;

  if(t1[0]->time < t2[0]->time) return -1;
  if(t1[0]->time > t2[0]->time) return  1;
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
  total = timer_read("main");

  /* open file */
  snprintf(f, 100, "%%s/timer.%%0%dd", (int) log10(nMPI_size())+1);
  snprintf(s, 1000, f, outdir, nMPI_rank());
  fp = fopen(s, "a");
  if(!fp) errorexits("could not open %s", s);
  fprintf(fp, "Timers at iteration %d, time %g\n",
          mesh->iteration, mesh->time);

  /* sort according to time spent */
  qsort(tdb, ntimers, sizeof(tdb[0]), timer_compar);

  /* output */
  for(i = 0; i < ntimers; i++)
  {
    t   = tdb[i];
    pct = 100. * t->time/total;
    fprintf(fp, "%-24s %5.1f  %10.3f %7ld\n", t->name, pct, t->time, t->n);
  }
  fprintf(fp,"\n");
  fclose(fp);
  return 0;
}
