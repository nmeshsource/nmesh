/* timing.c */
/* Wolfgang Tichy, 12/2020 */

#include "nmesh.h"
#include "amr.h"


tTiming Timing[1];



/* print global Timing struct */
void printTiming(void)
{
  printf("Timing->\n");
  printf("  mm1_speed = %g\n", Timing->mm1_speed);
  printf("  mm_speed = %g\n", Timing->mm_speed);
  printf("  myops = %g\n", Timing->myops);
  printf("  ops0 = %g\n", Timing->ops0);
  printf("  allops = %g\n", Timing->allops);
}


/* get time for one matrix mul. */
double time_mm_array0(tArray *At, tArray *B, tArray *AB)
{
  struct timespec tp0[1];
  struct timespec tp1[1];

  getRealTime(tp0);
  mm_array0(At, B, AB);
  getRealTime(tp1);

  return getTimeDiffIn_s(tp1, tp0);
}

/* do many matrix multiplications and distriubute them among
   OpenMP threads to determine the proc speed */
int timing_mm_speed(tMesh *mesh)
{
  struct timespec tp0[1];
  struct timespec tp1[1];
  int i;
  int nxmax = 80;
  int nruns = 100;
  int n[] = {nxmax,nxmax,nxmax};
  int nx = nxmax;
  int ny = nx-1;
  int nz = nx-2;
  double speednorm = 1.;
  double mm1_time;

  /* time before loop */
  getRealTime(tp0);

  /* time nruns times to have less jitter */
  mm1_time = 0.;
  NODELEVEL_Pragma(omp parallel for reduction(+:mm1_time))
  for(i=0; i<nruns; i++)
  {
    int k;
    tArray *At, *B, *AB;
    struct timespec t0[1];
    struct timespec t1[1];

    /* time before any work */
    getRealTime(t0);

    At = alloc_array2d(n[0], n[0]);
    B  = alloc_array(n);
    AB = alloc_array(n);

    /* put some numbers in At, B, AB */
    forarray(At, k) Arrd(At)[k] = k+1;
    forarray(B,  k) Arrd(B)[k]  = k*0.3;
    forarray(AB, k) Arrd(AB)[k] = k*9;

    nx = nxmax;
    redim_array(At, nx, nx, 1);
    redim_array(B,  nx, ny, nz);
    redim_array(AB, nx, ny, nz);

    //mm1_time += time_mm_array0(At, B, AB);
    mm_array0(At, B, AB);

    free_array(AB);
    free_array(B);
    free_array(At);

    /* time after work */
    getRealTime(t1);
    mm1_time += getTimeDiffIn_s(t1, t0);
  }

  /* time after loop */
  getRealTime(tp1);

  /* time and speed for one iteration of the loop */
  mm1_time /= nruns;
  Timing->mm1_speed = speednorm / mm1_time;
  PRFs(":\n");
  printf("  speed for 1 mm_array0 for %dx%d * %dx%d:    mm1_speed = %g\n",
         nx,nx, nx,ny*nz, Timing->mm1_speed);

  /* speed with which entire loop was done */
  Timing->mm_speed = speednorm / (getTimeDiffIn_s(tp1, tp0)/nruns);
  printf("  speed for %d mm_array0 using %4d thread(s):  mm_speed = %g\n",
         nruns, MAX_NTHREADS, Timing->mm_speed);

  return 0;
}

/* return Timing->mm_speed with non-zero floor */
double timing_get_mm_speed(tMesh *mesh)
{
  double myspeed  = Timing->mm_speed;
  double speedmin = 1e-50;

  /* in case we forgot to measure Timing->mm_speed, just set myspeed=1 */
  if(myspeed <= speedmin) myspeed = 1.;
  return myspeed;
}

/* get time from dat->info->load_TimeIn_s with non-zero floor */
double timing_get_elm_load_TimeIn_s(tElm *elm)
{
  double loadTmin = 1e-50;
  double et;
  tDat *dat = elm->dat;

  if(dat)
  {
    int ijk = elm_get_ijk(elm);
    double tw;

    /* read dat->info->load_TimeIn_s */
    et = dat->info->load_TimeIn_s;
    /* if we forgot to measure load_TimeIn_s of elm, just set et=loadTmin */
    if(et <= loadTmin) et = loadTmin;

    /* set timing weight for this elm */
    if(ijk==0) tw = 1.;
    else       tw = Timing->child1to7_weight;

    et = et * tw;
    return et;
  }
  else
  {
    return 0.;
  }
}

/* set number of operations myops that were done on this rank */
int timing_set_myops(tMesh *mesh)
{
  double myspeed = timing_get_mm_speed(mesh);
  struct list_head *pos;
  double myT = 0.;

  list_for_each(pos, &mesh->myelm_head)
  {
    tElm *elm = list_entry(pos, tElm, list);
    double et = timing_get_elm_load_TimeIn_s(elm);
    myT += et;
  }
  Timing->myops = myspeed * myT;
  return 0;
}

/* set total number of operations ops0 done on all ranks below current one */
int timing_set_myops_ops0_allops(tMesh *mesh)
{
  int size = nMPI_size();
  int rank = nMPI_rank();
  double ops1;

  /* set myops first */
  timing_set_myops(mesh);

  /* receive ops0 from rank-1, unless we are rank0 */
  if(rank > 0)
  {
    /* we use blocking MPI here */
    nMPI_Recv(&(Timing->ops0),1, nMPI_DOUBLE, rank-1, 123);
    /* This blocks until we actually get ops0 from previous rank.
       We do not want to go any further until we have ops0. */
  }

  /* send ops0+myops to rank+1 */
  ops1 = Timing->ops0 + Timing->myops;
  if(rank < size-1)
  {
    /* we use blocking MPI here */
    nMPI_Send(&(ops1),1, nMPI_DOUBLE, rank+1, 123);
    /* this blocks until ops1 is in some network buffer */
  }

  /* last rank knows total number of operations allops, so broadcast it */
  Timing->allops = ops1;
  nMPI_Bcast(&(Timing->allops),1, nMPI_DOUBLE, size-1);

  return 0;
}
