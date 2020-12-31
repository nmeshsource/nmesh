/* timing.c */
/* Wolfgang Tichy, 12/2020 */

#include "nmesh.h"
#include "amr.h"


tTiming Timing[1];


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
