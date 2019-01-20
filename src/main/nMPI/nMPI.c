/* mpi.c */
/* Wolfgang Tichy, 1/2019 */

#include "nmesh.h"
#include "nMPI.h"

/* include files for MPI */
#ifdef USEMPI
#include <mpi.h>
#endif


/* Wrappers for MPI_Init and MPI_Finalize */
int nMPI_Init(int *pargc, char ***pargv)
{
#ifdef USEMPI
  return MPI_Init(pargc, pargv);
#endif
  return 0;
}
int nMPI_Finalize(void)
{
#ifdef USEMPI
  return MPI_Finalize();
#endif
  return 0;
}

/* return MPI rank, if MPI is not compiled in return 0 */
int nMPI_rank(void)
{
  int rank=0;
#ifdef USEMPI
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
#endif
  return rank;
}

/* number of processes we are running with, or 1 if no MPI */
int nMPI_size(void)
{
  int size=1;
#ifdef USEMPI
  MPI_Comm_size(MPI_COMM_WORLD, &size);
#endif
  return size;
}

/* barrier */
int nMPI_barrier(void)
{
  int ret=0;
#ifdef USEMPI
  ret = MPI_Barrier(MPI_COMM_WORLD);
#endif
  return ret;
}
