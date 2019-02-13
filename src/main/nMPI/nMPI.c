/* nMPI.c */
/* Wolfgang Tichy, 1/2019 */

#include "nmesh.h"
#include "nMPI.h"



/* print some compile info */
int nMPI_print_compile_info(tMesh *mesh)
{
#ifdef USEMPI
  printf("MPI is compiled in.\n");
#else
  printf("MPI is not compiled in.\n");
#endif
  return 0;
}

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


/* exchange double buffers */
void nMPI_Isend_Irecv_double(double *sbuf, int ns, double *rbuf, int nr,
                             int rank_other, int s_tag, int r_tag,
                             nMPI_Req *s_req, nMPI_Req *r_req)
{
#ifdef USEMPI
  int errS, errR;

  PRF;printf(": %d to %d, ns=%d nr=%d s_tag=%d r_tag=%d\n",
             nMPI_rank(), rank_other, ns, nr, s_tag, r_tag);
  fflush(stdout);

  errS = MPI_Isend(sbuf, ns, MPI_DOUBLE, rank_other, s_tag,
                   MPI_COMM_WORLD, s_req);
  if(errS != MPI_SUCCESS) errorexit("MPI_Isend failed!\n");
  
  errR = MPI_Irecv(rbuf, nr, MPI_DOUBLE, rank_other, r_tag,
                   MPI_COMM_WORLD, r_req);
  if(errR != MPI_SUCCESS) errorexit("MPI_Irecv failed!\n");
#endif
}
