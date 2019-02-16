/* nMPI_defs.h */
/* (c) Wolfgang Tichy 2/2019 */
/* header file for global defs */

/* some MPI things we use in many places */
#ifdef USEMPI

#include <mpi.h>
#define WORLD MPI_COMM_WORLD
#define nMPI_Comm MPI_Comm
#define nMPI_Req  MPI_Request
#define nMPI_Stat MPI_Status

#else

#define WORLD 0
#define nMPI_Comm int
#define nMPI_Req  int
#define nMPI_Stat int

#endif
