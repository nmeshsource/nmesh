/* nMPI_defs.h */
/* (c) Wolfgang Tichy 2/2019 */
/* header file for global defs */

/* some MPI things we use in many places */
#ifdef USEMPI

#include <mpi.h>
#define nMPI_Req  MPI_Request
#define nMPI_Stat MPI_Status

#else

#define nMPI_Req  int
#define nMPI_Stat int

#endif
