/* nMPI_defs.h */
/* (c) Wolfgang Tichy 2/2019 */
/* header file for global defs */

/* some MPI things we use in many places */
#define WORLD main_comm

/* redefine some things */
#ifdef USEMPI

#include <mpi.h>
#define nMPI_COMM_WORLD MPI_COMM_WORLD
#define nMPI_Comm MPI_Comm
#define nMPI_Req  MPI_Request
#define nMPI_Stat MPI_Status

#else

#define nMPI_COMM_WORLD 0
#define nMPI_Comm int
#define nMPI_Req  int
#define nMPI_Stat int

#endif
