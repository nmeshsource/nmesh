/* nMPI_defs.h */
/* (c) Wolfgang Tichy 2/2019 */
/* header file for global defs */

/* some MPI things we use in many places */
#define WORLD main_comm /* use this, never nMPI_COMM_WORLD */

/* redefine some things */
#ifdef USEMPI

#include <mpi.h>
#define nMPI_COMM_WORLD MPI_COMM_WORLD
#define nMPI_Comm MPI_Comm
#define nMPI_Req  MPI_Request
#define nMPI_Stat MPI_Status
#define nMPI_Op   MPI_Op
#define nMPI_MAX  MPI_MAX
#define nMPI_MIN  MPI_MIN
#define nMPI_SUM  MPI_SUM
#define nMPI_Datatype MPI_Datatype
#define nMPI_INT      MPI_INT
#define nMPI_DOUBLE   MPI_DOUBLE

#else

#define nMPI_COMM_WORLD 0
#define nMPI_Comm int
#define nMPI_Req  int
#define nMPI_Stat int
#define nMPI_Op   int
#define nMPI_MAX  0
#define nMPI_MIN  0
#define nMPI_SUM  0
#define nMPI_Datatype int
#define nMPI_INT      0
#define nMPI_DOUBLE   0

#endif
