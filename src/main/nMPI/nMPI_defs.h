/* nMPI_defs.h */
/* (c) Wolfgang Tichy 2/2019 */
/* header file for global defs */

/* some MPI things we use in many places */
#define WORLD main_comm /* use this, never nMPI_COMM_WORLD */

/* redefine some things */
#ifdef USEMPI

#define Rank0_or_NoMPI Rank0
#include <mpi.h>
#define nMPI_COMM_WORLD MPI_COMM_WORLD
#define nMPI_COMM_NULL  MPI_COMM_NULL
#define nMPI_Comm MPI_Comm
#define nMPI_Req  MPI_Request
#define nMPI_Stat MPI_Status
#define nMPI_STATUS_IGNORE MPI_STATUS_IGNORE
#define nMPI_Op   MPI_Op
#define nMPI_MAX  MPI_MAX
#define nMPI_MIN  MPI_MIN
#define nMPI_SUM  MPI_SUM
#define nMPI_LAND MPI_LAND
#define nMPI_LOR  MPI_LOR
#define nMPI_Datatype MPI_Datatype
#define nMPI_CHAR     MPI_CHAR
#define nMPI_INT      MPI_INT
#define nMPI_DOUBLE   MPI_DOUBLE
#define nMPI_LONG     MPI_LONG
#define nMPI_MAXLOC     MPI_MAXLOC
#define nMPI_MINLOC     MPI_MINLOC
#define nMPI_DOUBLE_INT MPI_DOUBLE_INT
#define nMPI_TAG_UB MPI_TAG_UB

#else

#define Rank0_or_NoMPI 1
#define nMPI_COMM_WORLD 0
#define nMPI_COMM_NULL  0
#define nMPI_Comm int
#define nMPI_Req  int
#define nMPI_Stat int
#define nMPI_STATUS_IGNORE 0
#define nMPI_Op   int
#define nMPI_MAX  0
#define nMPI_MIN  0
#define nMPI_SUM  0
#define nMPI_LAND 0
#define nMPI_LOR  0
#define nMPI_Datatype int
#define nMPI_CHAR     0
#define nMPI_INT      0
#define nMPI_DOUBLE   0
#define nMPI_LONG     0
#define nMPI_MAXLOC     0
#define nMPI_MINLOC     0
#define nMPI_DOUBLE_INT 0
#define nMPI_TAG_UB 0

#endif
