/* nMPI.h */
/* (c) Wolfgang Tichy 1/2019 */
/* header file for nMPI local functions */


/* structure that holds global nMPI vars */
typedef struct {
  int comm_bits;     /* number of bits we use for MPI communicators */
  int ncomms;        /* number of MPI communicators: ncomms = 2^comm_bits */
  nMPI_Comm *comm;   /* list of MPI communicators */
  int tag_ub;        /* min of upper bound of MPI tags over all in comm */
  int tag_bits;      /* number of bits in tag_ub+1 */
  nMPI_Datatype TELM;  /* MPI_Datatype for tElm */
} tnMPIvars;


/* world comm from main */
extern nMPI_Comm main_comm;


/* nMPI.c */
