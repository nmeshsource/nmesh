/* nMPI.h */
/* (c) Wolfgang Tichy 1/2019 */
/* header file for nMPI local functions */


/* structure that holds global nMPI vars */
#define NCOMMS 256
typedef struct {
  nMPI_Comm comm[NCOMMS];
} tnMPIvars;


/* world comm from main */
extern nMPI_Comm main_comm;


/* nMPI.c */
