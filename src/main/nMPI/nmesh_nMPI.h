/* nmesh_nMPI.h */
/* (c) Wolfgang Tichy 1/2019 */
/* header file for global functions */


/* nMPI.c */
int nMPI_Init(int *pargc, char ***pargv);
int nMPI_Finalize(void);
int nMPI_rank(void);
int nMPI_size(void);
int nMPI_barrier(void);
