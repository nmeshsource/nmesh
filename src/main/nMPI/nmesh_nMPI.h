/* nmesh_nMPI.h */
/* (c) Wolfgang Tichy 1/2019 */
/* header file for global functions */


/* nMPI.c */
int nMPI_print_compile_info(tMesh *mesh);
int nMPI_Init(int *pargc, char ***pargv);
int nMPI_Finalize(void);
int nMPI_rank(void);
int nMPI_size(void);
int nMPI_barrier(void);
void nMPI_Isend_Irecv_double(double *sbuf, int ns, double *rbuf, int nr,
                             int rank_other, int s_tag, int r_tag,
                             nMPI_Req *s_req, nMPI_Req *r_req);
int nMPI_Waitall(int nreq, nMPI_Req *req, nMPI_Stat *stat) ;
