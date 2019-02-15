/* nmesh_nMPI.h */
/* (c) Wolfgang Tichy 1/2019 */
/* header file for global functions */


/* marco to do: if(Rank0) */
#define Rank0 (!nMPI_rank())

/* for MPI communication this can help */
typedef struct tCOM {
  int n_rq;             /* number of send/recv requests on each face */
  nMPI_Req *send_rq;    /* send request array */
  nMPI_Req *recv_rq;    /* recv req: recv_rq[2] is req. 2 */
  nMPI_Stat *send_stat; /* status array for sends */
  nMPI_Stat *recv_stat; /* status array recvs*/
  int *send_buflen;     /* length of send buffer */
  int *recv_buflen;     /* length of recv buffer */
  void **send_buf;      /* send buffer */
  void **recv_buf;      /* recv buffer */
  int entrysize;        /* size of entries in send_buf and recv_buf */
  int free_buf;         /* do we free buffers in free_com and realloc_com_reqs */
} tCom;



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
int nMPI_Waitall(int nreq, nMPI_Req *req, nMPI_Stat *stat);
tCom *alloc_com(int entrysize, int free_buf);
void free_com(tCom *com);
void realloc_com_reqs(tCom *com, int n_rq_new);
void print_com(tCom *com);
void put_buffers_in_com(tCom *com, int rq,
                        void *sbuf, int slen, void *rbuf, int rlen);
int append_buffers_to_com(tCom *com, void *sbuf,int slen, void *rbuf,int rlen);
void *get_com_send_buf(tCom *com, int rq);
void *get_com_recv_buf(tCom *com, int rq);
int nMPI_Waitall_com_send(tCom *com);
int nMPI_Waitall_com_recv(tCom *com);
int nMPI_Waitall_com(tCom *com);
void nMPI_Isend_Irecv_double_com(tCom *com, int rq,
                                 int rank_other, int s_tag, int r_tag);
void nMPI_Isend_double_com(tCom *com, int rq, int dest, int tag);
void nMPI_Irecv_double_com(tCom *com, int rq, int src, int tag);
