/* nmesh_nMPI.h */
/* (c) Wolfgang Tichy 1/2019 */
/* header file for global functions */


/* marco to do: if(Rank0) */
#define Rank0 (!nMPI_rank())

/* for MPI communication this can help */
typedef struct tCOM {
  int n_rq;             /* number of send/recv requests */
  int send_i;           /* counter, can be used in loops */
  int recv_i;           /* counter, can be used in loops */
  nMPI_Req *send_rq;    /* list of send requests */
  nMPI_Req *recv_rq;    /* recv req: recv_rq[2] is req. 2 */
  nMPI_Stat *send_stat; /* status array for sends */
  nMPI_Stat *recv_stat; /* status array recvs*/
  int *send_buflen;     /* length of send buffer */
  int *recv_buflen;     /* length of recv buffer */
  void **send_buf;      /* send buffer */
  void **recv_buf;      /* recv buffer */
  int entrysize;        /* size of entries in send_buf and recv_buf */
  int free_buf;         /* free buffers in free_com and realloc_com_reqs */
} tCom;



/* nMPI.c */
int nMPI_long_tag_to_commi_tag(long long_tag, int *commi, int *tag);
int nMPI_print_compile_info(tMesh *mesh);
int nMPIvars_init(tMesh *mesh);
int nMPIvars_finalize(tMesh *mesh);
nMPI_Comm nMPIvars_get_comm(int i);
int nMPIvars_get_ncomms(void);
int nMPIvars_get_ntags_limit(void);
int nMPI_Init(int *pargc, char ***pargv);
int nMPI_Finalize(void);
int nMPI_Abort(int errorcode);
int nMPI_rank(void);
int nMPI_size(void);
int nMPI_barrier(void);
int nMPI_Comm_get_attr(nMPI_Comm comm, int comm_keyval,
                       void *attribute_val, int *flag);
int nMPI_Comm_dup(nMPI_Comm comm, nMPI_Comm *newcomm);
int nMPI_Comm_free(nMPI_Comm *comm);
int nMPI_Send(const void *buf, int count, nMPI_Datatype datatype,
              int dest, int tag);
int nMPI_Recv(void *buf, int count, nMPI_Datatype datatype,
              int source, int tag);
int nMPI_Waitall(int nreq, nMPI_Req *req, nMPI_Stat *stat);
int nMPI_Allreduce(const void *sendbuf, void *recvbuf, int count,
                   nMPI_Datatype datatype, nMPI_Op op);
int nMPI_Reduce(const void *sendbuf, void *recvbuf, int count,
                nMPI_Datatype datatype, nMPI_Op op, int root);
int nMPI_Bcast(void *buffer, int count, nMPI_Datatype datatype, int root);
int nMPI_Ibcast(void *buffer, int count, nMPI_Datatype datatype,
                int root, nMPI_Req *request);
int nMPI_Test(nMPI_Req *request, int *flag, nMPI_Stat *stat);
int nMPI_Wait(nMPI_Req *req, nMPI_Stat *stat);

tCom *alloc_com(int entrysize, int free_buf);
void free_com(tCom *com);
void realloc_com_reqs(tCom *com, int n_rq_new);
void print_com(tCom *com);
void set_free_buf_in_com(tCom *com, int free_buf);
void put_buffers_in_com(tCom *com, int rq,
                        void *sbuf, int slen, void *rbuf, int rlen);
int append_buffers_to_com(tCom *com, void *sbuf,int slen, void *rbuf,int rlen);
void *get_com_send_buf(tCom *com, int rq);
void *get_com_recv_buf(tCom *com, int rq);
void set_com_send_buf(tCom *com, int rq, void *send_buf);
void set_com_recv_buf(tCom *com, int rq, void *recv_buf);
void set_com_counters(tCom *com, int si, int ri);
void inc_com_send_i(tCom *com);
void inc_com_recv_i(tCom *com);
void *get_com_send_i_buf_inc_i(tCom *com);
void *get_com_recv_i_buf_inc_i(tCom *com);
void *get_com_send_i_buf(tCom *com);
void *get_com_recv_i_buf(tCom *com);
void free_com_send_i_buf(tCom *com);
void free_com_recv_i_buf(tCom *com);
int nMPI_Waitall_com_send(tCom *com);
int nMPI_Waitall_com_recv(tCom *com);
int nMPI_Waitall_com(tCom *com);
int nMPI_Wait_com_send(tCom *com, int rq);
int nMPI_Wait_com_recv(tCom *com, int rq);
int nMPI_Test_com_send(tCom *com, int rq, int *flag);
int nMPI_Test_com_recv(tCom *com, int rq, int *flag);
int nMPI_Testall_com_send(tCom *com, int *flag);
int nMPI_Testall_com_recv(tCom *com, int *flag);
int nMPI_Testall_com(tCom *com, int *flag);
void nMPI_Isend_Irecv_double_com(tCom *com, int rq,
                                 int rank_other, int s_tag, int r_tag,
                                 nMPI_Comm s_comm, nMPI_Comm r_comm);
void nMPI_Isend_double_com(tCom *com, int rq, int dest, int tag,nMPI_Comm comm);
void nMPI_Irecv_double_com(tCom *com, int rq, int src, int tag,nMPI_Comm comm);



/**************************************************************************/
/* insert lists.h here because nMPI is the 2nd module,
   so all after will see the lists types and such */
/**************************************************************************/
#include "../main/lists.h"
