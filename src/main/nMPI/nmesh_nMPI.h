/* nmesh_nMPI.h */
/* (c) Wolfgang Tichy 1/2019 */
/* header file for global functions */


/* marco to do: if(Rank0) */
#define Rank0 (!nMPI_rank())


/* structure that holds global nMPI vars */
typedef struct {
  int comm_bits;     /* number of bits we use for MPI communicators */
  int ncomms;        /* number of MPI communicators: ncomms = 2^comm_bits */
  nMPI_Comm *comm;   /* list of MPI communicators */
  int tag_ub;        /* min of upper bound of MPI tags over all in comm */
  int tag_bits;      /* number of bits in tag_ub+1 */
  nMPI_Datatype TELM0; //MPI_Datatype for tElm from 0-offsetof(tElm, dat)
  nMPI_Datatype TEPLOC; //MPI_Datatype for tEploc
} tnMPIvars;


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
int nMPI_Type_contiguous(int count, nMPI_Datatype oldtype,
                         nMPI_Datatype *newtype);
int nMPI_Type_commit(nMPI_Datatype *datatype);
int nMPI_Type_free(nMPI_Datatype *datatype);
int nMPI_Send(const void *buf, int count, nMPI_Datatype datatype,
              int dest, int tag);
int nMPI_Recv(void *buf, int count, nMPI_Datatype datatype,
              int source, int tag);
int nMPI_Isend(const void *buf, int count, nMPI_Datatype datatype,
               int dest, int tag, nMPI_Comm comm, nMPI_Req *req);
int nMPI_Irecv(void *buf, int count, nMPI_Datatype datatype,
               int src, int tag, nMPI_Comm comm, nMPI_Req *req);
int nMPI_Waitall(int nreq, nMPI_Req *req, nMPI_Stat *stat);
int nMPI_Wait(nMPI_Req *req, nMPI_Stat *stat);
int nMPI_Allreduce(const void *sendbuf, void *recvbuf, int count,
                   nMPI_Datatype datatype, nMPI_Op op);
int nMPI_Reduce(const void *sendbuf, void *recvbuf, int count,
                nMPI_Datatype datatype, nMPI_Op op, int root);
int nMPI_Bcast(void *buffer, int count, nMPI_Datatype datatype, int root);
int nMPI_Ibcast(void *buffer, int count, nMPI_Datatype datatype,
                int root, nMPI_Req *request);
int nMPI_Test(nMPI_Req *request, int *flag, nMPI_Stat *stat);
int nMPI_Alloc_mem(size_t size, nMPI_Info info, void *baseptr);
int nMPI_Free_mem(void *base);
int nMPI_Win_allocate(size_t size, int disp_unit, nMPI_Info info,
                      nMPI_Comm comm, void *baseptr, nMPI_Win *win);
int nMPI_Win_create(void *base, size_t size, int disp_unit,
                    nMPI_Info info, nMPI_Comm comm, nMPI_Win *win);
int nMPI_Win_free(nMPI_Win *win);
int nMPI_Put(void *origin_addr, int origin_count,
             nMPI_Datatype origin_datatype,
             int target_rank, size_t target_disp, int target_count,
             nMPI_Datatype target_datatype, nMPI_Win win);
int nMPI_Get(void *origin_addr, int origin_count,
             nMPI_Datatype origin_datatype,
             int target_rank, size_t target_disp, int target_count,
             nMPI_Datatype target_datatype, nMPI_Win win);

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
void nMPI_Isend_Irecv_com(tCom *com, int rq, nMPI_Datatype datatype,
                          int rank_other, int s_tag, int r_tag,
                          nMPI_Comm s_comm, nMPI_Comm r_comm);
int nMPI_Isend_com(tCom *com, int rq, nMPI_Datatype datatype,
                   int dest, int tag, nMPI_Comm comm);
int nMPI_Irecv_com(tCom *com, int rq, nMPI_Datatype datatype,
                   int src, int tag, nMPI_Comm comm);
void nMPI_Isend_Irecv_double_com(tCom *com, int rq,
                                 int rank_other, int s_tag, int r_tag,
                                 nMPI_Comm s_comm, nMPI_Comm r_comm);
int nMPI_Isend_double_com(tCom *com, int rq, int dest, int tag, nMPI_Comm comm);
int nMPI_Irecv_double_com(tCom *com, int rq, int src, int tag, nMPI_Comm comm);



/**************************************************************************/
/* insert lists.h here because nMPI is the 2nd module,
   so all after will see the lists types and such */
/**************************************************************************/
#include "../main/lists.h"
