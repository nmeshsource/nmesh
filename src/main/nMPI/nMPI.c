/* nMPI.c */
/* Wolfgang Tichy, 1/2019 */

#include "nmesh.h"
#include "nMPI.h"

#define PR 1

/* print some compile info */
int nMPI_print_compile_info(tMesh *mesh)
{
#ifdef USEMPI
  printf("MPI is compiled in.\n");
#else
  printf("MPI is not compiled in.\n");
#endif
  return 0;
}

/* Wrappers for MPI_Init and MPI_Finalize */
int nMPI_Init(int *pargc, char ***pargv)
{
#ifdef USEMPI
  return MPI_Init(pargc, pargv);
#endif
  return 0;
}
int nMPI_Finalize(void)
{
#ifdef USEMPI
  return MPI_Finalize();
#endif
  return 0;
}

/* return MPI rank, if MPI is not compiled in return 0 */
int nMPI_rank(void)
{
  int rank=0;
#ifdef USEMPI
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
#endif
  return rank;
}

/* number of processes we are running with, or 1 if no MPI */
int nMPI_size(void)
{
  int size=1;
#ifdef USEMPI
  MPI_Comm_size(MPI_COMM_WORLD, &size);
#endif
  //size=2;
  return size;
}

/* barrier */
int nMPI_barrier(void)
{
  int ret=0;
#ifdef USEMPI
  ret = MPI_Barrier(MPI_COMM_WORLD);
#endif
  return ret;
}


/* exchange double buffers */
void nMPI_Isend_Irecv_double(double *sbuf, int ns, double *rbuf, int nr,
                             int rank_other, int s_tag, int r_tag,
                             nMPI_Req *s_req, nMPI_Req *r_req)
{
#ifdef USEMPI
  int errS, errR;

  PRF;printf(": %d to %d, ns=%d nr=%d s_tag=%d r_tag=%d\n",
             nMPI_rank(), rank_other, ns, nr, s_tag, r_tag);
  fflush(stdout);

  errS = MPI_Isend(sbuf, ns, MPI_DOUBLE, rank_other, s_tag,
                   MPI_COMM_WORLD, s_req);
  if(errS != MPI_SUCCESS) errorexit("MPI_Isend failed!\n");
  
  errR = MPI_Irecv(rbuf, nr, MPI_DOUBLE, rank_other, r_tag,
                   MPI_COMM_WORLD, r_req);
  if(errR != MPI_SUCCESS) errorexit("MPI_Irecv failed!\n");
#endif
}

/* check on requests */
int nMPI_Waitall(int nreq, nMPI_Req *req, nMPI_Stat *stat)
{
  int status = 0;
#ifdef USEMPI
  PRF;printf(": %d waiting for %d reqs to finish\n", nMPI_rank(), nreq);
  fflush(stdout);

  status = MPI_Waitall(nreq, req, stat);

  if(status == MPI_ERR_IN_STATUS)
    errorexiti("MPI_Waitall error after waiting for %d requests",nreq);
#endif
  return status;
}

void nMPI_Isend_double(double *buf, int blen, int dest, int tag, nMPI_Req *req)
{
#ifdef USEMPI
  PRF;printf(": %d to %d, blen=%d tag=%d\n", nMPI_rank(), dest, blen, tag);
  MPI_Isend(buf, blen, MPI_DOUBLE, dest, tag, MPI_COMM_WORLD, req);
#endif
}

void nMPI_Irecv_double(double *buf, int blen, int src, int tag, nMPI_Req *req)
{
#ifdef USEMPI
  PRF;printf(": %d from %d, blen=%d tag=%d\n", nMPI_rank(), src, blen, tag);
  MPI_Irecv(buf, blen, MPI_DOUBLE, src, tag, MPI_COMM_WORLD, req);
#endif
}


/**********************************************************************/
/* deal with MPI tCom struct */
/**********************************************************************/
/* create a com */
tCom *alloc_com(int entrysize, int free_buf)
{
  tCom *com;

  com = calloc(1, sizeof(tCom));
  if(!com) errorexit("out of memory for com");
  com->entrysize = entrysize;
  com->free_buf = free_buf;

  return com;
}

/* free com and all within it */
void free_com(tCom *com)
{
  int i;

  if(!com) return;

  free(com->send_rq);
  free(com->recv_rq);
  free(com->send_stat);
  free(com->recv_stat);
  if(com->free_buf)
    for(i=0; i<com->n_rq; i++)
    {
      free(com->send_buf[i]);
      free(com->recv_buf[i]);
    }
  free(com->send_buf);
  free(com->recv_buf);

  free(com);
}

/* alloc MPI requests */
void realloc_com_reqs(tCom *com, int n_rq_new)
{
  int n_rq;

  if(!com) return;
  n_rq = com->n_rq;

  /* free buffer contents */
  if( (n_rq_new < com->n_rq) && (com->free_buf) )
  {
    int i;
    for(i=n_rq_new; i<com->n_rq; i++)
    {
      free(com->send_buf[i]);
      free(com->recv_buf[i]);
    }
  }

  if(n_rq_new)
  {
    com->send_rq = realloc(com->send_rq, n_rq_new*sizeof(com->send_rq[0]));
    if(!com->send_rq) errorexit("out of memory for com->send_rq");
    com->recv_rq = realloc(com->recv_rq, n_rq_new*sizeof(com->recv_rq[0]));
    if(!com->recv_rq) errorexit("out of memory for com->recv_rq");

    com->send_stat = realloc(com->send_stat, n_rq_new*sizeof(com->send_stat[0]));
    if(!com->send_stat) errorexit("out of memory for com->send_stat");
    com->recv_stat = realloc(com->recv_stat, n_rq_new*sizeof(com->recv_stat[0]));
    if(!com->recv_stat) errorexit("out of memory for com->recv_stat");

    com->send_buflen = realloc(com->send_buflen, n_rq_new*sizeof(com->send_buflen[0]));
    if(!com->send_buflen) errorexit("out of memory for com->send_buflen");
    com->recv_buflen = realloc(com->recv_buflen, n_rq_new*sizeof(com->recv_buflen[0]));
    if(!com->recv_buflen) errorexit("out of memory for com->recv_buflen");

    /* realloc buffer lists */
    com->send_buf = realloc(com->send_buf, n_rq_new*com->entrysize);
    if(!com->send_buf) errorexit("out of memory for com->send_buf");
    com->recv_buf = realloc(com->recv_buf, n_rq_new*com->entrysize);
    if(!com->recv_buf) errorexit("out of memory for com->recv_buf");

    /* zero new stuff */
    if(n_rq_new > n_rq)
    {
      memset(com->send_rq + n_rq, 0, sizeof(com->send_rq[0])*(n_rq_new-n_rq));
      memset(com->send_stat + n_rq, 0, sizeof(com->send_stat[0])*(n_rq_new-n_rq));
      memset(com->send_buflen + n_rq, 0, sizeof(com->send_buflen[0])*(n_rq_new-n_rq));
      memset(com->send_buf + n_rq, 0, com->entrysize*(n_rq_new-n_rq));
      memset(com->recv_rq + n_rq, 0, sizeof(com->recv_rq[0])*(n_rq_new-n_rq));
      memset(com->recv_stat + n_rq, 0, sizeof(com->recv_stat[0])*(n_rq_new-n_rq));
      memset(com->recv_buflen + n_rq, 0, sizeof(com->recv_buflen[0])*(n_rq_new-n_rq));
      memset(com->recv_buf + n_rq, 0, com->entrysize*(n_rq_new-n_rq));
    }
  }
  else
  {
    free(com->send_rq);
    free(com->recv_rq);
    free(com->send_stat);
    free(com->recv_stat);
    free(com->send_buflen);
    free(com->recv_buflen);

    free(com->send_buf);
    free(com->recv_buf);
  }
  com->n_rq = n_rq_new;
}

/* point com buffers */
void put_buffers_in_com(tCom *com, int rq,
                        void *sbuf, int slen, void *rbuf, int rlen)
{
  com->send_buf[rq] = sbuf;
  com->send_buflen[rq] = slen;
  com->recv_buf[rq] = rbuf;
  com->recv_buflen[rq] = rlen;
}

/* get com buffer pointers */
void *get_com_send_buf(tCom *com, int rq)
{
  return com->send_buf[rq];
}
void *get_com_recv_buf(tCom *com, int rq)
{
  return com->recv_buf[rq];
}

/* wait for all send requests in com */
int nMPI_Waitall_com_send(tCom *com)
{
  return nMPI_Waitall(com->n_rq, com->send_rq, com->send_stat);
}

/* wait for all recv requests in com */
int nMPI_Waitall_com_recv(tCom *com)
{
  return nMPI_Waitall(com->n_rq, com->recv_rq, com->recv_stat);
}

/* wait for all requests in com */
int nMPI_Waitall_com(tCom *com)
{
  int status=0;
  PRF;printf("\n");
  status  = nMPI_Waitall_com_recv(com);
  status += nMPI_Waitall_com_send(com);
  return status;
}

/* do requst rq of com */
void nMPI_Isend_Irecv_double_com(tCom *com, int rq,
                                 int rank_other, int s_tag, int r_tag)
{
  nMPI_Isend_Irecv_double(com->send_buf[rq], com->send_buflen[rq],
                          com->recv_buf[rq], com->recv_buflen[rq],
                          rank_other, s_tag, r_tag,
                          &(com->send_rq[rq]), &(com->recv_rq[rq]));
}

void nMPI_Isend_double_com(tCom *com, int rq, int dest, int tag)
{
  nMPI_Isend_double(com->send_buf[rq], com->send_buflen[rq], dest, tag,
                    &(com->send_rq[rq]));
}

void nMPI_Irecv_double_com(tCom *com, int rq, int src, int tag)
{
  nMPI_Irecv_double(com->recv_buf[rq], com->recv_buflen[rq], src, tag,
                    &(com->recv_rq[rq]));
}
