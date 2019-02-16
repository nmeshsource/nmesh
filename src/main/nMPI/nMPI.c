/* nMPI.c */
/* Wolfgang Tichy, 1/2019 */

#include "nmesh.h"
#include "nMPI.h"

#define PR 1

/* my rank and size if MPI is not compiled in */
int noMPI_rank=0, noMPI_size=1;


/* print some compile info */
int nMPI_print_compile_info(tMesh *mesh)
{
#ifdef USEMPI
  printf("MPI is compiled in.\n");
#else
  printf("MPI is not compiled in. Posing as rank=%d and size=%d.\n",
         noMPI_rank, noMPI_size);
#endif
  return 0;
}

/* Wrappers for MPI_Init and MPI_Finalize */
int nMPI_Init(int *pargc, char ***pargv)
{
#ifdef USEMPI
  return MPI_Init(pargc, pargv);
#else
  /* for debugging we can start nmesh as:
     nmesh nam.par 2 5
     then it poses as rank=2 and size=5 */
  if(pargc[0]==4)
  {
    noMPI_rank = atoi(pargv[0][2]);
    noMPI_size = atoi(pargv[0][3]);
    if(noMPI_rank<0) noMPI_rank=0;
    if(noMPI_size<1) noMPI_size=1;
  }
  return 0;
#endif
}
int nMPI_Finalize(void)
{
  fclose(stderr);
  fclose(stdout);
#ifdef USEMPI
  return MPI_Finalize();
#endif
  return 0;
}

/* return MPI rank, if MPI is not compiled in return 0 */
int nMPI_rank(void)
{
  int rank = noMPI_rank;
#ifdef USEMPI
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
#endif
  return rank;
}

/* number of processes we are running with, or 1 if no MPI */
int nMPI_size(void)
{
  int size = noMPI_size;
#ifdef USEMPI
  MPI_Comm_size(MPI_COMM_WORLD, &size);
#endif
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

/* non-blocking send */
void nMPI_Isend_double(double *buf, int blen, int dest, int tag, nMPI_Req *req)
{
  PRF;printf(": %d to %d, blen=%d tag=%d\n", nMPI_rank(), dest, blen, tag);
#ifdef USEMPI
  MPI_Isend(buf, blen, MPI_DOUBLE, dest, tag, MPI_COMM_WORLD, req);
#endif
}

/* non-blocking recv */
void nMPI_Irecv_double(double *buf, int blen, int src, int tag, nMPI_Req *req)
{
  PRF;printf(": %d from %d, blen=%d tag=%d\n", nMPI_rank(), src, blen, tag);
#ifdef USEMPI
  MPI_Irecv(buf, blen, MPI_DOUBLE, src, tag, MPI_COMM_WORLD, req);
#endif
}

/* exchange double buffers */
void nMPI_Isend_Irecv_double(double *sbuf, int ns, double *rbuf, int nr,
                             int rank_other, int s_tag, int r_tag,
                             nMPI_Req *s_req, nMPI_Req *r_req)
{
#ifdef USEMPI
  int errS, errR;
#endif

  PRF;printf(": %d to %d, ns=%d nr=%d s_tag=%d r_tag=%d\n",
             nMPI_rank(), rank_other, ns, nr, s_tag, r_tag);
  fflush(stdout);
#ifdef USEMPI
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
  if(!nreq) return 0;
  PRF;printf(": %d waiting for %d requests to finish\n", nMPI_rank(), nreq);
#ifdef USEMPI
  fflush(stdout);

  status = MPI_Waitall(nreq, req, stat);

  if(status == MPI_ERR_IN_STATUS)
    errorexiti("MPI_Waitall error after waiting for %d requests",nreq);
#endif
  return status;
}

/* check on requests */
int nMPI_Wait(nMPI_Req *req, nMPI_Stat *stat)
{
  int status = 0;
  PRF;printf(": %d waiting for request to finish\n", nMPI_rank());
#ifdef USEMPI
  fflush(stdout);

  status = MPI_Wait(req, stat);

  if(status == MPI_ERR_IN_STATUS)
    errorexit("MPI_Wait error after waiting for request");
#endif
  return status;
}



/**********************************************************************/
/* deal with tCom struct for MPI */
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
  free(com->send_buflen);
  free(com->recv_buflen);
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
    int entrysize = com->entrysize; /* save to restore later after memset 0 */
    int free_buf  = com->free_buf;
    free(com->send_rq);
    free(com->recv_rq);
    free(com->send_stat);
    free(com->recv_stat);
    free(com->send_buflen);
    free(com->recv_buflen);

    free(com->send_buf);
    free(com->recv_buf);

    /* zero all in com */
    memset(com, 0, sizeof(tCom));
    com->entrysize = entrysize; /* restore */
    com->free_buf  = free_buf;
  }
  com->n_rq = n_rq_new;
}

void print_com(tCom *com)
{
  int n_rq = com->n_rq;
  printf("com: n_rq=%d send_i=%d recv_i=%d\n", n_rq, com->send_i, com->recv_i);
#ifndef USEMPI
  for(int i=0; i<n_rq; i++)
    printf("%d: send_rq=%d recv_rq=%d send_stat=%d recv_stat=%d\n",
           i, com->send_rq[i], com->recv_rq[i],
           com->send_stat[i], com->recv_stat[i]);
#endif
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

/* add new com entries and point its buffers to sbuf and rbuf */
int append_buffers_to_com(tCom *com, void *sbuf,int slen, void *rbuf,int rlen)
{
  int rq = com->n_rq;
  realloc_com_reqs(com, rq + 1);
  put_buffers_in_com(com, rq, sbuf,slen, rbuf,rlen);
  return rq;
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
/* set the counter in com */
void set_com_counters(tCom *com, int si, int ri)
{
  com->send_i = si;
  com->recv_i = ri;
}
/* increase the counters */
void inc_com_send_i(tCom *com)
{
  com->send_i += 1;
}
void inc_com_recv_i(tCom *com)
{
  com->recv_i += 1;
}

/* get com buffer pointers one after the other, i.e. increase i */
void *get_com_send_i_buf_inc_i(tCom *com)
{
  int i = com->send_i;
  if(i >= com->n_rq) errorexit("send_i>=n_rq");
  com->send_i = i+1;
  return get_com_send_buf(com, i);
}
void *get_com_recv_i_buf_inc_i(tCom *com)
{
  int i = com->recv_i;
  if(i >= com->n_rq) errorexit("recv_i>=n_rq");
  com->recv_i = i+1;
  //PRF;printf(": i=%d -> recv_i=%d\n", i, com->recv_i); fflush(stdout);
  return get_com_recv_buf(com, i);
}

/* get com buffer pointers at i */
void *get_com_send_i_buf(tCom *com)
{
  return com->send_buf[com->send_i];
}
void *get_com_recv_i_buf(tCom *com)
{
  return com->recv_buf[com->recv_i];
}

/* free buffers */
void free_com_send_i_buf(tCom *com)
{
  int i = com->send_i;
  free(com->send_buf[i]);
  com->send_buf[i] = NULL;
  com->send_buflen[i] = 0;
}
void free_com_recv_i_buf(tCom *com)
{
  int i = com->recv_i;
  //PRF;printf(": i=%d -> recv_i=%d com->recv_buf[i]=%p\n",
  //           i, com->recv_i, com->recv_buf[i]); fflush(stdout);
  free(com->recv_buf[i]);
  com->recv_buf[i] = NULL;
  com->recv_buflen[i] = 0;
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

/* wait for send request rq to finish */
int nMPI_Wait_com_send(tCom *com, int rq)
{
  return nMPI_Wait(&(com->send_rq[rq]), &(com->send_stat[rq]));
}
/* wait for recv request rq to finish */
int nMPI_Wait_com_recv(tCom *com, int rq)
{
  return nMPI_Wait(&(com->recv_rq[rq]), &(com->recv_stat[rq]));
}

/* do send and recv request rq of com */
void nMPI_Isend_Irecv_double_com(tCom *com, int rq,
                                 int rank_other, int s_tag, int r_tag)
{
  nMPI_Isend_Irecv_double(com->send_buf[rq], com->send_buflen[rq],
                          com->recv_buf[rq], com->recv_buflen[rq],
                          rank_other, s_tag, r_tag,
                          &(com->send_rq[rq]), &(com->recv_rq[rq]));
}
/* send only */
void nMPI_Isend_double_com(tCom *com, int rq, int dest, int tag)
{
  nMPI_Isend_double(com->send_buf[rq], com->send_buflen[rq], dest, tag,
                    &(com->send_rq[rq]));
}
/* recv only */
void nMPI_Irecv_double_com(tCom *com, int rq, int src, int tag)
{
  nMPI_Irecv_double(com->recv_buf[rq], com->recv_buflen[rq], src, tag,
                    &(com->recv_rq[rq]));
}
