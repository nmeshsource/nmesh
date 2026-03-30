/* nMPI.c */
/* Wolfgang Tichy, 1/2019 */

#include "nmesh.h"
#include "nMPI.h"

#define PR 0
#define PR0 if(PR){PRF;printf(": calling MPI function...\n");fflush(stdout);}
#define PR1 if(PR){PRF;printf(": done.\n");fflush(stdout);}


/* my rank and size if MPI is not compiled in, look at nMPI_Init for
   debugging help!!! */
int noMPI_rank=0, noMPI_size=1;

/* global vars of nMPI */
tnMPIvars nMPIvars[1];


/* init global struct */
int nMPIvars_init(tMesh *mesh)
{
  int comm_bits = Geti(Par("nMPI_communicator_bits"));
  int ncomms = 1<<comm_bits;
  char cname[nMPI_MAX_OBJECT_NAME];
  int tag_ub, tag_bits;
  int i;

  /* get mem for communicators */
  nMPIvars->comm = calloc(ncomms, sizeof(nMPIvars->comm[0]));
  nMPIvars->ncomms = ncomms;
  nMPIvars->comm_bits = comm_bits;

  /* make communicators with nMPI_Comm_dup */
  for(i=0; i<ncomms; i++)
  {
    nMPIvars->comm[i] = nMPI_COMM_NULL;
    nMPI_Comm_dup(main_comm, &(nMPIvars->comm[i]));
    snprintf(cname,nMPI_MAX_OBJECT_NAME, "c%d", i); //name for dupd comm
    nMPI_Comm_set_name(nMPIvars->comm[i], cname);
  }

  /* get max number of MPI tags */
  tag_ub = INT_MAX;
  for(i=0; i<ncomms; i++)
  {
    int flag;
    void *v;

    nMPI_Comm_get_attr(nMPIvars->comm[i], nMPI_TAG_UB, &v, &flag);
    if(flag)
    {
      int ub = *((int *)v); /* upper bound */
      /* store min of ub */
      if(ub < tag_ub)
        tag_ub = ub;
    }
  }
  nMPIvars->tag_ub = tag_ub;

  /* set tag_bits */
  for(tag_bits=0; (tag_ub)>>tag_bits; tag_bits++);
  nMPIvars->tag_bits = tag_bits;

  /* set some MPI_Datatypes */
  nMPI_Type_contiguous(sizeof(tElm0), nMPI_CHAR, &(nMPIvars->TELM0));
  nMPI_Type_commit(&(nMPIvars->TELM0));
  nMPI_Type_contiguous(sizeof(tEploc), nMPI_CHAR, &(nMPIvars->TEPLOC));
  nMPI_Type_commit(&(nMPIvars->TEPLOC));

  return 0;
}

/* free global struct */
int nMPIvars_finalize(tMesh *mesh)
{
  int i;

  /* delete communicators */
  for(i=0; i<nMPIvars->ncomms; i++)
    if(nMPIvars->comm[i] != nMPI_COMM_NULL)
      nMPI_Comm_free(&(nMPIvars->comm[i]));

  /* free mem for communicators */
  free(nMPIvars->comm);

  /* free MPI_Datatypes we have created */
  nMPI_Type_free(&(nMPIvars->TELM0));
  nMPI_Type_free(&(nMPIvars->TEPLOC));

  return 0;
}

/* return nMPIvars->comm[i] */
nMPI_Comm nMPIvars_get_comm(int i)
{
  if(i>=nMPIvars->ncomms || i<0) return nMPI_COMM_NULL;
  return nMPIvars->comm[i];
}
int nMPIvars_get_ncomms(void)
{
  return nMPIvars->ncomms;
}
int nMPIvars_get_tag_ub(void)
{
  return nMPIvars->tag_ub;
}

/* take a long tag number and split it into 2 pieces:
   1. index commi of one of our communicators in nMPIvars->comm
   2. smaller tag number */
int nMPI_long_tag_to_commi_tag(long long_tag, int *commi, int *tag)
{
  int tag_bits = nMPIvars->tag_bits;
  long ci = long_tag>>tag_bits;
  long ntag = long_tag - (ci<<tag_bits);

  *commi = ci;
  *tag   = ntag;

  if(ci < nMPIvars->ncomms) return 0;
  return -1;
}


/* print some compile info */
int nMPI_print_compile_info(tMesh *mesh)
{
  PRFs(":\n");
#ifdef USEMPI
  printf(" MPI is compiled in. MPI_VERSION = %d\n", MPI_VERSION);
#else
  printf(" MPI is not compiled in. Posing as rank=%d and size=%d.\n",
         noMPI_rank, noMPI_size);
#endif
  printf(" nMPI_rank() = %d\n", nMPI_rank());
  printf(" nMPI_size() = %d\n", nMPI_size());
  printf(" nMPIvars->ncomms = %d\n", nMPIvars_get_ncomms());
  printf(" nMPIvars->comm_bits = %d\n", nMPIvars->comm_bits);
  printf(" nMPIvars->tag_ub = %d\n", nMPIvars_get_tag_ub());
  printf(" nMPIvars->tag_bits = %d\n", nMPIvars->tag_bits);

  return 0;
}


/************************************************************************/
/* funcs to deal with MPI errors */
/************************************************************************/

/* WT util: print MPI error */
void nMPI_print_error(int errcode)
{
#ifdef USEMPI
  char errbuffer[MPI_MAX_ERROR_STRING];
  int errlen;
  MPI_Error_string(errcode, errbuffer, &errlen);
  printf("MPI-error %d: %s\n", errcode, errbuffer);
#endif
}

/* WT util: print MPI_Status struct */
void nMPI_print_Stat(nMPI_Stat stat)
{
#ifdef USEMPI
  printf("ERR,SRC,TAG=%d,%d,%d",
         stat.MPI_ERROR, stat.MPI_SOURCE, stat.MPI_TAG);
#endif
}

/* WT util: print errors in an array of MPI_Status except skip_code */
void nMPI_print_error_MPI_Stat_array_skip(int nreq, nMPI_Stat *stat,
                                          int skip_code)
{
#ifdef USEMPI
  int i;
  for(i=0; i<nreq; i++)
  {
    nMPI_Stat st = stat[i];
    if( (st.MPI_ERROR != skip_code) && (st.MPI_ERROR != MPI_SUCCESS) )
    {
      printf("  MPI-error %d in request %d with src %d and tag %d\n",
             st.MPI_ERROR, i, st.MPI_SOURCE, st.MPI_TAG);
      printf("    ");
      nMPI_print_error(st.MPI_ERROR);
    }
  }
#endif
}
/* WT util: print all errors in an array of MPI_Status */
void nMPI_print_error_MPI_Stat_array(int nreq, nMPI_Stat *stat)
{
#ifdef USEMPI
  nMPI_print_error_MPI_Stat_array_skip(nreq, stat, MPI_SUCCESS);
#endif
}

/* Exit if we have an MPI error. This works only if the MPI error handler
   is set to MPI_ERRORS_RETURN */
void nMPI_check_error(const char *file, int line, const char *func, int stat)
{
#ifdef USEMPI
  if(stat != MPI_SUCCESS)
  {
    printf("%s:%d: nMPI-Error %d in %s\n", file, line, stat, func);
    printf("  ");
    nMPI_print_error(stat);
    errorexit("Exiting nmesh due to nMPI-Error");
  }
#endif
}

/* WT util: print MPI_Comm */
void nMPI_print_Comm_name(nMPI_Comm comm)
{
#ifdef USEMPI
  char comm_name[MPI_MAX_OBJECT_NAME+1];
  int resultlen;
  MPI_Comm_get_name(comm, comm_name, &resultlen);
  printf("%s", comm_name);
#endif
}


/********************************************************************/
/* Wrappers for MPI functions */
/********************************************************************/

/* Wrappers for MPI_Init and MPI_Finalize */
int nMPI_Init(int *pargc, char ***pargv)
{
  int ret=0;
#ifdef USEMPI
  int size=1, rank=0;
#if defined(USEOMP) && !defined(PLAIN_MPI_INIT)
  int required = MPI_THREAD_FUNNELED; /* only masterthread makes MPI calls */
  int provided;
  ret = MPI_Init_thread(pargc, pargv, required, &provided);
  PR1;
  if(provided < required) /* exit if MPI cannot do threads */
    errorexit("MPI_Init_thread: provided < required. "
              "MPI library doesn't support USEOMP.");
#else
  ret = MPI_Init(pargc, pargv); /* no OpenMP, no threads */
  PR1;
#endif
  /* perf programs like HPCtoolkit want calls to MPI_Comm_rank */
  MPI_Comm_size(MPI_COMM_WORLD, &size);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
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
#endif
  return ret;
}
int nMPI_Finalize(void)
{
  int ret=0;
#ifdef USEMPI
  PR0;
  ret = MPI_Finalize();
  PR1;
#endif
  /* We should never close stderr & stdout. If we close anything at all it
     would be the files they have been redirected to... */
  //fclose(stderr);
  //fclose(stdout);
  return ret;
}

/* abort MPI in case of errorexit */
int nMPI_Abort(int errorcode)
{
  fflush(stderr);
  fflush(stdout);
#ifdef USEMPI
  PR0;
  return MPI_Abort(WORLD, errorcode);
  PR1;
#endif
  return errorcode;
}

/* return MPI rank, if MPI is not compiled in return 0 */
int nMPI_rank(void)
{
  int rank = noMPI_rank;
#ifdef USEMPI
  //PR0;
  MPI_Comm_rank(WORLD, &rank);
  //PR1;
#endif
  return rank;
}

/* number of processes we are running with, or 1 if no MPI */
int nMPI_size(void)
{
  int size = noMPI_size;
#ifdef USEMPI
  PR0;
  MPI_Comm_size(WORLD, &size);
  PR1;
#endif
  return size;
}

/* barrier */
int nMPI_barrier(void)
{
  int ret=0;
#ifdef USEMPI
  PR0;
  ret = MPI_Barrier(WORLD);
  PR1;
#endif
  return ret;
}

/* get MPI attribute */
int nMPI_Comm_get_attr(nMPI_Comm comm, int comm_keyval,
                       void *attribute_val, int *flag)
{
  int ret=0;
  *flag=0;
#ifdef USEMPI
  PR0;
  ret = MPI_Comm_get_attr(comm, comm_keyval, attribute_val, flag);
  PR1;
#endif
  return ret;
}

/* duplicate a communicator */
int nMPI_Comm_dup(nMPI_Comm comm, nMPI_Comm *newcomm)
{
  int ret=0;
#ifdef USEMPI
  PR0;
  ret = MPI_Comm_dup(comm, newcomm);
  PR1;
#endif
  return ret;
}

/* free a communicator */
int nMPI_Comm_free(nMPI_Comm *comm)
{
  int ret=0;
#ifdef USEMPI
  PR0;
  ret = MPI_Comm_free(comm);
  PR1;
#endif
  return ret;
}

/* set a communicator name */
int nMPI_Comm_set_name(nMPI_Comm comm, const char *comm_name)
{
  int ret=0;
#ifdef USEMPI
  PR0;
  ret = MPI_Comm_set_name(comm, comm_name);
  PR1;
#endif
  return ret;
}

/* get a communicator name */
int nMPI_Comm_get_name(nMPI_Comm comm, char *comm_name, int *resultlen)
{
  int ret=0;
#ifdef USEMPI
  PR0;
  ret = MPI_Comm_get_name(comm, comm_name, resultlen);
  PR1;
#endif
  return ret;
}

/* set error handler */
int nMPI_Comm_set_errhandler(nMPI_Comm comm, nMPI_Errhandler errhandler)
{
  int ret=0;
#ifdef USEMPI
  PR0;
  ret = MPI_Comm_set_errhandler(comm, errhandler);
  PR1;
#endif
  return ret;
}

/* Create a new contiguous datatype */
int nMPI_Type_contiguous(int count, nMPI_Datatype oldtype,
                         nMPI_Datatype *newtype)
{
  int ret=0;
#ifdef USEMPI
  PR0;
  ret = MPI_Type_contiguous(count, oldtype, newtype);
  PR1;
#endif
  return ret;
}
/* Make a new datatype available */
int nMPI_Type_commit(nMPI_Datatype *datatype)
{
  int ret=0;
#ifdef USEMPI
  PR0;
  ret = MPI_Type_commit(datatype);
  PR1;
#endif
  return ret;
}
/* Free a datatype made with nMPI_Type_commit */
int nMPI_Type_free(nMPI_Datatype *datatype)
{
  int ret=0;
#ifdef USEMPI
  PR0;
  ret = MPI_Type_free(datatype);
  PR1;
#endif
  return ret;
}

/* general blocking send */
int nMPI_Send(const void *buf, int count, nMPI_Datatype datatype,
              int dest, int tag)
{
  int stat = 0;
  if(PR)
  {
    PRF;printf(": %d to %d, count=%d tag=%d\n", nMPI_rank(), dest, count, tag);
  }
#ifdef USEMPI
  PR0;
  stat = MPI_Send(buf,count, datatype, dest, tag, WORLD);
  PR1;
#endif
  return stat;
}

/* general blocking recv */
int nMPI_Recv(void *buf, int count, nMPI_Datatype datatype,
              int source, int tag)
{
  int stat = 0;
  if(PR)
  {
    PRF;printf(": %d from %d, count=%d tag=%d\n",
               nMPI_rank(), source, count, tag);
  }
#ifdef USEMPI
  PR0;
  stat = MPI_Recv(buf,count, datatype, source, tag, WORLD, MPI_STATUS_IGNORE);
  PR1;
#endif
  return stat;
}

/* general non-blocking send */
int nMPI_Isend(const void *buf, int count, nMPI_Datatype datatype,
               int dest, int tag, nMPI_Comm comm, nMPI_Req *req)
{
  int stat = 0;
  if(PR)
  {
    PRF;printf(": %d to %d, count=%d tag=%d\n", nMPI_rank(), dest, count, tag);
  }
  if(!count) errorexit("MPI_Isend allows count=0, but we do not!");
#ifdef USEMPI
  PR0;
  stat = MPI_Isend(buf, count, datatype, dest, tag, comm, req);
  if(stat != MPI_SUCCESS) {PRF;printf(": failed: %d\n", stat);}
  PR1;
#endif
  return stat;
}

/* general non-blocking recv */
int nMPI_Irecv(void *buf, int count, nMPI_Datatype datatype,
               int src, int tag, nMPI_Comm comm, nMPI_Req *req)
{
  int stat = 0;
  if(PR)
  {
    PRF;printf(": %d from %d, count=%d tag=%d\n", nMPI_rank(), src, count, tag);
  }
  if(!count) errorexit("MPI_Irecv allows count=0, but we do not!");
#ifdef USEMPI
  PR0;
  stat = MPI_Irecv(buf, count, datatype, src, tag, comm, req);
  if(stat != MPI_SUCCESS) {PRF;printf(": failed: %d\n", stat);}
  PR1;
#endif
  return stat;
}

/* exchange general buffers */
int nMPI_Isend_Irecv(void *sbuf, int ns, void *rbuf, int nr,
                     nMPI_Datatype datatype,
                     int rank_other, int s_tag, int r_tag,
                     nMPI_Comm s_comm, nMPI_Comm r_comm,
                     nMPI_Req *s_req, nMPI_Req *r_req)
{
  int stat = 0;
#ifdef USEMPI
  int errS, errR;
#endif
  if(PR)
  {
    PRF;printf(": %d to %d, ns=%d nr=%d s_tag=%d r_tag=%d\n",
                nMPI_rank(), rank_other, ns, nr, s_tag, r_tag);
    fflush(stdout);
  }
  //for(int i=0; i<ns; i++) printf(" %g", sbuf[i]);
  //printf("\n");
  //fflush(stdout);
  if(!ns || !nr)
    errorexit("MPI_Isend, MPI_Irecv allow ns=0, nr=0, but we do not!");
#ifdef USEMPI
  PR0;
  errS = MPI_Isend(sbuf, ns, datatype, rank_other, s_tag, s_comm, s_req);
  if(errS != MPI_SUCCESS)
  {
    PRF;printf(": MPI_Isend failed: %d\n", errS);
    return errS;
  }
  
  errR = MPI_Irecv(rbuf, nr, datatype, rank_other, r_tag, r_comm, r_req);
  if(errR != MPI_SUCCESS)
  {
    PRF;printf(": MPI_Irecv failed: %d\n", errR);
    return errR;
  }
  PR1;
  stat = MPI_SUCCESS;
#endif
  return stat;
}

/* non-blocking send for double */
int nMPI_Isend_double(double *buf, int blen, int dest, int tag,
                      nMPI_Comm comm, nMPI_Req *req)
{
  return nMPI_Isend(buf,blen, nMPI_DOUBLE, dest, tag, comm, req);
}

/* non-blocking recv for double */
int nMPI_Irecv_double(double *buf, int blen, int src, int tag,
                      nMPI_Comm comm, nMPI_Req *req)
{
  return nMPI_Irecv(buf,blen, nMPI_DOUBLE, src, tag, comm, req);
}

/* exchange double buffers */
int nMPI_Isend_Irecv_double(double *sbuf, int ns, double *rbuf, int nr,
                            int rank_other, int s_tag, int r_tag,
                            nMPI_Comm s_comm, nMPI_Comm r_comm,
                            nMPI_Req *s_req, nMPI_Req *r_req)
{
  return nMPI_Isend_Irecv(sbuf,ns, rbuf,nr, nMPI_DOUBLE, rank_other,
                          s_tag, r_tag, s_comm, r_comm, s_req, r_req);
}


/* check on requests */
int nMPI_Waitall(int nreq, nMPI_Req *req, nMPI_Stat *stat)
{
  int status = 0;
  if(!nreq) return nMPI_SUCCESS;
  if(PR)
  {
    int rank = nMPI_rank();
    PRF;printf(": %d waiting for %d requests to finish\n", rank, nreq);
    fflush(stdout);
  }
#ifdef USEMPI
  status = MPI_Waitall(nreq, req, stat);
  if(PR)
  {
    PRF;printf(": done waiting, status=%d\n", status);
    fflush(stdout);
  }
  if(status == MPI_ERR_IN_STATUS)
  {
    PRF;printf(": error after waiting for %d requests\n", nreq);
    nMPI_print_error(status);
    nMPI_print_error_MPI_Stat_array_skip(nreq, stat, MPI_ERR_PENDING);
  }
#endif
  return status;
}

/* check on requests */
int nMPI_Wait(nMPI_Req *req, nMPI_Stat *stat)
{
  int status = 0;
  if(PR)
  {
    int rank = nMPI_rank();
    PRF;printf(": %d waiting for request to finish\n", rank);
    fflush(stdout);
  }
#ifdef USEMPI
  status = MPI_Wait(req, stat);
  if(PR)
  {
    PRF;printf(": done waiting, status=%d\n", status);
    fflush(stdout);
  }
  if(status == MPI_ERR_IN_STATUS)
    {PRF;printf(": error after waiting for request\n");}
#endif
  return status;
}

/* reduction */
int nMPI_Allreduce(const void *sendbuf, void *recvbuf, int count,
                   nMPI_Datatype datatype, nMPI_Op op)
{
  int status = 0;
#ifdef USEMPI
  PR0;
  status = MPI_Allreduce(sendbuf, recvbuf, count, datatype, op, WORLD);
  PR1;
#endif
  return status;
}
int nMPI_Reduce(const void *sendbuf, void *recvbuf, int count,
                nMPI_Datatype datatype, nMPI_Op op, int root)
{
  int status = 0;
#ifdef USEMPI
  PR0;
  status = MPI_Reduce(sendbuf, recvbuf, count, datatype, op, root, WORLD);
  PR1;
#endif
  return status;
}

/* non-blocking reduction */
int nMPI_Iallreduce(const void *sendbuf, void *recvbuf, int count,
                    nMPI_Datatype datatype, nMPI_Op op, nMPI_Req *request)
{
  int status = 0;
#ifdef USEMPI
  PR0;
#if MPI_VERSION >= 3
  status = MPI_Iallreduce(sendbuf, recvbuf, count, datatype, op, WORLD,
                          request);
#else
  status = MPI_Allreduce(sendbuf, recvbuf, count, datatype, op, WORLD);
  *request = MPI_REQUEST_NULL;
#endif
  PR1;
#endif
  return status;
}


/* All send stuff in sendbuf to root, root recvs from all into recvbuf.
   recvbuf is ignored for all nonroot processes. */
int nMPI_Gather(const void *sendbuf, int sendcount, nMPI_Datatype sendtype,
                void *recvbuf, int recvcount, nMPI_Datatype recvtype, int root)
{
  int status = 0;
#ifdef USEMPI
  PR0;
  status = MPI_Gather(sendbuf, sendcount, sendtype,
                      recvbuf, recvcount, recvtype, root, WORLD);
  PR1;
#endif
  return status;
}

/* like nMPI_Gather, but returns immediately */
int nMPI_Igather(const void *sendbuf, int sendcount, nMPI_Datatype sendtype,
                 void *recvbuf, int recvcount, nMPI_Datatype recvtype,
                 int root, nMPI_Req *request)
{
  int status = 0;
#ifdef USEMPI
  PR0;
  status = MPI_Igather(sendbuf, sendcount, sendtype,
                       recvbuf, recvcount, recvtype, root, WORLD, request);
  PR1;
#endif
  return status;
}


/* blocking broadcast from rank root to all others in MPI_Comm WORLD */
int nMPI_Bcast(void *buffer, int count, nMPI_Datatype datatype, int root)
{
  int status = 0;
  if(PR)
  {
    PRF;printf(": count=%d root=%d\n", count, root);
  }
#ifdef USEMPI
  PR0;
  status = MPI_Bcast(buffer, count, datatype, root, WORLD);
  PR1;
#endif
  return status;
}

/* non-blocking broadcast from rank root to all others in MPI_Comm WORLD */
int nMPI_Ibcast(void *buffer, int count, nMPI_Datatype datatype,
                int root, nMPI_Req *request)
{
  int status = 0;
#ifdef USEMPI
  PR0;
#if MPI_VERSION >= 3
  status = MPI_Ibcast(buffer, count, datatype, root, WORLD, request);
#else
  status = MPI_Bcast(buffer, count, datatype, root, WORLD);
  *request = MPI_REQUEST_NULL;
#endif
  PR1;
#endif
  return status;
}


/* Tests for the completion of all previously initiated requests */
int nMPI_Testall(int nreq, nMPI_Req *req, int *flag, nMPI_Stat *stat)
{
  int status = 0;
  if(!nreq) return nMPI_SUCCESS;
  if(PR)
  {
    PRF;printf(": %d testing for %d requests to finish\n", nMPI_rank(), nreq);
    fflush(stdout);
  }
#ifdef USEMPI
  PR0;
  status = MPI_Testall(nreq, req, flag, stat);
  PR1;
  if(status == MPI_ERR_IN_STATUS)
  {
    PRF;printf(": error after testing %d requests\n",nreq);
    nMPI_print_error(status);
    nMPI_print_error_MPI_Stat_array_skip(nreq, stat, MPI_ERR_PENDING);
  }
#else
  *flag = 1;
#endif
  return status;
}

/* Test for the completion of a specific send or receive. */
int nMPI_Test(nMPI_Req *request, int *flag, nMPI_Stat *stat)
{
  int status = 0;
#ifdef USEMPI
  PR0;
  status = MPI_Test(request, flag, stat);
  PR1;
#else
  *flag = 1;
#endif
  return status;
}

/* alloc size bytes of memory */
int nMPI_Alloc_mem(size_t size, nMPI_Info info, void *baseptr)
{
  int status = 0;
#ifdef USEMPI
  PR0;
  status = MPI_Alloc_mem(size, info, baseptr);
  PR1;
#else
  // *baseptr = calloc(1, size);
  *( (char **) baseptr ) = calloc(1, size);
#endif
  return status;
}

/* free memory allocated with nMPI_Alloc_mem */
int nMPI_Free_mem(void *base)
{
  int status = 0;
#ifdef USEMPI
  PR0;
  status = MPI_Free_mem(base);
  PR1;
#else
  free(base);
#endif
  return status;
}

/* alloc and create window of size bytes */
int nMPI_Win_allocate(size_t size, int disp_unit, nMPI_Info info,
                      nMPI_Comm comm, void *baseptr, nMPI_Win *win)
{
  int status = 0;
#ifdef USEMPI
  PR0;
#if MPI_VERSION >= 2
  status = MPI_Win_allocate(size, disp_unit, info, comm, baseptr, win);
#else
  status = MPI_Alloc_mem(size, info, baseptr);
  status = MPI_Win_create(*( (char **) baseptr ), size, disp_unit,
                          info, comm, win);
#endif
  PR1;
#else
  // *baseptr = calloc(1, size);
  *( (char **) baseptr ) = calloc(1, size);
#endif
  return status;
}

/* create MPI window for RMA */
int nMPI_Win_create(void *base, size_t size, int disp_unit,
                    nMPI_Info info, nMPI_Comm comm, nMPI_Win *win)
{
  int status = 0;
#ifdef USEMPI
  PR0;
  status = MPI_Win_create(base, size,disp_unit, info, comm, win);
  PR1;
#endif
  return status;
}

/* remove a previously created window */
int nMPI_Win_free(nMPI_Win *win)
{
  int status = 0;
#ifdef USEMPI
  PR0;
  status = MPI_Win_free(win);
  PR1;
#endif
  return status;
}

/* Copies data from the origin memory to the target */
int nMPI_Put(void *origin_addr, int origin_count,
             nMPI_Datatype origin_datatype,
             int target_rank, size_t target_disp, int target_count,
             nMPI_Datatype target_datatype, nMPI_Win win)
{
  int status = 0;
#ifdef USEMPI
  PR0;
  status = MPI_Put(origin_addr, origin_count, origin_datatype, target_rank,
                   target_disp, target_count, target_datatype, win);
  PR1;
#endif
  return status;
}

/* Get data from a memory window on a remote process */
int nMPI_Get(void *origin_addr, int origin_count,
             nMPI_Datatype origin_datatype,
             int target_rank, size_t target_disp, int target_count,
             nMPI_Datatype target_datatype, nMPI_Win win)
{
  int status = 0;
#ifdef USEMPI
  PR0;
  status = MPI_Get(origin_addr, origin_count, origin_datatype, target_rank,
                   target_disp, target_count, target_datatype, win);
  PR1;
#endif
  return status;
}


/**********************************************************************/
/* deal with tCom struct for MPI */
/**********************************************************************/
/* create a com:
   if free_buf=1: send_buf[i] recv_buf[i] will be freed for all i
                  in realloc_com_reqs or free_com
   Note: com->entrysize=entrysize is not used by anything so far */
tCom *alloc_com(int entrysize, int free_buf)
{
  tCom *com;

  if(PR) PRFs(":\n");

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

  if(PR) PRFs(":\n");

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

  if(PR) { PRF;printf(": n_rq_new=%d\n", n_rq_new); }

  if(!com) return;
  n_rq = com->n_rq;

  /* free buffer contents */
  if( (n_rq_new < n_rq) && (com->free_buf) )
  {
    int i;
    for(i=n_rq_new; i<n_rq; i++)
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
    com->send_buf = realloc(com->send_buf, n_rq_new*sizeof(void *));
    if(!com->send_buf) errorexit("out of memory for com->send_buf");
    com->recv_buf = realloc(com->recv_buf, n_rq_new*sizeof(void *));
    if(!com->recv_buf) errorexit("out of memory for com->recv_buf");

    /* zero new stuff */
    if(n_rq_new > n_rq)
    {
      memset(com->send_rq + n_rq, 0, sizeof(com->send_rq[0])*(n_rq_new-n_rq));
      memset(com->send_stat + n_rq, 0, sizeof(com->send_stat[0])*(n_rq_new-n_rq));
      memset(com->send_buflen + n_rq, 0, sizeof(com->send_buflen[0])*(n_rq_new-n_rq));
      memset(com->send_buf + n_rq, 0, sizeof(void *)*(n_rq_new-n_rq));
      memset(com->recv_rq + n_rq, 0, sizeof(com->recv_rq[0])*(n_rq_new-n_rq));
      memset(com->recv_stat + n_rq, 0, sizeof(com->recv_stat[0])*(n_rq_new-n_rq));
      memset(com->recv_buflen + n_rq, 0, sizeof(com->recv_buflen[0])*(n_rq_new-n_rq));
      memset(com->recv_buf + n_rq, 0, sizeof(void *)*(n_rq_new-n_rq));
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

/* print one entry in com only */
void print_com_entry(tCom *com, int i)
{
#ifdef USEMPI
  printf("%d: send: buf=%d,%p stat.",
         i, com->send_buflen[i], (void *) com->send_buf[i]);
  nMPI_print_Stat(com->send_stat[i]);
  printf(" rq=%lx", (ulong) com->send_rq[i]);
  printf("  recv: buf=%d,%p stat.",
         com->recv_buflen[i], (void *) com->recv_buf[i]);
  nMPI_print_Stat(com->recv_stat[i]);
  printf(" rq=%lx\n", (ulong) com->recv_rq[i]);
#else
  printf("%d: send: buf=%d,%p stat=%d rq=%d",
         i, com->send_buflen[i], (void *) com->send_buf[i],
         com->send_stat[i], com->send_rq[i]);
  printf("  recv: buf=%d,%p stat=%d rq=%d\n",
         com->recv_buflen[i], (void *) com->recv_buf[i],
         com->recv_stat[i], com->recv_rq[i]);
#endif
}
/* print com header and entry i */
void print_com_at(tCom *com, int i)
{
  int n_rq = com->n_rq;

  if(PR) PRFs(":\n");

  printf("com%p: n_rq=%d send_i=%d recv_i=%d\n",
         (void *) com, n_rq, com->send_i, com->recv_i);
  print_com_entry(com, i);
}
/* print all in com */
void print_com(tCom *com)
{
  int n_rq = com->n_rq;

  if(PR) PRFs(":\n");

  printf("com%p: n_rq=%d send_i=%d recv_i=%d\n",
         (void *) com, n_rq, com->send_i, com->recv_i);

  for(int i=0; i<n_rq; i++) print_com_entry(com, i);
}

/* set free_buf flag in com */
void set_free_buf_in_com(tCom *com, int free_buf)
{
  if(PR) { PRF;printf(": free_buf=%d\n", free_buf); }

  com->free_buf  = free_buf;
}

/* point com buffers */
void put_buffers_in_com(tCom *com, int rq,
                        void *sbuf, int slen, void *rbuf, int rlen)
{
  if(PR) PRFs(":\n");

  com->send_buf[rq] = sbuf;
  com->send_buflen[rq] = slen;
  com->recv_buf[rq] = rbuf;
  com->recv_buflen[rq] = rlen;
}

/* add new com entries and point its buffers to sbuf and rbuf */
int append_buffers_to_com(tCom *com, void *sbuf,int slen, void *rbuf,int rlen)
{
  int rq = com->n_rq;

  if(PR) PRFs(":\n");

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

/* set com buffer pointers */
void set_com_send_buf(tCom *com, int rq, void *send_buf)
{
  com->send_buf[rq] = send_buf;
}
void set_com_recv_buf(tCom *com, int rq, void *recv_buf)
{
  com->recv_buf[rq] = recv_buf;
}

/* set the counter in com */
void set_com_counters(tCom *com, int si, int ri)
{
  com->send_i = si;
  com->recv_i = ri;
}
/* increase the counters, or wrap back to 0 */
void inc_com_send_i(tCom *com)
{
  int i = com->send_i + 1;
  if(i < com->n_rq)
    com->send_i = i;
  else
    com->send_i = 0;
}
void inc_com_recv_i(tCom *com)
{
  int i = com->recv_i + 1;
  if(i < com->n_rq)
    com->recv_i = i;
  else
    com->recv_i = 0;
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
  if(com->send_buf) return com->send_buf[com->send_i];
  else              return NULL;
}
void *get_com_recv_i_buf(tCom *com)
{
  if(com->recv_buf) return com->recv_buf[com->recv_i];
  else              return NULL;
}

/* free buffers */
void free_com_send_i_buf(tCom *com)
{
  int i = com->send_i;
  if(com->send_buf)
  {
    free(com->send_buf[i]);
    com->send_buf[i] = NULL;
    com->send_buflen[i] = 0;
  }
}
void free_com_recv_i_buf(tCom *com)
{
  int i = com->recv_i;
  //PRF;printf(": i=%d -> recv_i=%d com->recv_buf[i]=%p\n",
  //           i, com->recv_i, (void *) com->recv_buf[i]); fflush(stdout);
  if(com->recv_buf)
  {
    free(com->recv_buf[i]);
    com->recv_buf[i] = NULL;
    com->recv_buflen[i] = 0;
  }
}

/* wait for all send requests in com */
int nMPI_Waitall_com_send(tCom *com)
{
  int stat = nMPI_Waitall(com->n_rq, com->send_rq, com->send_stat);
#ifdef USEMPI
  if(stat != MPI_SUCCESS)
  {
    PRF;printf(": nMPI_Waitall failed: %d\n", stat);
    nMPI_print_error(stat);
  }
#endif
  return stat;
}

/* wait for all recv requests in com */
int nMPI_Waitall_com_recv(tCom *com)
{
  int stat =  nMPI_Waitall(com->n_rq, com->recv_rq, com->recv_stat);
#ifdef USEMPI
  if(stat != MPI_SUCCESS)
  {
    PRF;printf(": nMPI_Waitall failed: %d\n", stat);
    nMPI_print_error(stat);
  }
#endif
  return stat;
}

/* wait for all requests in com */
int nMPI_Waitall_com(tCom *com)
{
  int stat1, stat2;
  if(PR)
  {
    PRF;printf("\n");
  }
  stat1 = nMPI_Waitall_com_send(com);
  if(stat1 != nMPI_SUCCESS) return stat1;
  stat2 = nMPI_Waitall_com_recv(com);
  return stat2;
}

/* wait for send request rq to finish */
int nMPI_Wait_com_send(tCom *com, int rq)
{
  int stat;
  if(PR)
  {
    PRFs(": ");print_com_at(com, rq);
    PRF;printf(": request number rq=%d\n", rq);
  }
  stat = nMPI_Wait(&(com->send_rq[rq]), &(com->send_stat[rq]));
#ifdef USEMPI
  if(stat != MPI_SUCCESS)
  {
    PRFs(": ");print_com_at(com, rq);
    PRF;printf(": request number rq=%d\n", rq);
    PRF;printf(": nMPI_Wait failed: %d\n", stat);
  }
#endif
  return stat;
}
/* wait for recv request rq to finish */
int nMPI_Wait_com_recv(tCom *com, int rq)
{
  int stat;
  if(PR)
  {
    PRFs(": ");print_com_at(com, rq);
    PRF;printf(": request number rq=%d\n", rq);
  }
  stat = nMPI_Wait(&(com->recv_rq[rq]), &(com->recv_stat[rq]));
#ifdef USEMPI
  if(stat != MPI_SUCCESS)
  {
    PRFs(": ");print_com_at(com, rq);
    PRF;printf(": request number rq=%d\n", rq);
    PRF;printf(": nMPI_Wait failed: %d\n", stat);
  }
#endif
  return stat;
}

/* test if send request rq has finished */
int nMPI_Test_com_send(tCom *com, int rq, int *flag)
{
  int stat;
  if(PR)
  {
    PRFs(": ");print_com_at(com, rq);
    PRF;printf(": request number rq=%d\n", rq);
  }
  stat = nMPI_Test(&(com->send_rq[rq]), flag, &(com->send_stat[rq]));
#ifdef USEMPI
  if(stat != MPI_SUCCESS) {PRF;printf(": nMPI_Test failed: %d\n", stat);}
#endif
  return stat;
}
/* test if recv request rq has finished */
int nMPI_Test_com_recv(tCom *com, int rq, int *flag)
{
  int stat;
  if(PR)
  {
    PRFs(": ");print_com_at(com, rq);
    PRF;printf(": request number rq=%d\n", rq);
  }
  stat = nMPI_Test(&(com->recv_rq[rq]), flag, &(com->recv_stat[rq]));
#ifdef USEMPI
  if(stat != MPI_SUCCESS) {PRF;printf(": nMPI_Test failed: %d\n", stat);}
#endif
  return stat;
}

/* wait for all send requests in com */
int nMPI_Testall_com_send(tCom *com, int *flag)
{
  int stat = nMPI_Testall(com->n_rq, com->send_rq, flag, com->send_stat);
#ifdef USEMPI
  if(stat != MPI_SUCCESS) {PRF;printf(": nMPI_Testall failed: %d\n", stat);}
#endif
  return stat;
}
/* wait for all recv requests in com */
int nMPI_Testall_com_recv(tCom *com, int *flag)
{
  int stat = nMPI_Testall(com->n_rq, com->recv_rq, flag, com->recv_stat);
#ifdef USEMPI
  if(stat != MPI_SUCCESS) {PRF;printf(": nMPI_Testall failed: %d\n", stat);}
#endif
  return stat;
}
/* wait for all requests in com */
int nMPI_Testall_com(tCom *com, int *flag)
{
  int rflag, sflag;
  int stat1, stat2;
  if(PR)
  {
    PRF;printf("\n");
  }
  stat1 = nMPI_Testall_com_send(com, &sflag);
  if(stat1 != nMPI_SUCCESS) return stat1;
  stat2 = nMPI_Testall_com_recv(com, &rflag);
  *flag = rflag && sflag;
  return stat2;
}

/* do send and recv request rq of com */
int nMPI_Isend_Irecv_com(tCom *com, int rq, nMPI_Datatype datatype,
                         int rank_other, int s_tag, int r_tag,
                         nMPI_Comm s_comm, nMPI_Comm r_comm)
{
  if(PR)
  {
    PRFs(": ");print_com_at(com, rq);
    PRF;
    printf(": rq=%d rank_other=%d s_tag=%d r_tag=%d",
           rq, rank_other, s_tag, r_tag);
    printf(" s_comm=");nMPI_print_Comm_name(s_comm);
    printf(" r_comm=");nMPI_print_Comm_name(r_comm);
    printf("\n");
  }
  return nMPI_Isend_Irecv(com->send_buf[rq], com->send_buflen[rq],
                          com->recv_buf[rq], com->recv_buflen[rq],
                          datatype, rank_other, s_tag, r_tag, s_comm, r_comm,
                          &(com->send_rq[rq]), &(com->recv_rq[rq]));
}
/* send only */
int nMPI_Isend_com(tCom *com, int rq, nMPI_Datatype datatype,
                   int dest, int tag, nMPI_Comm comm)
{
  int stat = nMPI_Isend(com->send_buf[rq], com->send_buflen[rq], datatype,
                        dest, tag, comm, &(com->send_rq[rq]));
  return stat;
}
/* recv only */
int nMPI_Irecv_com(tCom *com, int rq, nMPI_Datatype datatype,
                   int src, int tag, nMPI_Comm comm)
{
  int stat = nMPI_Irecv(com->recv_buf[rq], com->recv_buflen[rq], datatype,
                        src, tag, comm, &(com->recv_rq[rq]));
  return stat;
}


/* do send and recv request rq of com for double */
int nMPI_Isend_Irecv_double_com(tCom *com, int rq,
                                int rank_other, int s_tag, int r_tag,
                                nMPI_Comm s_comm, nMPI_Comm r_comm)
{
  return nMPI_Isend_Irecv_com(com,rq, nMPI_DOUBLE,
                              rank_other, s_tag, r_tag, s_comm, r_comm);
}
/* send only */
int nMPI_Isend_double_com(tCom *com, int rq, int dest, int tag, nMPI_Comm comm)
{
  return nMPI_Isend_com(com,rq, nMPI_DOUBLE, dest, tag, comm);
}
/* recv only */
int nMPI_Irecv_double_com(tCom *com, int rq, int src, int tag, nMPI_Comm comm)
{
  return nMPI_Irecv_com(com,rq, nMPI_DOUBLE, src, tag, comm);
}
