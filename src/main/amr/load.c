/* load.c */
/* Wolfgang Tichy, 2/2019 */

#include "nmesh.h"
#include "amr.h"

#define PR 0


/**********************************************************************/
/* functions to move nodes between procs */
/**********************************************************************/

/* compute desired rank */
int desiredrank(long nid, long nnodes, int size)
{
  double N = nnodes;
  double s = size;
  double nperproc = N/s;
  double desrank = nid/nperproc;
  return desrank;
}

/* simplistic load balancing */
void simple_load_balance(tMesh *mesh)
{
  long nnodes = mesh->nln;
  long nid;
  tNlist *elem;
  tNode *node;
  int size = nMPI_size();
  int desrank;
/*
fornodelist(mesh->lns, elem)
{
node = elem->node;
tDat *dat = node->dat;
if(dat) dat->info->load_TimeIn_s = 1./(node->nid+1);
}
Yo(10);
printmesh(mesh);
load_balance_nodeload(mesh);
Yo(20);
printmesh(mesh);
return;
*/
  tCom *scom = alloc_com(sizeof(double), 1);
  tCom *rcom = alloc_com(sizeof(double), 1);

  PRF;printf(": nnodes=%ld\n", nnodes);

  /* free surfaces & indc since they will change now anyway */
  evolve_free_communication_structs(mesh);

  /* fill MPI send and recv buffers */
  fornodelist(mesh->lns, elem)
  {
    node = elem->node;
    nid = node->nid;
    desrank = desiredrank(nid, nnodes, size);
    if(node->datrank != desrank)
      move_node_to_rank(node, desrank, scom, rcom, 1);
  }
  nMPI_Waitall_com_send(scom);
  free_com(scom);  /* free scom with all its buffers */
  nMPI_Waitall_com_recv(rcom);

  /* get var data out of recv buffer */
  set_com_counters(rcom, 0,0);
  fornodelist(mesh->lns, elem)
  {
    node = elem->node;
    nid = node->nid;
    desrank = desiredrank(nid, nnodes, size);
    if(node->datrank != desrank)
      move_node_to_rank(node, desrank, scom, rcom, 0);
  }

  free_com(rcom);

  update_mesh_myln_node_nid(mesh);
  PRF;printf(": --> %d on this proc\n", total_nnodes_in_myln(mesh->myln));

  /* now that nodes are elsewhere re-init surfaces & indc */
  evolve_init_communication_structs(mesh);
}

/* return: number of variables and number of doubles inside dat */
/* NOTE: this does not take the extra space in some vars into account */
int nvars_ndoubles_in_dat(tDat *dat, int *ndoubles)
{
  tMesh *mesh;
  int nvars, vi, sizeofinfo;

  if(!dat)
  {
    *ndoubles = 0;
    return 0;
  }
  mesh = dat->node->pat->mesh;

  /* find amount of data in vars */
  for(nvars=0, *ndoubles=0, vi=0; vi<dat->nv; vi++)
    if(dat->v[vi])
    {
      *ndoubles += (dat->v[vi]->N) * (MeshVarType(mesh, vi)!=AUXVAR);
      if(PR) { PRF;printf(": vi=%d ndoubles=%d\n", vi, *ndoubles); }
      nvars++;
    }

  /* add amount in dat->info */
  sizeofinfo = sizeof(dat->info);
  *ndoubles += (sizeofinfo + sizeof(double)-1)/sizeof(double);

  return nvars;
}

/* pack all (non-aux) vars of dat into a buffer which contains:
  |nvars||varind1|npoints1|<--data1-->||varind2|npoints2|<--data2-->||...
  here nvars is the number of variables that we send.
  NOTE: this does not take the extra space in some vars into account.
  We also append dat->info to the end of this buffer:
  |sizeofinfo||info|
  The buffer has to be freed by caller later */
double *buffer_with_all_needed_dat_vars(tDat *dat, int *buflen)
{
  tMesh *mesh = dat->node->pat->mesh;
  int len;
  double *buf;
  int vi, datlen, nvars, bi, N;
  int sizeofinfo;

  /* find amount of data */
  nvars = nvars_ndoubles_in_dat(dat, &datlen);
  if(PR) { PRF;printf(": nvars=%d datlen=%d\n", nvars, datlen); }

  /* alloc buffer */
  len = 1 + nvars*2; /* for nvars, varinds, npoints */
  len += datlen;     /* for double data */
  len += 1;          /* for sizeofinfo */
  buf = calloc(len, sizeof(double));

  /* fill buffer with data in variables */
  buf[0] = nvars;
  for(bi=1, vi=0; vi<dat->nv; vi++)
  {
    /* add to buffer if eneabled and not auxiliary var */
    if(dat->v[vi])
    {
      N = (dat->v[vi]->N) * (MeshVarType(mesh, vi)!=AUXVAR);
      buf[bi++] = vi;
      buf[bi++] = N;
      memcpy(buf+bi, dat->v[vi]->d, N * sizeof(double));
      bi += N;
      if(PR) { PRF;printf(": vi=%d bi=%d\n", vi, bi); }
    }
  }

  /* now append the stuff in dat->info */
  sizeofinfo = sizeof(dat->info);
  buf[bi++] = sizeofinfo;
  memcpy(buf+bi, dat->info, sizeofinfo);
  bi += (sizeofinfo + sizeof(double)-1)/sizeof(double);

  /* return pointer to this buffer, and its length */
  *buflen = bi;
  return buf;
}
/* put buffer back into dat and enable all vars needed */
int write_buffer_into_dat_vars(tDat *dat, double *buf)
{
  tNode *node = dat->node;
  int i, vi, nvars, bi, N;
  int sizeofinfo;

  if(!buf) return 0;

  /* write buffer into vars */
  nvars = buf[0];
  for(bi=1, i=0; i<nvars; i++)
  {
    vi = buf[bi++];
    N  = buf[bi++];
    enablevarcomp_innode(node, vi);
    memcpy(dat->v[vi]->d, buf+bi, N * sizeof(double));
    bi += N;
    if(PR) { PRF;printf(": vi=%d bi=%d\n", vi, bi); }
  }

  /* now put the end of the buffer in dat->info */
  sizeofinfo = buf[bi++];
  memcpy(dat->info, buf+bi, sizeofinfo);

  return bi;
}

/* move data on one node between 2 ranks: the buffer we send/recv contains
  |nvars||varind1|npoints1|<--data1-->||varind2|npoints2|<--data2-->||... */
/* NOTE: this does not take the extra space in some vars into account */
void move_node_to_rank(tNode *node, int desrank,
                       tCom *scom, tCom *rcom, int setbufs)
{
  tDat *dat = node->dat;
  int slen, rlen;
  double *sbuf, *rbuf;
  int rank = nMPI_rank();
  int other, rq;

  TIMER_START;

  if(setbufs) /* setup buffers and fill them */
  {
    if(PR) { PRF;printf(": nid%ld datrank%d rank%d desrank%d\n",
                        node->nid, node->datrank, rank, desrank); }
    if(rank == node->datrank)
    {
      /* alloc and fill buffer */
      sbuf = buffer_with_all_needed_dat_vars(dat, &slen);
      if(PR)
      {
        PRF;printf(": sbuf =");
        for(int i=0; i<min2(3, slen); i++) printf(" %g", sbuf[i]);
        printf("\n");
        if(node->dat) printf("nv=%d\n", node->dat->nv);
      }
      other = desrank;
      /* put buffers in com */
      rq = append_buffers_to_com(scom, sbuf,slen, NULL,0);
      //print_com(scom);

      /* send */
      nMPI_Isend_double_com(scom, rq, other, node->nid, WORLD);
    }
    if(rank == desrank)
    {
      tMesh *mesh = node->pat->mesh;

      other = node->datrank;
      /* receive only if datrank is in valid range */
      if(other >= 0)
      {
        /* alloc buffer */
        rlen = 1 + (mesh->nvdb) * (2 + node->np); /* size to hold all vars */
        rlen += 1 + sizeof(node->dat->info); // size to hold info and its len
        rbuf = calloc(rlen, sizeof(double));
        /* put buffers in com */
        rq = append_buffers_to_com(rcom, NULL,0, rbuf,rlen);
        //print_com(rcom);
        /* receive */
        nMPI_Irecv_double_com(rcom, rq, other, node->nid, WORLD);
      }

      /* allocate space already and init some stuff */
      if(node->dat) errorexit("destination node should not have dat yet");
      node->dat = alloc_dat(node);

      if(0)
      {
        PRF;printf(": nid%ld rank%d node->dat=%p\n",
                   node->nid, rank, (void *) node->dat);
      }
      if(PR) { PRF;printf(": calling coordinates_init_node\n"); }
      coordinates_init_node(node);
    }
    if(PR) fflush(stdout);
  }
  else /* retrieve data from buffers */
  {
    if(rank == desrank)
    {
      /* now unpack the buffers */
      rbuf = get_com_recv_i_buf(rcom);
      write_buffer_into_dat_vars(node->dat, rbuf);
      /* we can free the buffer here already */
      free_com_recv_i_buf(rcom);
      inc_com_recv_i(rcom);
    }
    else
    {
      free_dat(node->dat);
      node->dat = NULL;
    }
    node->datrank = desrank;
    if(PR) fflush(stdout);
  }

  TIMER_STOP;
}

/* move all nodes in list to rank */
void move_nodelist_to_rank(tNlist *list, int desrank)
{
  tNlist *elem, *list0;
  tNode *node = list->node;
  //tMesh *mesh = node->pat->mesh;
  tCom *scom = alloc_com(sizeof(double), 1);
  tCom *rcom = alloc_com(sizeof(double), 1);

  if(PR) { PRF;printf(": desrank=%d\n", desrank); }

  /* find element 0 in list */
  list0 = first_nodelist(list);

  /* fill MPI send and recv buffers */
  fornodelist(list0, elem)
  {
    node = elem->node;
    if(node->datrank != desrank)
      move_node_to_rank(node, desrank, scom, rcom, 1);
  }
  nMPI_Waitall_com_send(scom);
  free_com(scom);  /* free scom with all its buffers */
  nMPI_Waitall_com_recv(rcom);

  /* get var data out of recv buffer */
  set_com_counters(rcom, 0,0);
  fornodelist(list0, elem)
  {
    node = elem->node;
    if(node->datrank != desrank)
      move_node_to_rank(node, desrank, scom, rcom, 0);
  }

  free_com(rcom);
}


/************************************************************************/
/* functions to time what's happening on a node */
/************************************************************************/

/* reset load timer on this node */
void loadtimer_reset(tNode *node)
{
  tDat *dat = node->dat;
  if(dat)
  {
    dat->info->load_timer_running = 0;
    dat->info->load_TimeIn_s      = 0.;
  }
}

/* start load timer on this node */
void loadtimer_start(tNode *node)
{
  tDat *dat = node->dat;
  if(!dat) errorexit("I can only time nodes that I have");

  /* save current time if timer is not running yet */
  if( !(dat->info->load_timer_running) )
  {
    dat->info->load_timer_running = 1;
    getRealTime(dat->info->load_start);
  }
  else
  {
    errorexit("load timer has already been started");
  }
}

/* start load timer and add time spent to counter dat->info->load_TimeIn_s*/
void loadtimer_stop(tNode *node)
{
  tDat *dat = node->dat;
  if(!dat) errorexit("I can only time nodes that I have");

  /* get time difference if timer was running */
  if( (dat->info->load_timer_running) )
  {
    struct timespec tp[1];

    getRealTime(tp);
    dat->info->load_TimeIn_s += getTimeDiffIn_s(tp, dat->info->load_start);

    dat->info->load_timer_running = 0;
  }
  else
  {
    errorexit("load timer is not running");
  }
}


/************************************************************************/
/* functions for load balancing based on node timings in
   dat->info->load_TimeIn_s */
/************************************************************************/

/* Set array nodeload[] that contains the measured time each node used,
   Also return the sum over nodeload[] */
double load_set_nodeload_array(tMesh *mesh, double *nodeload)
{
  tNlist *elem;
  double loadsum = 0.;

  /* we assume that mesh->lns has the same order for all ranks */
  fornodelist(mesh->lns, elem)
  {
    tNode *node = elem->node;
    tDat *dat = node->dat;
    int datrank = node->datrank;
    long nid = node->nid;
    double load = 0.;

    /* we need to broadcast nodeload from my nodes to all other ranks */
    if(dat) load = node->dat->info->load_TimeIn_s;

    /* in case we forgot to measure the times, just set load to a very small
       uniform number: */
    if(load == 0.) load = 1e-50;

    nMPI_Bcast(&load,1, nMPI_DOUBLE, datrank);
    nodeload[nid] = load;
    loadsum += load;
  }

  PRFs(":\n");
  for(long l=0; l<mesh->nln; l++)
    printf("nodeload[%ld]=%g ", l, nodeload[l]);

  return loadsum;
}

/* sum over leaf nodes (starting at ln0) until loadsum reaches the value
   desired_load:
   Return: leaf after the one where we reached desired_load.
   In: ln0, desired_load[], nodeload.  Out: actual_load */
tNlist *inc_leaf_until_desired_loadsum(tNlist *ln0, double desired_loadsum,
                                       const double nodeload[],
                                       double *actual_loadsum)
{
  tNlist *elem;
  double sum = 0.;

  fornodelist(ln0, elem)
  {
    tNode *node = elem->node;
    long nid = node->nid;

    sum += nodeload[nid];

    if(sum >= desired_loadsum) break;
  }

  *actual_loadsum = sum;
  if(elem) return elem->next;
  else     return NULL;
}

/* Set rank_start[i] array.
   It contains the first leaf that rank i should have */
void load_set_desired_rank_start(tMesh *mesh, double desired_loadsum,
                                 const double nodeload[],
                                 tNlist **rank_start)
{
  int size = nMPI_size();
  tNlist *ln0, *ln1;
  double actual_loadsum;
  int rank;

  rank = 0;
  for(ln0 = mesh->lns; ln0; ln0 = ln1)
  {
    if(rank >= size) break;
    rank_start[rank++] = ln0;
    ln1 = inc_leaf_until_desired_loadsum(ln0, desired_loadsum,
                                         nodeload, &actual_loadsum);
  }
  PRFs(":\n");
  for(rank=0; rank<size; rank++)
    printf("rank_start[%d]=nid%ld ", rank, rank_start[rank]->node->nid);
}

/* compute desired rank */
int load_get_desiredrank(long nid, tNlist **rank_start, int size)
{
  int rank;

  for(rank = 0; rank<size-1; rank++)
  {
    tNlist *ln1 = rank_start[rank+1];
    tNode *node = ln1->node;
    long nid1 = node->nid;
    //PRF;printf(": rank=%d nid=%ld nid1=%ld\n", rank, nid, nid1);
    /* we assume that leaf node nids are assigned in ascending order */
    if(nid < nid1) break;
  }
  if(rank >= size)
    errorexit("could not find the rank that should have nid");

  return rank;
}

/* load balancing based on measured node loads */
void load_balance_nodeload(tMesh *mesh)
{
  long nnodes = mesh->nln;
  int size = nMPI_size();
  double totalload;
  double *nodeload    = dmalloc(nnodes);
  tNlist **rank_start = calloc(size, sizeof(rank_start[0]));
  int desrank;
  long nid;
  tNlist *elem;
  tNode *node;
  tCom *scom = alloc_com(sizeof(double), 1);
  tCom *rcom = alloc_com(sizeof(double), 1);

  PRF;printf(": nnodes=%ld ", nnodes);
  if(!nodeload || !rank_start)
  {
    free_com(rcom);
    free_com(scom);
    free(rank_start);
    free(nodeload);
    printf("  WARNING: quitting ");PRF;
    printf(" due to lack of memory!!!\n");
    /* do fallback? */
    //simple_load_balance(mesh);
    return;
  }

  /* get measured load for each node (in nodeload) and also totalload */
  totalload = load_set_nodeload_array(mesh, nodeload);
  printf("totalload=%g\n", totalload);

  /* set array rank_start with 1st desired leaf node for each rank */
  load_set_desired_rank_start(mesh, totalload/size, nodeload, rank_start);

  /* free surfaces & indc since they will change now anyway */
  evolve_free_communication_structs(mesh);

  /* fill MPI send and recv buffers */
  fornodelist(mesh->lns, elem)
  {
    node = elem->node;
    nid = node->nid;
    desrank = load_get_desiredrank(nid, rank_start, size);
    if(node->datrank != desrank)
      move_node_to_rank(node, desrank, scom, rcom, 1);
  }
  nMPI_Waitall_com_send(scom);
  free_com(scom);  /* free scom with all its buffers */
  nMPI_Waitall_com_recv(rcom);

  /* get var data out of recv buffer */
  set_com_counters(rcom, 0,0);
  fornodelist(mesh->lns, elem)
  {
    node = elem->node;
    nid = node->nid;
    desrank = load_get_desiredrank(nid, rank_start, size);
    if(node->datrank != desrank)
      move_node_to_rank(node, desrank, scom, rcom, 0);
  }

  free_com(rcom);

  update_mesh_myln_node_nid(mesh);
  PRF;printf(": --> %d on this proc\n", total_nnodes_in_myln(mesh->myln));

  /* now that nodes are elsewhere re-init surfaces & indc */
  evolve_init_communication_structs(mesh);

  /* free temp arrays */
  free(rank_start);
  free(nodeload);
}
