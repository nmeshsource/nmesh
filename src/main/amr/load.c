/* load.c */
/* Wolfgang Tichy, 2/2019 */

#include "nmesh.h"
#include "amr.h"

#define PR 0


/* use timings */
extern tTiming Timing[1];


/* object we pass around to figure out the desired rank of a node */
typedef struct tLOADINFO {
  long nid;
  long nnodes;
  int size;
  tNlist **rank_start; /* rank_start[i] is one of mesh->lns for rank i,
                          set in load_set_desired_rank_start */
} tLoadinfo;


/**********************************************************************/
/* functions to move nodes between procs */
/**********************************************************************/

/* compute desired rank */
int desiredrank_simple(tLoadinfo *li)
{
  double N = li->nnodes;
  double s = li->size;
  double nperproc = N/s;
  double desrank = li->nid/nperproc;
  return desrank;
}

/* simplistic load balancing using desiredrank_simple */
void simple_load_balance(tMesh *mesh)
{
/*
tNlist *elem;
fornodelist(mesh->lns, elem)
{
tNode *node = elem->node;
tDat *dat = node->dat;
if(dat) dat->info->load_TimeIn_s = 1./(node->nid+1);
}
//Yo(10);
//printmesh(mesh);
load_balance(mesh, LOADBAL_NODETIMES);
//Yo(20);
//printmesh(mesh);
return;
*/
  load_balance(mesh, LOADBAL_SIMPLE);
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

/* send around timing results and save them in array speed[] */
double load_set_speed_array(tMesh *mesh, double *speed)
{
  int size = nMPI_size();
  double myspeed = Timing->mm_speed;
  double speedmin = 1e-50;
  int rank;

  /* in case we forgot to measure Timing->mm_speed, just set myspeed=1 */
  if(myspeed <= speedmin) myspeed = 1.;

  for(rank=0; rank<size; rank++)
  {
    if(rank == nMPI_rank()) speed[rank] = myspeed;
    nMPI_Bcast(&(speed[rank]),1, nMPI_DOUBLE, rank);
  }

  PRFs(":\n");
  for(rank=0; rank<size; rank++)
    printf("speed[%d]=%g ", rank, speed[rank]);

  return myspeed;
}


/* Set array nodeload[] that contains the measured time each node used,
   Also return the sum over nodeload[] */
double load_set_nodeload_array(tMesh *mesh, const double *speed,
                               double *nodeload)
{
  int rank = nMPI_rank();
  double myspeed;
  tNlist *elem;
  double loadmin = 1e-50;
  double loadsum = 0.;

  if(speed) myspeed = speed[rank];
  else      myspeed = 1.;

  /* we assume that mesh->lns has the same order for all ranks */
  fornodelist(mesh->lns, elem)
  {
    tNode *node = elem->node;
    tDat *dat = node->dat;
    int datrank = node->datrank;
    long nid = node->nid;
    double load = loadmin;

    /* we need to broadcast nodeload from my nodes to all other ranks */
    if(dat) load = node->dat->info->load_TimeIn_s * myspeed;

    /* in case we forgot to measure the times, just set load to a very small
       uniform number: */
    if(load < loadmin) load = loadmin;

    nMPI_Bcast(&load,1, nMPI_DOUBLE, datrank);
    nodeload[nid] = load;
    loadsum += load;
  }

  //PRFs(":\n");
  //for(long l=0; l<mesh->nln; l++)
  //  printf("nodeload[%ld]=%g ", l, nodeload[l]);

  return loadsum;
}

/* sum over leaf nodes (starting at ln0) until loadsum reaches the value
   desired_load:
   Return: leaf after the one where we reached desired_load.
   In: ln0, desired_load[], nodeload.  Out: actual_load */
tNlist *inc_leaf_until_desired_loadsum(tNlist *ln0, const double *speed,
                                       int rank, double desired_loadsum,
                                       const double nodeload[],
                                       double *actual_loadsum)
{
  tNlist *elem;
  double rankspeed, sum;

  if(speed) rankspeed = speed[rank];
  else      rankspeed = 1.;

  sum = 0.;
  fornodelist(ln0, elem)
  {
    tNode *node = elem->node;
    long nid = node->nid;

    sum += nodeload[nid] / rankspeed;

    if(sum >= desired_loadsum) break;
  }

  *actual_loadsum = sum;
  if(elem) return elem->next;
  else     return NULL;
}

/* Set rank_start[i] array.
   It contains the first leaf that rank i should have */
void load_set_desired_rank_start(tMesh *mesh, const double *speed,
                                 double desired_loadsum,
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
    ln1 = inc_leaf_until_desired_loadsum(ln0, speed, rank, desired_loadsum,
                                         nodeload, &actual_loadsum);
  }
  //PRFs(":\n");
  //for(rank=0; rank<size; rank++)
  //  printf("rank_start[%d]=nid%ld ", rank, rank_start[rank]->node->nid);
}

/* compute desired rank */
int load_get_desiredrank(tLoadinfo *li)
{
  int rank;

  for(rank = 0; rank < li->size-1; rank++)
  {
    tNlist *ln1 = li->rank_start[rank+1];
    tNode *node = ln1->node;
    long nid1 = node->nid;
    //PRF;printf(": rank=%d nid=%ld nid1=%ld\n", rank, li->nid, nid1);
    /* we assume that leaf node nids are assigned in ascending order */
    if(li->nid < nid1) break;
  }
  if(rank >= li->size)
    errorexit("could not find the rank that should have li->nid");

  return rank;
}


/* load balancing, where we can choose the balancing strategy:
   strategy = LOADBAL_SIMPLE, LOADBAL_NODETIMES, ... */
void load_balance(tMesh *mesh, int strategy)
{
  long nnodes = mesh->nln;
  int size = nMPI_size();
  tLoadinfo li[1];
  int (*desiredrank)(tLoadinfo *li); /* func. pointer for distrib. strategy*/
  int desrank;
  tNlist *elem;
  tNode *node;
  tCom *scom;
  tCom *rcom;
  double totalload;            //only for LOADBAL_NODETIMES
  double *nodeload = NULL;     //only for LOADBAL_NODETIMES
  tNlist **rank_start = NULL;  //only for LOADBAL_NODETIMES
  double *speed = NULL;        //only for LOADBAL_NODETIMES_SPEEDS

  PRF;printf(": strategy=%d nnodes=%ld ", strategy, nnodes);

  /* set const part of li needed for all strategies */
  li->nnodes = nnodes;
  li->size   = size;

  /* set up stuff for each strategy */
  switch(strategy)
  {
  case LOADBAL_NODETIMES_SPEEDS:

    speed = calloc(size, sizeof(speed[0]));
    if(!speed)
    {
      free(speed);
      printf("  WARNING: quitting ");PRF;
      printf(" due to lack of memory!!!\n");
      /* do fallback? */
      //load_balance(mesh, LOADBAL_SIMPLE);
      return;
    }
    load_set_speed_array(mesh, speed);

  case LOADBAL_NODETIMES:

    /* load balancing based on measured node loads */
    nodeload   = dmalloc(nnodes);
    rank_start = calloc(size, sizeof(rank_start[0]));
    if(!nodeload || !rank_start)
    {
      free(rank_start);
      free(nodeload);
      printf("  WARNING: quitting ");PRF;
      printf(" due to lack of memory!!!\n");
      /* do fallback? */
      //load_balance(mesh, LOADBAL_SIMPLE);
      return;
    }

    /* get measured load for each node (in nodeload) and also totalload */
    totalload = load_set_nodeload_array(mesh, speed, nodeload);
    printf("totalload=%g\n", totalload);

    /* set array rank_start with 1st desired leaf node for each rank */
    load_set_desired_rank_start(mesh, speed, totalload/size, nodeload,
                                rank_start);

    /* set const part of li, and pick function to calc. desired rank */
    li->rank_start = rank_start;
    desiredrank = load_get_desiredrank;

    break;

  case LOADBAL_SIMPLE:

    /* pick function to calc. desired rank */
    desiredrank = desiredrank_simple;
    printf("\n");

    break;

  default:
    errorexit("unknown strategy");
  }

  /* for MPI data transfers */
  scom = alloc_com(sizeof(double), 1);
  rcom = alloc_com(sizeof(double), 1);

  /* free surfaces & indc since they will change now anyway */
  evolve_free_communication_structs(mesh);

  /* fill MPI send and recv buffers */
  fornodelist(mesh->lns, elem)
  {
    node = elem->node;
    li->nid = node->nid;
    desrank = desiredrank(li);
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
    li->nid = node->nid;
    desrank = desiredrank(li);
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
  free(speed);
}

/* function that can be scheduled in LOADBALANCING */
int load_balance_if_needed(tMesh *mesh)
{
  int amr_load_balance = Par("amr_load_balance");

  if(Getv(amr_load_balance, "timingbased"))
  {
    load_balance(mesh, LOADBAL_NODETIMES);
  }
  else if(Getv(amr_load_balance, "simple"))
  {
    load_balance(mesh, LOADBAL_SIMPLE);
  }
  return 0;
}
