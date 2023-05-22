/* load.c */
/* Wolfgang Tichy, 2/2019 */

#include "nmesh.h"
#include "amr.h"

#define PR 0


/* use timings, some MPI datatypes, amr vars */
extern tTiming Timing[1];
extern tnMPIvars nMPIvars[1];


/**********************************************************************/
/* functions to move nodes between procs */
/**********************************************************************/

/* simplistic load balancing */
void simple_load_balance(tMesh *mesh)
{
  load_balance(mesh, 1);
}


/* load balancing, where we can choose the balancing strategy:
   strategy = LOADBAL_SIMPLE ... */
void load_balance(tMesh *mesh, int strategy)
{
  /* for now we keep sibling 1-7 together with sibling 0 */
//FIXME: put the next line back in:
//  Timing->sibl1to7_weight = 0.;

  /* when we move elms much of in mesh->nbmesh will become wrong */
  amr_remove_mesh_nbelm(mesh);

  /* move elms bewteen ranks */
  load_balance_elms(mesh);

  /* update mesh->myelm */
  alloc_and_set_mesh_myelm(mesh);

  /* set fnb */
  amr_erase_all_elm_fnb(mesh);
  amr_elm_nbinfo_to_elm_fnb(mesh);
  //amr_elm_nbinfo_set_nnbinfo_mesh(mesh, 1); //make nnbinfo positive
  //NOTE:  ^--this call meses up nbinfo!!! ===> FIXME: remove these lines

  /* set elm->n and elm->pt_typ for the elms of mesh->nbmesh */
//FIXME: put the next line back in:
  amr_get_nbelm_elmheaders(mesh);
}

/* function that can be scheduled in LOADBALANCING */
int load_balance_if_needed(tMesh *mesh)
{
  int amr_load_balance = Par("amr_load_balance");

  if(Getv(amr_load_balance, "yes"))
  {
    load_balance(mesh, 1);
  }
  return 0;
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
      int vt = MeshVarType(mesh, vi);
      *ndoubles += (dat->v[vi]->N) * (vt!=AUXVAR);
      if(PR) { PRF;printf(": vi=%d ndoubles=%d\n", vi, *ndoubles); }
      nvars++;
    }
    //else
    //{
    //  if(PR) { PRF;printf(": vi=%d no dat\n", vi); }
    //}

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
      int vt = MeshVarType(mesh, vi);
      N = (dat->v[vi]->N) * (vt!=AUXVAR);
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
    if(N > dat->v[vi]->N) redim_array(dat->v[vi], N,1,1); //CHECK
    memcpy(dat->v[vi]->d, buf+bi, N * sizeof(double));
    bi += N;
    if(PR) { PRF;printf(": vi=%d bi=%d\n", vi, bi); }
  }

  /* now put the end of the buffer in dat->info */
  sizeofinfo = buf[bi++];
  memcpy(dat->info, buf+bi, sizeofinfo);

  /* make sure var amr_elm_nbinfo has the correct dimensions */
  amr_elm_nbinfo_redim_according_to_nnbinfo(node);

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
                        Node_eid(node), node->datrank, rank, desrank);
             fflush(stdout); }
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
      nMPI_Isend_double_com(scom, rq, other, Node_eid(node), WORLD);
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
        nMPI_Irecv_double_com(rcom, rq, other, Node_eid(node), WORLD);
      }

      /* allocate space already and init some stuff */
      if(node->dat) errorexit("destination node should not have dat yet");
      node->dat = alloc_dat(node);

      if(0)
      {
        PRF;printf(": nid%ld rank%d node->dat=%p\n",
                   Node_eid(node), rank, (void *) node->dat);
      }
      if(PR) { PRF;printf(": calling coordinates_init_node\n"); }
      coordinates_init_node(node);
    }
    if(PR) fflush(stdout);
  }
  else /* retrieve data from buffers */
  {
    if(PR) { PRF;printf(": nid%ld datrank%d rank%d desrank%d\n",
                        Node_eid(node), node->datrank, rank, desrank);
             fflush(stdout); }
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
  double myspeed = timing_get_mm_speed(mesh);
  double avspeed;
  int rank;

  for(rank=0; rank<size; rank++)
  {
    if(rank == nMPI_rank()) speed[rank] = myspeed;
    nMPI_Bcast(&(speed[rank]),1, nMPI_DOUBLE, rank);
  }

  avspeed = 0.;
  for(rank=0; rank<size; rank++) avspeed += speed[rank];
  avspeed /= size;

  //PRFs(":\n");
  //for(rank=0; rank<size; rank++)
  //  printf("speed[%d]=%g ", rank, speed[rank]);

  return avspeed;
}


/*************************************************************************/
/* load balance for elms */
/*************************************************************************/

/* fill in tElmfl myfl[1] with my first and last elm */
void get_my_Elmfl(tMesh *mesh, tElmfl myfl[1])
{
  long nelms = mesh->nmyelm;

  myfl->nelms = nelms;
  if(nelms <= 0) return;
  myfl->elm_fl[0] = *(mesh->myelm[0]);        /* shallow copies */
  myfl->elm_fl[1] = *(mesh->myelm[nelms-1]);
}

/* exchange first and last elms with rank+1 and rank-1 */
void get_nbr_rank_info(tMesh *mesh)
{
  int size = nMPI_size();
  int rank = nMPI_rank();
  tCom *com;
  int rq;
  tElmfl myfl[1];
  tElmfl *fl_m1 = mesh->nbr->fl_m1;
  tElmfl *fl_p1 = mesh->nbr->fl_p1;

  /* get my first and last elm from mesh->myelm */
  get_my_Elmfl(mesh, myfl);

  /* for MPI data transfers */
  //FIXME: why sizeof(double)???, should 2nd arg be 0 or 1????
  //com = alloc_com(sizeof(double), 0);
  com = alloc_com(sizeof(char), 0);

  //alloc_com is stupid!!! Its 1st arg should always be sizeof(void *)

  /* send myfl to rank-1 and also receive fl_m1 from rank-1 */
  if(rank>0)
  {
    rq = append_buffers_to_com(com, myfl,sizeof(myfl[0]),
                                    fl_m1,sizeof(fl_m1[0]));
    nMPI_Isend_Irecv_com(com, rq, nMPI_CHAR, rank-1, -1,+1, WORLD, WORLD);
  }

  /* send myfl to rank+1 and also receiv fl_p1 from rank+1 */
  if(rank < size-1)
  {
    rq = append_buffers_to_com(com, myfl,sizeof(myfl[0]),
                                    fl_p1,sizeof(fl_p1[0]));
    nMPI_Isend_Irecv_com(com, rq, nMPI_CHAR, rank+1, +1,-1, WORLD, WORLD);
  }

  /* wait until all sent and received */
  nMPI_Waitall_com(com);
  free_com(com);
}




/* Move data (dat) for elms that have been moved.
   Here we assume:
   * no dat has been moved yet
   * the list in mesh->myelm_head has been updated to contain only the elms
     we want on this rank (but some have no dat yet)
   * the list in mesh->myelm contains the elms that we had before the update,
     and thus has all the elms for which we still have dat
   * the elms whose dat that each proc sends away have the receiver rank in
     elm->dat->info->desrank
   * the elms whose dat each proc recvs have the sender rank in elm->datrank
*/
void load_exchange_dat_after_moving_elms(tMesh *mesh)
{
  int rank = nMPI_rank();
  struct list_head *pos;
  tCom *scom, *rcom;
  long i;

  PRFs(":\n");

  /* for MPI data transfers */
  scom = alloc_com(sizeof(double), 1);
  rcom = alloc_com(sizeof(double), 1);

  /* loop over mesh->myelm and send dat from elms that have a
     dat->info->desrank different from my rank */
  for(i=0; i<mesh->nmyelm; i++)
  {
    int desrank;
    tElm *elm = mesh->myelm[i];
    tDat *dat = elm->dat;
    if(!dat) errorexit("this elm must have dat");

    desrank = dat->info->desrank;

    /* fill MPI send buffers */
    if(rank != desrank)
      move_node_to_rank(elm, desrank, scom, rcom, 1);
  }

  /* loop over mesh->myelm_head list and recv dat from elms that have a
     datrank different from my rank */
  list_for_each(pos, &mesh->myelm_head)
  {
    tElm *elm = list_entry(pos, tElm, list);
    int datrank = elm->datrank;
    tDat *dat = elm->dat;

    /* redundant error checks */
    if(!dat)
    { if(datrank == rank) errorexit("dat=NULL but datrank=rank"); }
    else
    { if(datrank != rank) errorexit("dat!=NULL but datrank!=rank"); }

    /* setup MPI recv buffers */
    if(datrank != rank)
      move_node_to_rank(elm, rank, scom, rcom, 1);
  }

  /* wait for MPI sends and recvs */
  nMPI_Waitall_com_send(scom);
  free_com(scom);  /* free scom with all its buffers */
  nMPI_Waitall_com_recv(rcom);

  /* get var data out of recv buffers */
  set_com_counters(rcom, 0,0);
  list_for_each(pos, &mesh->myelm_head)
  {
    tElm *elm = list_entry(pos, tElm, list);
    if(elm->datrank != rank)
      move_node_to_rank(elm, rank, scom, rcom, 0);
  }
  free_com(rcom);  /* free rcom with all its buffers */

  /* free all elms in mesh->myelm, that now are on another rank */
  for(i=0; i<mesh->nmyelm; i++)
  {
    int desrank;
    tElm *elm = mesh->myelm[i];
    tDat *dat = elm->dat;
    if(!dat) errorexit("this elm must have dat");

    desrank = dat->info->desrank;

    /* remove elms that we don't have anymore */
    if(rank != desrank)
    {
      free_elm(elm);
      mesh->myelm[i] = NULL;
    }
  }
}

/* Info about load balancing of elms:
   ==================================
   *in each elm we have:
    t_elm = elm->dat->info->load_TimeIn_s

   *on one rank r we have:
    T_r = time to do something (e.g. a RK step) on all elms on this rank
    T_r = \sum_elm t_elm
    ops_r = # of operations actually done on this rank
          = v_r * T_r -> Timing->myops
            here: v_r = Timing->mm_speed -> speed[r]

   *total ops:
    allops = \sum_r ops_r -> Timing->allops

   *starting point of ops on each rank r:
    ops0_r = \sum_{ri=0}^{r-1} ops_ri -> ops0

   *distribution weights:
    w_r = v_r / (\sum_ri v_ri)

   *let's give each rank this many ops
    ops_bal_r = w_r * allops

   *we construct:
    ops_bal_sum[r] = \sum_{ri=0}^r ops_bal_ri

   *On rank r we want to have all elms that have an
      ops(elm) := ops0 + v_r \sum_{elm'=1st_elm}^elm t_elm'
    between ops_bal_sum[r-1] and ops_bal_sum[r].
   *All elms that are not within this range will be sent to the rank
    where they fall into this range. All elms on other ranks that fall in
    this range will be revcd from the other rank.

   *We first exchange the numbers to be sent and recvd (ns_elms[r] and
    nr_elms[r]) per rank r.
   *Once we have ns_elms[r] and nr_elms[r] we then exchange the elms via MPI.
   *After this we use load_exchange_dat_after_moving_elms (which uses the
    old move_node_to_rank) to exchange the dat between the ranks.
*/

/* comparison function for load_desired_rank */
int load_cmp_ops_bal_sum(const void *key, const void *ar, void *arg)
{
  double ops_bal_sum = *((const double *) ar);
  double ops_elm_sum = *((const double *) key);
  double max = *((double *) arg);
  double diff = (ops_elm_sum - ops_bal_sum)/max;
  int ret = diff * (INT_MAX/4);
  //printf("ops_bal_sum=%g ops_elm_sum=%g max=%g  diff=%g ret=%d\n",
  //       ops_bal_sum, ops_elm_sum, max, diff, ret);
  return ret;
}

/* find rank we want to put elm on
   In: size,ops_bal_sum  =>  Out: ops_elm_sum */
int load_desired_rank(int size, const double *ops_bal_sum, double ops_elm_sum)
{
  const void *val;
  double max = ops_bal_sum[size-1];
  size_t off, num;
  off = 0;
  num = size;

  //printf("off=%zu num=%zu\n", off, num);
  //printf("base: ops_bal_sum[0]=%g\n", ops_bal_sum[0]);
  ////printf("base: ops_bal_sum[1]=%g\n", ops_bal_sum[1]);
  //printf("key: ops_elm_sum=%g\n", ops_elm_sum);
  val = bisectionsearch(&ops_elm_sum, ops_bal_sum, &off, &num,
                        sizeof(ops_bal_sum[0]), load_cmp_ops_bal_sum, &max);
  ////printf("val=%p off=%zu num=%zu\n", val, off, num);
  if(val)
  {
    if(num == 2) return off+1;
    else         return off;
  }
  else
  {
    if(off == 0) return 0;
    else         return size-1;
  }
}

/* main load balancing function for elms:
   -We first determine how many elms (ns_elms[rk]) we send to another rank rk.
   -Then we tell the other ranks about it, and find out how many we will
    receive from them.
   -Then we pack the els we want to send into s_elms[rk], and remove them
    from our mesh->myelm_head list.
   -Then we send and recv the elm-headers (tElm0 type).
   -Finally we need to exchange the dat, for all elms that are now supposed
    to be elsewhere. */
void load_balance_elms(tMesh *mesh)
{
  int size = nMPI_size();
  int rank = nMPI_rank();
  int rk, torank;
  double avspeed, myspeed;
  double *speed = NULL;
  struct list_head *pos, *sav;
  double ops0, allops;
  double *ops_bal_sum = NULL;
  double myT = 0.;
  double w;
  tCom *scom, *rcom;
  ulong *ns_elms, *nr_elms; // number of elms to send or recv for each rank
  tElm0 **s_elms, **r_elms; // s_elms[3][7] elm7 to be sent to rank3 */
  //ulong nkeep; // number of elms we keep on this rank

  /* get how ops are currently distributed */
  timing_set_myops_ops0_allops(mesh);
  //printTiming();
  ops0   = Timing->ops0;
  //myops  = Timing->myops;
  allops = Timing->allops;

  /* get speeds on all ranks */
  ops_bal_sum = calloc(size, sizeof(ops_bal_sum[0]));
  speed       = calloc(size, sizeof(speed[0]));
  if(!speed || !ops_bal_sum)
    errorexit("no memory for speed or ops_bal_sum");
  avspeed = load_set_speed_array(mesh, speed);
  myspeed = speed[rank];

  /* ops needed for load balance */
  w = speed[0]/(avspeed*size); /* weight for rank0 */
  ops_bal_sum[0] = w * allops; /* ops rank0 should have for balance */
  for(rk=1; rk<size; rk++)
  {
    w = speed[rk]/(avspeed*size); /* weight for rank rk */
    /* sum over ops that rank 0 to rank rk should have */
    ops_bal_sum[rk] = ops_bal_sum[rk-1] + w * allops;
  }

  /* we do not need speed array any longer */
  free(speed);

  /* memory for number of elms we send to or recv from each rank */
  ns_elms = calloc(size, sizeof(ns_elms[0]));
  nr_elms = calloc(size, sizeof(nr_elms[0]));
  if(!ns_elms || !nr_elms) errorexit("no memory for ns_elms or nr_elms");

  ///* get boundaries op0 and op1 into which ops_bal has to fall
  //   within allops */
  //op1 = ops_bal_sum[rank];
  //if(rank>0) op0 = op1 - ops_bal_sum[rank-1];
  //else       op0 = 0.;

  /* find all elms that are not within my boundaries */
  torank = -1;
  list_for_each(pos, &mesh->myelm_head)
  {
    double et;
    int desrank;
    tElm *elm = list_entry(pos, tElm, list);
    tDat *dat = elm->dat;
    if(!dat) errorexit("this elm must have dat");

    et = timing_get_elm_load_TimeIn_s(elm);
    myT += et;
    desrank = load_desired_rank(size, ops_bal_sum, ops0 + myT*myspeed);
    //printelm(elm);
    //printf("desrank=%d\n", desrank);
    if(desrank != torank)
      torank = desrank;
    ns_elms[torank] += 1;
    dat->info->desrank = torank; /* store rank where this elm should go */
  }
  /* I don't send to myself, so zero ns_elms[rank] */
  //nkeep = ns_elms[rank];
  ns_elms[rank] = 0;

  /* from here on ops_bal_sum is no longer needed */
  free(ops_bal_sum);

  /* tell rank rk that I will send it ns_elms[rk] elms, and
     recv from rank rk how many (nr_elms[rk]) I will get */
  scom = alloc_com(sizeof(long), 0);
  rcom = alloc_com(sizeof(long), 0);
  for(rk=0; rk<size; rk++)
    if(rk != rank)
    {
      int rq;

      /* tell that I send ns_elms to others */
      rq = append_buffers_to_com(scom, &(ns_elms[rk]),1, NULL,0);
      nMPI_Isend_com(scom, rq, nMPI_UNSIGNED_LONG, rk, 100, WORLD);

      /* recv nr_elms from others */
      rq = append_buffers_to_com(rcom, NULL,0, &(nr_elms[rk]),1);
      nMPI_Irecv_com(rcom, rq, nMPI_UNSIGNED_LONG, rk, 100, WORLD);
    }

  /* alloc s_elms[rk] */
  s_elms = rows_calloc(size, ns_elms, sizeof(tElm0));

  /* set s_elms that has all elms that are not within my boundaries */
  torank = -1;
  list_for_each_safe(pos, sav, &mesh->myelm_head)
  {
    int desrank;
    ulong i;
    tElm *elm = list_entry(pos, tElm, list);
    tDat *dat = elm->dat;
    if(!dat) errorexit("this elm must have dat");

    desrank = dat->info->desrank;
    if(desrank != torank)
    {
      torank = desrank;
      i = 0;
    }
    if(torank != rank)
    {
      memcpy(&(s_elms[torank][i]), elm, sizeof(tElm0));
      i++;
      /* Remove this elm from list. But elm is still in mesh->myelm. */
      list_del(&elm->list);
    }
  }

  /* now send all the elms that we have told about, to the other ranks */
  for(rk=0; rk<size; rk++)
    if(rk != rank)
    {
      if(ns_elms[rk])
      {
        int rq;
        rq=append_buffers_to_com(scom, &(s_elms[rk][0]),ns_elms[rk], NULL,0);
        nMPI_Isend_com(scom, rq, nMPIvars->TELM0, rk, 200, WORLD);
      }
    }

  /* wait for recvs in rcom */
  nMPI_Waitall_com_recv(rcom);

  /* make new rcom */
  free_com(rcom);
  rcom = alloc_com(sizeof(long), 0);

  //Yo(100);
  //for(rk=0; rk<size; rk++)
  //  printf("rk=%d ns_elms[rk]=%zu nr_elms[rk]=%zu\n", rk, ns_elms[rk], nr_elms[rk]);
  //fflush(stdout);

  /* alloc r_elms[rk] */
  r_elms = rows_calloc(size, nr_elms, sizeof(tElm0));

  /* recv all the elms that others have sent, i.e. receive all that I have
     been told about */
  for(rk=0; rk<size; rk++)
    if(rk != rank)
    {
      if(nr_elms[rk])
      {
        int rq;
        rq=append_buffers_to_com(rcom, NULL,0, &(r_elms[rk][0]),nr_elms[rk]);
        //nMPI_Irecv_com(rcom, rq, nMPIvars->TELM0, rk, 100, WORLD);
        nMPI_Irecv_com(rcom, rq, nMPIvars->TELM0, rk, 200, WORLD);
      }
    }

  /* wait for recvs in rcom */
  nMPI_Waitall_com_recv(rcom);
  free_com(rcom);

  /* wait for sends in scom, then free all send related stuff */
  nMPI_Waitall_com_send(scom);
  free_com(scom);
  free(ns_elms);
  rows_free(s_elms, size);

  /* insert r_elms after current end of list */
  for(rk=rank+1; rk<size; rk++)
  {
    ulong i;
    for(i=0; i<nr_elms[rk]; i++)
    {
      tElm *elm = alloc_elm_init_pat(mesh, r_elms[rk][i].eploc->p); /* fresh elm */
      memcpy(elm, &(r_elms[rk][i]), sizeof(tElm0)); /* init elm from r_elms[i] */
      /* now add elm to the end of list in mesh */
      list_add_tail(&elm->list, &mesh->myelm_head);
    }
  }

  /* insert r_elms before current beginning of list */
  pos = &mesh->myelm_head; /* position where we insert */
  for(rk=0; rk<rank; rk++)
  {
    ulong i;
    for(i=0; i<nr_elms[rk]; i++)
    {
      tElm *elm = alloc_elm_init_pat(mesh, r_elms[rk][i].eploc->p); /* fresh elm */
      memcpy(elm, &(r_elms[rk][i]), sizeof(tElm0)); /* init elm from r_elms[i] */

      /* now add elm after pos */
      list_add(&elm->list, pos);
      pos = &elm->list; /* move insert position by one */
    }
  }

  /* free all received elm0 */
  free(nr_elms);
  rows_free(r_elms, size);

  /* free surfaces & indc since they will change now anyway */
  evolve_free_communication_structs(mesh);

  //Yo(200);
  //PRF;printf(": %zu in mesh->myelm_head\n", list_count_nodes(&mesh->myelm_head));
  //printf("  mesh->nmyelm=%lu\n", mesh->nmyelm);

  /* move dat to correct ranks now */
  load_exchange_dat_after_moving_elms(mesh);

  //alloc_and_set_mesh_myelm(mesh);
  //NOTE: update_mesh_myln_node_nid call causes an update of mesh->myelm

  //FIXME: adapt  update_mesh_myln_node_nid
  update_mesh_myln_node_nid(mesh);
  PRF;printf(": --> %lu on this proc\n", mesh->nmyelm);
  //printf(": --> %zu on this proc\n", list_count_nodes(&mesh->myelm_head));


  //FIXME: call function that set's up elm->fnb and such...
  //       maybe also update_mesh_myln_node_nid ???

  /* now that nodes are elsewhere re-init surfaces & indc */
  evolve_init_communication_structs(mesh);
}
