/* surface.c */
/* Wolfgang Tichy, 2/2019 */

#include "nmesh.h"
#include "amr.h"

#define PR 0

/* use amr and DGglobals */
extern tAMR amr[1];
extern tDGglobals DGglobals[1];

/* functions to exchange surface data */

/* Note:
  Currently the exchange happens only between nearest neighbors.
  I.e. in e.g. the xy-plane data can be exchanged only between C
  and the neighbors 0,1,..., BUT not between C and NE
        ---------
        | . . . |            C is current node
  NW    | .3. . |  NE        0,1,2,3,4,5 are 5 nearest neighbors
        | . . . |            NW, NE, ... are next-to-nearest neighbors
-------------------------    the ". . ." indicate grid points inside
| . . . | . . . | . . . |
| .0. . | .C. . | .1. . |    C can exchange data with 0,1,2,3,4,5
| . . . | . . . | . . . |    C CANNOT exchange data with e.g. NW or NE, ...
-------------------------
        | . . . |
  SW    | .2. . | SE
        | . . . |
        ---------
*/

/**********************************************************************/
/* allocate and fill surfaces for vars that need it */
/**********************************************************************/
/* empty surface that we need to fill in */
tSurface *alloc_empty_surface(int nnb)
{
  tSurface *s = calloc(1, sizeof(*s));

  s->nnbsurf = nnb;
  s->nbsurf = calloc(nnb, sizeof(s->nbsurf[0]));
  s->allocd_nbsurf = calloc(nnb, sizeof(s->allocd_nbsurf[0]));

  return s;
}

/* free all we need to in a surface */
void free_surface(tSurface *s)
{
  int i;

  if(!s) return;

  /* free content of lists */

  /* free mysurf only if is allocd */
  if(s->allocd_mysurf) free_array(s->mysurf);

  /* free ajsurf only if is allocd */
  if(s->allocd_ajsurf) free_array(s->ajsurf);

  /* free nbsurf[i] only if it is allocd  */
  for(i=0; i<s->nnbsurf; i++)
    if(s->allocd_nbsurf[i]) free_array(s->nbsurf[i]);

  /* free lists */
  free(s->nbsurf);
  free(s->allocd_nbsurf);

  /* free surface */
  free(s);
}

/* free surfaces on node */
void free_all_surfaces(tNode *node)
{
  tDat *dat = node->dat;
  int vi,f;

  if(!dat) return;

  /* free all surfaces */
  for(f=0; f<6; f++)
  {
    for(vi=0; vi<dat->nv; vi++)
    {
      free_surface(dat->s[f][vi]);
      dat->s[f][vi] = NULL;
    }
  }
}

/* free all surfaces on all nodes in the mesh */
void free_all_myln_surfaces(tMesh *mesh)
{
  if(PR) PRFs(":\n");
  TIMER_START;
  formylnodes(mesh)
  {
    tNode *node = MyLnode;
    free_all_surfaces(node);
  }
  TIMER_STOP;
}


/* initialize a surface for var vi at face with nnb neighbors */
tSurface *init_surface(tNode *node, int face, int vi)
{
  int dir = face/2;
  int zones;
  tDat *dat = node->dat;
  int i, nfnb;
  tSurface *s;
  int n[3];
  int alloc_mysurf;

  /* do nothing if no data on this node */
  if(!dat) return NULL;
  nfnb = node->nfnb[face];

  /* do nothing if face has no neighbors */
  if(!nfnb) return NULL;

  /* do nothing if var is not enabled */
  if(!dat->v[vi]) return NULL;

  /* do nothing if ghost zone width is 0 for this var */
  zones = MeshVarSurfacezones(node->pat->mesh, vi);
  if(zones==0) return NULL;

  /* prep. */
  s = alloc_empty_surface(node->nfnb[face]);
  s->dat = dat;
  s->face = face;
  s->vi = vi;

  /* set n */
  for(i=0; i<3; i++) n[i] = node->n[i];
  n[dir] = zones;

  /* decide if we allocate mem for mysurf */
  if( (node->n[dir] == zones) && (dat->info->use_fv == 0) )
    alloc_mysurf = 0;
  else
    alloc_mysurf = 1;

  /* allocate or point my surface array */
  if(alloc_mysurf) s->mysurf = alloc_array(n);
  else             s->mysurf = dat->v[vi];
  s->allocd_mysurf = alloc_mysurf;

  return s;
}

/* init all sufaces of a node */
int init_all_surfaces(tNode *node)
{
  tDat *dat = node->dat;
  int face, vi, cnt;

  if(!dat) return 0;

  cnt=0;
  for(face=0; face<6; face++)
  {
    for(vi=0; vi<node->dat->nv; vi++)
    {
      free_surface(dat->s[face][vi]);
      dat->s[face][vi] = init_surface(node, face, vi);
      if(dat->s[face][vi]) cnt++;
    }
  }
  return cnt;
}

/* init all surfaces on all nodes in the mesh */
void init_all_myln_surfaces(tMesh *mesh)
{
  if(PR) PRFs(":\n");
  TIMER_START;
  formylnodes(mesh)
  {
    tNode *node = MyLnode;
    init_all_surfaces(node);
  }
  TIMER_STOP;
}

/**********************************************************************/
/* get data into these surfaces */
/**********************************************************************/
/* set mysurf array from my own data */
void set_mysurf(tSurface *s)
{
  tDat *dat = s->dat;
  tNode *node = dat->node;
  int f = s->face;
  int dir = f/2;
  int vi = s->vi;

  /* do nothing if this surface is NULL */
  if(!s) return;

  if(s->allocd_mysurf)
  {
    int p;
    int zones = MeshVarSurfacezones(node->pat->mesh, vi);

    if(zones<1) zones = 1;
    p = (node->n[dir] - zones) * (f%2); /* starting plane of surface */

    copy_array_planes(zones, dat->v[vi], dir, p, s->mysurf, 0);
  }
  else
  {
    s->mysurf = dat->v[vi]; //FIXME: this should be set only in init_surface
  }
}

/* set all mysurf of a node */
int set_all_mysurf(tNode *node)
{
  tDat *dat = node->dat;
  int face, vi, cnt;

  if(!dat) return 0;

  cnt=0;
  for(face=0; face<6; face++)
    for(vi=0; vi<node->dat->nv; vi++)
    {
      tSurface *s = dat->s[face][vi];
      if(s)
      {
        set_mysurf(s);
        cnt++;
      }
    }
  return cnt;
}

/* set all surfaces on all nodes in the mesh */
void set_all_myln_mysurf(tMesh *mesh)
{
  formylnodes(mesh)
  {
    tNode *node = MyLnode;
    set_all_mysurf(node);
  }
}


/***************************************************************************/
/* put mysurf into buffers and start MPI send/recv to get data into nbsurf */
/***************************************************************************/

/* count number of vars that have surfaces to be exchanged and set myN,
   input: node,my_f, nb,nb_f   output: nvars, vind, my_n, nb_n */
void find_nvars_vind_n_nbn(tNode *node, int my_f, tNode *nb, int nb_f,
                           int *nvars, int *vind, int my_n[3], int nb_n[3])
{
  tDat *dat = node->dat;
  int todo=1;
  int vi, i;
  int my_dir = my_f/2;
  int nb_dir = nb_f/2;

  /* set *vind to var "X" that has no surfaces */
  *vind=0;

  /* count number of vars that have surfaces to be exchanged and set myN */
  for(*nvars=0, vi=0; vi<dat->nv; vi++)
  {
    int zones = MeshVarSurfacezones(node->pat->mesh, vi);
    if(zones && dat->v[vi])
    {
      if(todo)
      {
        /* set nb_n */
        for(i=0; i<3; i++)
        {
          nb_n[i] = nb->n[i];
          my_n[i] = node->n[i];
        }
        nb_n[nb_dir] = my_n[my_dir] = zones;
        *vind = vi;
        todo=0;
      }
      (*nvars)++;
    }
  }
}

/* get all surfaces from neighbor with index ni at face */
void request_surfaces_exchange_for_all_vars(tNode *node, int face, int ni)
{
  tNode *nb = node->fnb[face][ni];
  tDat *dat = node->dat;
  int nb_f, nb_ni;
  int found, vi;

  /* do nothing if this node is on other proc */
  if(!dat) return;

  /* find face nb_f of nb that faces me */
  found = locate_facenb_in_fnbs(nb, node, &nb_f, &nb_ni);
  if(!found)
  {
    printf("node:\n");
    printelm(node);
    //printnfaces(node);
    printbfaces(node->pat);
    printf("neighbor node:\n");
    printelm(nb);
    //printnfaces(nb);
    printbfaces(nb->pat);
    errorexit("couldn't find nb face!!!");
  }

  /* is nb local? */
  if(nb->dat)
  {
    /* nb is local so just point s->nbsurf[ni] to its data */
    for(vi=0; vi<node->dat->nv; vi++)
    {
      tSurface *my_s = dat->s[face][vi];
      tSurface *nb_s = nb->dat->s[nb_f][vi];
      if(!my_s) continue;
      if(nb_s) my_s->nbsurf[ni] = nb_s->mysurf;
    }
  }
  else
  {
    /* nb is on other process so use MPI to exchange data */
    long s_ltag, r_ltag;
    int rq, nb_rank, s_tag, r_tag, lid, nb_lid, ci;
    nMPI_Comm s_comm, r_comm;
    tCom *com = dat->com[face];
    int nb_n[3], nb_N;
    int my_n[3], my_N;
    double *sbuf, *rbuf;  /* buffers for MPI */
    tArray *nbsurfs_f_ni; /* segmented array containing nbsurf[ni] for all vi */
    tSurface *s;
    int zones, nvars, cnt, si;

    /* count number of vars that have surfaces to exchange and set my_N */
    find_nvars_vind_n_nbn(node,face, nb,nb_f, &nvars, &vi, my_n, nb_n);
    nb_N = nb_n[0] * nb_n[1] * nb_n[2];
    my_N = my_n[0] * my_n[1] * my_n[2];

    /* use MPI to recv nb->dat->s[nb_f][vi]->mysurf in s->nbsurf[ni],
       and also send s->mysurf to nb->dat->s[nb_f][vi]->nbsurf[nb_ni] */
    nb_rank = nb->datrank;
    //s_tag = (node->nid)*6 + face;
    //r_tag = (nb->nid)*6 + nb_f;
    //r_comm = nb->comm;
    //s_comm = node->comm;
    lid = calc_node_lid(node);
    nb_lid = calc_node_lid(nb);
    //s_ltag = (lid*64 + ni)*6 + face;
    //r_ltag = (nb_lid*64 + nb_ni)*6 + nb_f;
    /* We do not use nb_ni and ni any more since the elms in mesh->nbelm
       have only the fnb needed for communication and not all of them.
       Thus ni and nb_ni do not necessarily agree on different ranks. */
    /* ALSO: We are sending the same data to all nbs at face face. Thus we
       can as well use the same s_tag! One day we should optimize this,
       and send the same data only once even if several nbs want it. */
    /* see also func request_indc_exchange_for_vl */
    s_ltag = lid*6 + face;
    r_ltag = nb_lid*6 + nb_f;
    nMPI_long_tag_to_commi_tag(s_ltag, &ci, &s_tag);
    s_comm = nMPIvars_get_comm(ci);
    nMPI_long_tag_to_commi_tag(r_ltag, &ci, &r_tag);
    r_comm = nMPIvars_get_comm(ci);

    /* alloc one segmented array for nbsurf of all nvars variables needed */
    s = dat->s[face][vi]; /* surface of 1st var that needs it on this face */
    if(!s) errorexit("1st var that needs surfaces, has no surface mem!");
    if(!s->nbsurf[ni])
      nbsurfs_f_ni = alloc_array_with_segs(nb_n, 0, nvars);
    else
      nbsurfs_f_ni = s->nbsurf[ni];

    /* alloc send and recv buffers */
    sbuf = calloc(nvars * my_N, sizeof(double));
    rbuf = nbsurfs_f_ni->d; /* use segmented array as rbuf */
    // it is good to use a long segmented array as rbuf
    // with dat->s[face][vi]->nbsurf[ni] pointing to the segments

    /* save buffers in com */
    rq = append_buffers_to_com(com, sbuf, nvars*my_N, rbuf,nvars*nb_N);

    /* fill send buffer */
    for(cnt=0, si=0, vi=0; vi<node->dat->nv; vi++)
    {
      s = dat->s[face][vi];
      zones = MeshVarSurfacezones(node->pat->mesh, vi);
      /* do nothing if var has no zones to exchange */
      if(s && zones && dat->v[vi])
      {
        /* point surface data to nbsurfs_f_ni=rbuf to later recv
           neighbor data */
        if(!s->nbsurf[ni])
        {
          s->nbsurf[ni] = get_array_seg(nbsurfs_f_ni, si);
          si++;                     // inc segment index
          s->allocd_nbsurf[ni] = 1; // flag that we allocd
        }

        /* save MPI request number in the array */
        s->nbsurf[ni]->info = rq;

        /* fill buffer for MPI exchange: sbuf[] = s->mysurf->d[] */
        memcpy(sbuf+cnt, s->mysurf->d, my_N * sizeof(sbuf[0]));
        cnt += my_N;
      }
    }
    /* now call MPI */
    if(0)
    {
      PRF;printf(": face=%d ni=%d nb_f=%d nb_ni=%d\n", face, ni, nb_f, nb_ni);
      printf("node=");printelm(node);
      printf("nb  =");printelm(nb);
      PRFs(": ");print_com_at(com, rq);
      printf(": rq=%d nb_rank=%d s_tag=%d r_tag=%d",
             rq, nb_rank, s_tag, r_tag);
      printf(" s_comm=");nMPI_print_Comm_name(s_comm);
      printf(" r_comm=");nMPI_print_Comm_name(r_comm);
      printf("\n");
    }
    MCK(
    nMPI_Isend_Irecv_double_com(com, rq, nb_rank, s_tag,r_tag, s_comm,r_comm)
    );
  }
}


/* put nbsurf from all faces and variables for this node in buffers */
void request_all_surfaces_exchange(tNode *node)
{
  int face, ni;
  tDat *dat = node->dat;

  /* do nothing if this node is on other proc */
  if(!dat) return;

  /* free req, send/recv arrays, before we start any send/recv */
  free_dat_reqs_after_Waitall_com_send(node);
  // NOTE: This Waitall in leads to a deadlock when it is called repeatedly
  // from the same node, as I would like to do in RK4.

  for(face=0; face<6; face++)
  {
    /* FIXME: set sendbuffer sbuf here already because
       request_surfaces_exchange_for_all_vars sends the some sbuf for all ni */

    for(ni=0; ni<node->nfnb[face]; ni++)
    {
      request_surfaces_exchange_for_all_vars(node, face, ni);
    }
  }

  /* signal that nb. surfaces on this node have not been received yet */
  dat->surfs_set = 0;
}

/* request surface exchanges on all my nodes in the mesh
   Note: The user has to call request_all_myln_surfaces_exchange! It calls
   request_all_surfaces_exchange on all nodes. If we call
   request_all_surfaces_exchange(n1) for only node n1, MPI deadlocks because
   the other nodes are not sending to n1 or receiving from n1 */
void request_all_myln_surfaces_exchange(tMesh *mesh)
{
  TIMER_START;

  /* If we want threads in this loop, we need MPI_Init_thread with
     MPI_THREAD_MULTIPLE, instead of MPI_THREAD_FUNNELED. */
  formylnodes_noomp(mesh)
  {
    tNode *node = MyLnode;
    request_all_surfaces_exchange(node);
  }
  TIMER_STOP;
}


/**********************************************************************/
/* get the nbsurf data out of the MPI buffers */
/**********************************************************************/

/* get all surfaces from neighbor with index ni at face */
void get_surfaces_for_all_vars(tNode *node, int face, int ni)
{
  tNode *nb = node->fnb[face][ni];
  tDat *dat = node->dat;
  tCom *com = dat->com[face];
  int nb_f, nb_ni;
  int found, vi;
  int nb_n[3], my_n[3];
  int nvars, rq;
  //int nb_N, zones, cnt;
  //double *rbuf;

  /* do nothing if this node is on other proc */
  if(!dat) return;

  /* do nothing if com is empty */
  if(com->n_rq == 0) return;

  /* find face nb_f of nb that faces me */
  found = locate_facenb_in_fnbs(nb, node, &nb_f, &nb_ni);
  if(!found) errorexit("couldn't find nb face!!!");

  /* if nb is local we have already exchanged info  */
  if(nb->dat) return;

  /* nb is on other process, we have used MPI to exchange data */

  /* count number of vars that exchanged surfaces and set my_N */
  find_nvars_vind_n_nbn(node,face, nb,nb_f, &nvars, &vi, my_n, nb_n);
  //nb_N = nb_n[0] * nb_n[1] * nb_n[2];
  /* do nothing if there are no vars that exchanged surfaces */
  if(!nvars) return;

  //PRF;printf(": nvars=%d\n", nvars);

  /* get MPI request number */
  rq = dat->s[face][vi]->nbsurf[ni]->info;

  /* wait for MPI buffer */
  loadtimer_pause(node);  /* we don't want to time MPI_Wait */
  MCK( nMPI_Wait_com_recv(com, rq) );
  loadtimer_resume(node); /* but we time everything else */

  /* find our recv buffer, and set it to NULL.
     We do this because s->nbsurf[ni]->d already points there, and because
     then realloc_com_reqs or free_com will not free s->nbsurf[ni]->d */
  set_com_recv_buf(com, rq, NULL);

//  /* get data out of recv buffer */
//  for(cnt=0, vi=0; vi<node->dat->nv; vi++)
//  {
//    tSurface *s = dat->s[face][vi];
//
//    zones = MeshVarSurfacezones(node->pat->mesh, vi);
//    if(!zones) continue; /* do nothing if var has no zones to exchange */
//    /* do nothing if var has no zones to exchange */
//    if(zones && dat->v[vi])
//    {
//      /* get neighbor data from buffer */
//      memcpy(s->nbsurf[ni]->d, rbuf+cnt, nb_N * sizeof(rbuf[0]));
//      cnt += nb_N;
//    }
//  }
}


/* free req and send arrays after all has been sent */
void free_dat_reqs_after_Waitall_com_send(tNode *node)
{
  int face;
  tDat *dat = node->dat;

  if(!dat) return;

  for(face=0; face<6; face++)
  {
    /* to be sure, wait again for all recvs */
    MCK( nMPI_Waitall_com_recv(dat->com[face]) );

    /* wait until all has been sent, then free all buffers for this face */
    MCK( nMPI_Waitall_com_send(dat->com[face]) );
    realloc_com_reqs(dat->com[face], 0); /* free req and send arrays */
  }
}

/* get nbsurf from all faces and variables for this node out of buffers */
void get_all_surfaces(tNode *node)
{
  int face, ni;
  tDat *dat = node->dat;

  if(!dat) return;

  /* set surfaces only once */
  if(dat->surfs_set) return;

  for(face=0; face<6; face++)
  {
    /* get nbsurf for each neighbor */
    for(ni=0; ni<node->nfnb[face]; ni++)
    {
      get_surfaces_for_all_vars(node, face, ni);
    }

    /* set ajsurf on this nodeface via interpolation */
    set_ajsurf_forall_vars(node, face);
    /* FIXME: if we use formylnodes_noomp to call get_all_surfaces, it might
       be better to later call set_all_myln_ajsurf instead of using the
       set_ajsurf_forall_vars above. */

    /* After set_ajsurf_forall_vars we could free nbsurf already */
    //FIXME: to conserve memory we should free nbsurf here!!!
    //free_nbsurf_only_forall_vars(node, face);
  }

  /* signal that nb. surfaces on this node have been received now */
  dat->surfs_set = 1;
}

/* get nbsurf for all nodes out of buffers and free the buffers */
void get_all_myln_surfaces(tMesh *mesh)
{
  TIMER_START;

  /* If we want threads in this loop, we need MPI_Init_thread with
     MPI_THREAD_MULTIPLE, instead of MPI_THREAD_FUNNELED. */
  formylnodes_noomp(mesh)
  {
    tNode *node = MyLnode;

    loadtimer_start(node);  /* time interp in set_ajsurf_forall_vars */
    get_all_surfaces(node);
    loadtimer_stop(node);
  }

  /* postpone Waitall until we have finished all nodefaces. This could have
     been already called in get_all_surfaces to free mem earlier.*/
  formylnodes_noomp(mesh)
  {
    tNode *node = MyLnode;
    free_dat_reqs_after_Waitall_com_send(node);
  }

  TIMER_STOP;
}


/* Issue MPI_Testall calls to test node->dat->com for completion.
   We do this only to cause MPI progression!
   sendrecv=0 tests only sends
   sendrecv=1 tests only recvs
   sendrecv=2 tests both sends and recvs */
int causeMPIprogress_all_myln_surfaces(tMesh *mesh, int sendrecv)
{
  int flag = 1;

  TIMER_START;

  /* If we want threads in this loop, we need MPI_Init_thread with
     MPI_THREAD_MULTIPLE, instead of MPI_THREAD_FUNNELED. */
  formylnodes_noomp(mesh)
  {
    tNode *node = MyLnode;
    tDat *dat = node->dat;
    int face;

    if(dat)
      for(face=0; face<6; face++)
      {
        tCom *com = node->dat->com[face];
        int fl;

        if(com)
        {
          switch(sendrecv)
          {
          case 0:
            MCK( nMPI_Testall_com_send(com, &fl) );
            break;
          case 1:
            MCK( nMPI_Testall_com_recv(com, &fl) );
            break;
          default:
            MCK( nMPI_Testall_com(com, &fl) );
          }
          flag = flag && fl;
        }
      }
  }

  TIMER_STOP;

  return flag;
}


/**********************************************************************/
/* set adjacent surface in ajsurf using data in nbsurf */
/**********************************************************************/
/* do two nodes have same bounding boxes, orthogonal to dir? */
int same_bbox_normal_to_dir(tNode *node1, tNode *node2, int dir)
{
  int d, samebb = 1;

  for(d=0; d<3; d++)
    if(d!=dir)
    {
      int ff = d*2;
      if(!dequal(node1->bbox[ff],   node2->bbox[ff]) ||
         !dequal(node1->bbox[ff+1], node2->bbox[ff+1])) { samebb=0; break; }
    }
  return samebb;
}

/* do two nodes have same point number n and pt_typ, orthogonal to dir? */
int same_n_and_pt_typ_normal_to_dir(tNode *node1, tNode *node2, int dir)
{
  int d, same_n_typ = 1;

  for(d=0; d<3; d++)
    if(d!=dir)
    {
      if( (node1->n[d] != node2->n[d]) ||
          (node1->pt_typ[d] != node2->pt_typ[d]) )
      {
        same_n_typ = 0;
        break;
      }
    }
  return same_n_typ;
}

/* do two nodes have same point number n, orthogonal to norm1 and norm2? */
int same_node_n_normal_to_norm1_2(tNode *node1, int norm1,
                                  tNode *node2, int norm2)
{
  int d1 = Dir1_norm(norm1);
  int d2 = Dir2_norm(norm1);
  int od1 = Dir1_norm(norm2);
  int od2 = Dir2_norm(norm2);

  if(node1->n[d1] != node2->n[od1])  return 0;
  if(node1->n[d2] != node2->n[od2])  return 0;

  return 1;
}


/* return the first surface tht is not NULL */
tSurface *first_nonNULL_surf_in_dat(tDat *dat, int f)
{
  int vi;
  tSurface *s = NULL;
  for(vi=0; vi<dat->nv; vi++)
  {
    s = dat->s[f][vi];
    if(s) break;
  }
  return s;
}


/* set ajsurf array from data in nbsurf on face f for all vars */
void set_ajsurf_forall_vars(tNode *node, int f)
{
  tMesh *mesh = Elm_mesh(node);
  tPat *pat = node->pat;
  int nnb = node->nfnb[f];
  int dir = f/2;
  //int p = (node->n[dir] - 1) * (f%2); /* plane of surface */
  tSurface *s1;
  int *s1_n;
  tNode *nb;
  int vi, ni, found, nb_f, nb_ni, nb_dir;
  int Cp_is_set;
  tArray *Cp[2], **Ip, **Res;
  tArray *(*Cb)[2];
  tDat *dat = node->dat;
  char str[100];
  int WENOorder;
  double (*interp1d_fv)(int k, double x, int np,
                        const double *x_p, const double *w_interp);
  if(!dat) return;

  TIMER_START;

  /* if there is only 1 neighbor we may not need to interpolate */
  if(nnb == 1)
  {
    nb = node->fnb[f][0];
    /* if we have only one neighbor on the same level in the same patch
       we may not need interpolation */
    if(nb->pat == pat)
    {
      if(Elm_l(nb) == Elm_l(node))
      {
        int same_n_t = same_n_and_pt_typ_normal_to_dir(nb, node, dir);
        /* if number of points is the same we can copy or just point
           ajsurf to nbsurf[0] */
        if(same_n_t)
        {
          for(vi=0; vi<dat->nv; vi++)
          {
            tSurface *s = dat->s[f][vi];
            /* do nothing if this surface is NULL */
            if(s) s->ajsurf = s->nbsurf[0];
          }
          goto end_set_ajsurf_forall_vars;
        }
      }
    }
    else /* neighbor is from diff patch */
    {
      tSurface *s0 = first_nonNULL_surf_in_dat(dat, f);

      found = locate_facenb_in_fnbs(nb, node, &nb_f, &nb_ni);
      if(!found) errorexit("couldn't find nb face!!!");

      /* if nb has only node as neighbor we may not need to interpolate */
      if(1 && nb->nfnb[nb_f] == 1 && s0->nbsurf[0]->N == s0->mysurf->N)
      {
        int *n = node->n;
        int pl = (n[dir]-1)*(f%2);
        int d1 = Dir1_norm(dir);
        int d2 = Dir2_norm(dir);
        int n1 = n[d1];
        int n2 = n[d2];
        int *nb_n = nb->n;
        int odir = nb_f/2;
        int opl = (nb_n[odir]-1)*(nb_f%2);
        int od1 = Dir1_norm(odir);
        int od2 = Dir2_norm(odir);
        int on1 = nb_n[od1];
        int on2 = nb_n[od2];

        if( ( node->pt_typ[d1] == nb->pt_typ[od1] &&
              node->pt_typ[d2] == nb->pt_typ[od2] ) ||
            ( node->pt_typ[d2] == nb->pt_typ[od1] &&
              node->pt_typ[d1] == nb->pt_typ[od2] ) )
        {
          //tMesh *mesh = node->pat->mesh;
          int ix = Ind("x");
          double *px[] = { Vard(node,ix), Vard(node,ix+1), Vard(node,ix+2) };
          int i,j,k, ind;
          int inb[3][3];
          double x0[3], dist[3], mx0;

          /* Find distance to closest nb point for 3 node points.
             [Note: if n[]={1,1,1}:
               -all 3 points are the same
               -all 3 lie in the node center, not at the surface!!!
               -distance to nb points cannot be zero!!!
                ==> interpolation will be used, even though we could copy
               +However, this is likely a rare case. With refinement most nodes
                will not be at patch boundaries. So don't worry for now.]
          */
          /* 3 points: */
          /* point 1 */
          ijk_inplaneN(dir, i,j,k, 0,0,pl);
          ind = Ind_n(i,j,k, n);
          x0[0] = px[0][ind];
          x0[1] = px[1][ind];
          x0[2] = px[2][ind];
          dist[0] = nearest_corner_of_xyz_inplaneN(nb, odir, opl, inb[0], x0);
          mx0 = magnitude_xyz(x0);
          if(mx0>0.) dist[0] = dist[0]/mx0; /* normalize dist */
          //printf("i,j,k: %d %d %d inb[0]: %d %d %d\n",
          //       i,j,k, inb[0][0],inb[0][1],inb[0][2]);

          /* point 2 */
          ijk_inplaneN(dir, i,j,k, n1-1,0,pl);
          ind = Ind_n(i,j,k, n);
          x0[0] = px[0][ind];
          x0[1] = px[1][ind];
          x0[2] = px[2][ind];
          dist[1] = nearest_corner_of_xyz_inplaneN(nb, odir, opl, inb[1], x0);
          mx0 = magnitude_xyz(x0);
          if(mx0>0.) dist[1] = dist[1]/mx0; /* normalize dist */
          //printf("i,j,k: %d %d %d inb[1]: %d %d %d\n",
          //       i,j,k, inb[1][0],inb[1][1],inb[1][2]);

          /* point 3 */
          ijk_inplaneN(dir, i,j,k, 0,n2-1,pl);
          ind = Ind_n(i,j,k, n);
          x0[0] = px[0][ind];
          x0[1] = px[1][ind];
          x0[2] = px[2][ind];
          dist[2] = nearest_corner_of_xyz_inplaneN(nb, odir, opl, inb[2], x0);
          mx0 = magnitude_xyz(x0);
          if(mx0>0.) dist[2] = dist[2]/mx0; /* normalize dist */
          //printf("i,j,k: %d %d %d inb[2]: %d %d %d\n",
          //       i,j,k, inb[2][0],inb[2][1],inb[2][2]);

          //pr3v("dist", dist);
          //printf("f%d i,j,k: %d %d %d dir%d pl%d nb_f%d odir%d opl%d\n",
          //       f, i,j,k, dir,pl, nb_f,odir,opl);
          if(PR)
          {
            //PRF;printf(": node->nid=%ld nb->nid=%ld  f%d nb_f=%d\n",
            //           node->nid, nb->nid, f, nb_f);
            PRF;printf(": %s ", nodename(node,str,99));
            printf("%s  f%d nb_f%d:", nodename(nb,str,99), f, nb_f);
          }

          /* if all 3 points coincide with grid points of nb we can copy */
          if(dist[0]<=dequaleps && dist[1]<=dequaleps && dist[2]<=dequaleps)
          {
            /* if axis are aligned in both nodes */
            if(inb[0][od1] == 0 && inb[0][od2] == 0)
            {
              if(inb[1][od1] == on1-1 && inb[1][od2] == 0)
              {
                if(inb[2][od1] == 0 && inb[2][od2] == on2-1)
                {
                  if(on1==n1 && on2==n2)
                  {
                    /* just copy by pointing to same array */
                    if(PR) printf(" point ajsurf to nbsurf[0]\n");
                    for(vi=0; vi<dat->nv; vi++)
                    {
                      tSurface *s = dat->s[f][vi];
                      /* do nothing if this surface is NULL */
                      if(s) s->ajsurf = s->nbsurf[0];
                    }
                    goto end_set_ajsurf_forall_vars;
                  }
                }
              }
              if(inb[1][od1] == 0 && inb[1][od2] == on2-1)
              {
                if(inb[2][od1] == on1-1 && inb[2][od2] == 0)
                {
                  if(on1==n2 && on2==n1)
                  {
                    /* copy with two axis interchanged */
                    if(PR) printf(" copy from nbsurf[0] with two axis interchanged\n");
                    copy_ajsurf_from_nbsurf0(node,f,nb_f, 1,0,0);
                    goto end_set_ajsurf_forall_vars;
                  }
                }
              }
            }

            /* if axis1 reversed */
            if(inb[0][od1] == on1-1 && inb[0][od2] == 0)
            {
              if(inb[1][od1] == 0 && inb[1][od2] == 0)
              {
                if(inb[2][od1] == on1-1 && inb[2][od2] == on2-1)
                {
                  if(on1==n1 && on2==n2)
                  {
                    /* copy with axis1 reversed */
                    if(PR) printf(" copy from nbsurf[0] with axis1 reversed\n");
                    copy_ajsurf_from_nbsurf0(node,f,nb_f, 0,1,0);
                    goto end_set_ajsurf_forall_vars;
                  }
                }
              }
              if(inb[1][od1] == on1-1 && inb[1][od2] == on2-1)
              {
                if(inb[2][od1] == 0 && inb[2][od2] == 0)
                {
                  if(on1==n2 && on2==n1)
                  {
                    /* copy with two axis interchanged and axis1 reversed */
                    if(PR) printf(" copy from nbsurf[0] with two axis interchanged and axis1 reversed\n");
                    copy_ajsurf_from_nbsurf0(node,f,nb_f, 1,1,0);
                    goto end_set_ajsurf_forall_vars;
                  }
                }
              }
            }

            /* if axis2 reversed */
            if(inb[0][od1] == 0 && inb[0][od2] == on2-1)
            {
              if(inb[1][od1] == on1-1 && inb[1][od2] == on2-1)
              {
                if(inb[2][od1] == 0 && inb[2][od2] == 0)
                {
                  if(on1==n1 && on2==n2)
                  {
                    /* copy with axis2 reversed */
                    if(PR) printf(" copy from nbsurf[0] with axis2 reversed\n");
                    copy_ajsurf_from_nbsurf0(node,f,nb_f, 0,0,1);
                    goto end_set_ajsurf_forall_vars;
                  }
                }
              }
              if(inb[1][od1] == 0 && inb[1][od2] == 0)
              {
                if(inb[2][od1] == on1-1 && inb[2][od2] == on2-1)
                {
                  if(on1==n2 && on2==n1)
                  {
                    /* copy with two axis interchanged and axis2 reversed */
                    if(PR) printf(" copy from nbsurf[0] with two axis interchanged and axis2 reversed\n");
                    copy_ajsurf_from_nbsurf0(node,f,nb_f, 1,0,1);
                    goto end_set_ajsurf_forall_vars;
                  }
                }
              }
            }

            /* if axis1 & 2 reversed */
            if(inb[0][od1] == on1-1 && inb[0][od2] == on2-1)
            {
              if(inb[1][od1] == 0 && inb[1][od2] == on2-1)
              {
                if(inb[2][od1] == on1-1 && inb[2][od2] == 0)
                {
                  if(on1==n1 && on2==n2)
                  {
                    /* copy with both axis reversed */
                    if(PR) printf(" copy from nbsurf[0] with both axis reversed\n");
                    copy_ajsurf_from_nbsurf0(node,f,nb_f, 0,1,1);
                    goto end_set_ajsurf_forall_vars;
                  }
                }
              }
              if(inb[1][od1] == on1-1 && inb[1][od2] == 0)
              {
                if(inb[2][od1] == 0 && inb[2][od2] == on2-1)
                {
                  if(on1==n2 && on2==n1)
                  {
                    /* copy with two axis interchanged and both reversed */
                    if(PR) printf(" copy from nbsurf[0] with both axis interchanged and both reversed\n");
                    copy_ajsurf_from_nbsurf0(node,f,nb_f, 1,1,1);
                    goto end_set_ajsurf_forall_vars;
                  }
                }
              }
            }
          }
          /* nothing matches if we get here and we need to interpolate */
        }
      }
      /* go on to interpolation ... */
    }
  }

  /* Ok if we get here we need interpolation */
  s1 = first_nonNULL_surf_in_dat(dat, f);
  /* do nothing if all sufaces are NULL */
  if(!s1) goto end_set_ajsurf_forall_vars;
  s1_n = s1->mysurf->n;

  /* array memory to store points of mysurf in X coords */
  Cp[0] = alloc_array(s1_n);
  Cp[1] = alloc_array(s1_n);
  Cp_is_set = 0;

  /* use data from all neighbors to interpolate into ajsurf */
  /* storage for nb coords and interp. results */
  Ip = calloc(nnb, sizeof(Ip[0]));
  Cb = calloc(nnb, sizeof(Cb[0]));
  Res = calloc(nnb, sizeof(Res[0]));
  if(!Ip || !Cb || !Res) errorexit("no memory for Ip, Cb, Res");

  /* 1. set neighbor coords within node surface */
  for(ni=0; ni<nnb; ni++)
  {
    nb = node->fnb[f][ni];
    found = locate_facenb_in_fnbs(nb, node, &nb_f, &nb_ni);
    if(!found) errorexit("couldn't find nb face!!!");
    nb_dir = nb_f/2;

    /* array memory to store points of neighbors in Xb coords */
    Ip[ni] = alloc_array(s1_n);
    Cb[ni][0] = alloc_array(s1_n);
    Cb[ni][1] = alloc_array(s1_n);
    Res[ni] = alloc_array(s1_n);

    /* all X,Y,Z coords are within same patch */
    if(nb->pat == pat)
    {
      if(!Cp_is_set)
      {
        /* points of mysurf */
        fill_2arrays_with_nodepoints(node, dir, Cp);
        /* convert Cp from Xb to X coords for node,
           these X are spread over the neighbor nodes */
        array_Xplane_of_Xb(node, dir, Cp, Cp);
        Cp_is_set = 1;
      }
      /* find points inside neigh. -> mask is returned in Ip */
      array_find_Xplane_in_node(nb,nb_dir, Cp, Ip[ni]);

      /* convert Cp to neighbor's internal basis coords */
      array_Xbplane_of_X(nb,nb_dir, Cb[ni], Cp);
    }
    else
    {
      /* get node points in neighbor coords */
      tBface *bface0 = pat->bfaces[f];
      int ioC0 = bface0->ioC0_0;
      tArray *oC[2];

      /* check if the neighbor coords are saved in ioC0... */
      if(ioC0>0)
      {
        oC[0] = VarA(node,ioC0+f);
        oC[1] = VarA(node,ioC0+6+f);
        /* find points inside neigh. -> mask is returned in Ip */
        mark_points_in_nb_f(node,f,Cp, nb,nb_f,oC, Ip[ni]);
        errorexit("mark_points_in_nb_f checks only 2 out of 3 coords!!!");
        if(Cp[0]->N != oC[0]->N)
          errorexit("Cp[0]->N != oC[0]->N");
      }
      else /* compute oC from Cp */
      {
        /* point oC to Cb[ni], to save oC in Cb[ni] */
        oC[0] = Cb[ni][0];
        oC[1] = Cb[ni][1];
        /* compute oC in node plane and find points inside neigh.
           -> mask is returned in Ip */
        array_find_nbXface_of_Xface(node,f, nb,nb_f, oC, Ip[ni]);
      }

      //printarray_int(Ip[ni]);
      if(0 && f==1 && Node_eid(nb)==7)
      {
        printf("ni=%d:\n", ni);
        printelm(node);
        printbfaces_on_f(node->pat, f);
        printelm(nb);
        printbfaces_on_f(nb->pat, nb_f);
        printarray_int(Ip[ni]);
        printarray(Cp[0]);
        printarray(oC[0]);
        printarray(Cp[1]);
        printarray(oC[1]);
      }

      /* convert oC to neighbor's internal basis coords */
      array_Xbplane_of_X(nb,nb_dir, Cb[ni], oC);
    }
  }

  /* 2. use interpolation to get vars from neighbors to node */
  /* choose 1d interpolator for basis_interp2d_toIpoints in fv case: */
  switch(DGglobals->fv_surface_interp_mode)
  {
  case FV_2DINTERP_LINEAR:
    interp1d_fv = basis_pw_linear;
    break;
  case FV_2DINTERP_PARAB:
    interp1d_fv = basis_pw_parab;
    break;
  case FV_2DINTERP_WENO:
    interp1d_fv = NULL;
    WENOorder = Geti(amr->WENO_interp_order);
    break;
  default:
    errorexiti("illegal value: DGglobals->fv_surface_interp_mode = %d",
               DGglobals->fv_surface_interp_mode);
  }
  /* now loop over all vars and interpolate */
  for(vi=0; vi<dat->nv; vi++)
  {
    tSurface *s = dat->s[f][vi];
    if(s)
    {
      int k, cnt;

      /* unless it has mem already, get mem for ajsurf */
      if(s->allocd_ajsurf==0 || s->ajsurf==NULL)
      {
        /* get mem for ajsurf, needs to be freed later */
        s->ajsurf = alloc_array(s->mysurf->n);
        s->allocd_ajsurf = 1;
      }

      /* interpolation within each nb, save results in Res */
      for(ni=0; ni<nnb; ni++)
      {
        int od1, od2;

        nb = node->fnb[f][ni];
        found = locate_facenb_in_fnbs(nb, node, &nb_f, &nb_ni);
        if(!found) errorexit("couldn't find nb face!!!");
        nb_dir = nb_f/2;
        /* For now we assume that pt_typ=P_UNIFORM means fin.vol., and that
           we thus use e.g. linear instead of Lagrange interpolation!
           Later we may want to check nb->dat->info->use_fv, but this needs
           to be send via MPI... */
        od1 = Dir1_norm(nb_dir);
        od2 = Dir2_norm(nb_dir);
        if(nb->pt_typ[od1]==P_UNIFORM && nb->pt_typ[od2]==P_UNIFORM)
        {
          if(DGglobals->fv_surface_interp_mode == FV_2DINTERP_WENO)
            interp2d_toIpoints(nb, s->nbsurf[ni], nb_dir,0, Cb[ni],Ip[ni],
                               WENOorder,INTERP_WENO,1., Res[ni]);
          else
            basis_interp2d_toIpoints(nb, s->nbsurf[ni], nb_dir,0,
                                     Cb[ni],Ip[ni], Res[ni], interp1d_fv);
        }
        else /* Lagrange interpolation is the default */
        {
          basis_interp2d_toIpoints(nb, s->nbsurf[ni], nb_dir,0,
                                   Cb[ni],Ip[ni], Res[ni], Lagrange_of_x);
        }

if(0 && Node_eid(nb)==17 && Node_eid(nb)==64 && vi==35)
{
//tMesh *mesh = node->pat->mesh;
int ll;
printelm(node);
printelm(nb);
printf("nb_f=%d nb_ni=%d nb_dir=%d    f=%d ni=%d\n",
        nb_f, nb_ni, nb_dir, f,ni);
printf("%s\n", VarName(vi));

forarray(Ip[ni], ll)
{
if(Ip[ni]->i[ll]>=0)
{
int iX = Ind("X");
double *pX[] = { Vard(node, iX), Vard(node, iX+1), Vard(node, iX+2) };
double nbCb[] = { Cb[ni][0]->d[ll], Cb[ni][1]->d[ll]};
double nbC[2], nbX[3], nbx[3], X[3], x[3];
int nbd1 = Dir1_norm(nb_dir);
int nbd2 = Dir2_norm(nb_dir);
X_of_Xb_indir(nb, nbd1, nbCb[0], &(nbC[0]));
X_of_Xb_indir(nb, nbd2, nbCb[1], &(nbC[1]));

X_from_C_on_face(nb->pat, nb_dir*2 , nbC, nbX);
set_xyz(nb->pat,0,-1, nbX, nbx);

// coords of point ll in face f:
int i,j, ijk;
switch(dir)
{
case 0:
k = kOfInd_n(ll, s->ajsurf->n);
j = jOfInd_n_k(ll, s->ajsurf->n, k);
i = (f%2)* (node->n[0]-1);
break;
default:
errorexit("arrgh. I thought we have only dir0 in our tests!!!");
}
ijk = Ind_n(i,j,k, node->n);
X[0] = pX[0][ijk];
X[1] = pX[1][ijk];
X[2] = pX[2][ijk];
set_xyz(node->pat,0,-1, X, x);

double res = Res[ni]->d[ll];
double mys = s->mysurf->d[ll];
if(fabs((res-mys)/mys) >1e-4)
{
printf("node->pat->p=%d ", node->pat->p);
pr3v("  X", X);
pr3v("  x", x);
printf("\n");
printf("nb->pat->p=%d   ", nb->pat->p);
pr3v("nbX", nbX);
pr3v("nbx", nbx);
printf("\n");
printf("ll=%d  i,j,k=%d %d %d   Ip=%d\n", ll, i,j,k, Ip[ni]->i[ll]);
printf("Res[ni]->d[ll]=%g s->mysurf->d[ll]=%g\n",
Res[ni]->d[ll], s->mysurf->d[ll]);
//exit(8);
}
}
}
//printarray(s->mysurf);
//printarray(s->nbsurf[ni]);
//exit(9);
}

      }

      /* take average of results from different nb */
      forarray(s->ajsurf, k)
      {
        cnt = 0;
        s->ajsurf->d[k] = 0.;
        for(ni=0; ni<nnb; ni++)
        {
          if(Ip[ni]->i[k] >= 0) /* if nb ni has the point */
          {
            cnt++; /* count num of neigh. who have this point */
            s->ajsurf->d[k] += Res[ni]->d[k];
          }
        }
        /* average if there was more than one nb with this point */
        if(cnt>1) s->ajsurf->d[k] /= cnt;
      }
    } /* end: if(s) */
  }

  /* 3. free temp arrays for coords */
  for(ni=0; ni<nnb; ni++)
  {
    free_array(Res[ni]);
    free_array(Cb[ni][1]);
    free_array(Cb[ni][0]);
    free_array(Ip[ni]);
  }
  free(Res);
  free(Cb);
  free(Ip);
  free_array(Cp[1]);
  free_array(Cp[0]);

end_set_ajsurf_forall_vars:
  TIMER_STOP;
}


/* copy from nbsurf[0] to ajsurf.  */
void copy_ajsurf_from_nbsurf0(tNode *node, int f, int nb_f,
                              int intrch, int rev1, int rev2)
{
  tDat *dat = node->dat;
  int dir = f/2;
  //int d1 = Dir1_norm(dir);
  //int d2 = Dir2_norm(dir);
  int vi;
  int i,j,k;
  int odir = nb_f/2;
  //int od1 = Dir1_norm(odir);
  //int od2 = Dir2_norm(odir);
  //int on1 = nb_n[od1];
  //int on2 = nb_n[od2];
  int oi,oj,ok;

  for(vi=0; vi<dat->nv; vi++)
  {
    tSurface *s = dat->s[f][vi];
    if(s)
    {
      int nzones = MeshVarSurfacezones(node->pat->mesh, vi);
      int zone, *n, *nn;

      /* unless it has mem already, get mem for ajsurf */
      if(s->allocd_ajsurf==0 || s->ajsurf==NULL)
      {
        /* get mem for ajsurf, needs to be freed later */
        s->ajsurf = alloc_array(s->mysurf->n);
        s->allocd_ajsurf = 1;
      }

      n  = s->ajsurf->n;
      nn = s->nbsurf[0]->n;

      /* now copy s->nbsurf[0] into s->ajsurf */
      for(zone=0; zone<nzones; zone++)
      {
        forplaneN(dir, i,j,k, n, zone)
        {
          int ind1 = Ind_n(i,j,k, s->ajsurf->n);
          int I1 = i1_norm(i,j,k, dir);
          int I2 = i2_norm(i,j,k, dir);
          int i2,j2,k2, ind2;

          switch(odir)
          {
          case 0:  oi = zone;  oj = I1;    ok = I2;    break;
          case 1:  oi = I1;    oj = zone;  ok = I2;    break;
          case 2:  oi = I1;    oj = I2;    ok = zone;  break;
          default: errorexit("odir must be 0,1,2");
          }

          i2 = oi;  j2 = oj;  k2 = ok;

          switch(odir)
          {
          case 0:
            if(!intrch)
            {
              if(rev1) { j2 = nn[1]-1 - oj; }
              if(rev2) { k2 = nn[2]-1 - ok; }
            }
            else
            {
              j2 = ok;  k2 = oj;
              if(rev1) { j2 = nn[1]-1 - ok; }
              if(rev2) { k2 = nn[2]-1 - oj; }
            }
            break;
          case 1:
            if(!intrch)
            {
              if(rev1) { i2 = nn[0]-1 - oi; }
              if(rev2) { k2 = nn[2]-1 - ok; }
            }
            else
            {
              i2 = ok;  k2 = oi;
              if(rev1) { i2 = nn[0]-1 - ok; }
              if(rev2) { k2 = nn[2]-1 - oi; }
            }
            break;
          case 2:
            if(!intrch)
            {
              if(rev1) { i2 = nn[0]-1 - oi; }
              if(rev2) { j2 = nn[1]-1 - oj; }
            }
            else
            {
              i2 = oj;  j2 = oi;
              if(rev1) { i2 = nn[0]-1 - oj; }
              if(rev2) { j2 = nn[1]-1 - oi; }
            }
            break;
          }
          //printf("   ajsurf: %d %d %d  nbsurf: %d %d %d\n", i,j,k, i2,j2,k2);
          ind2 = Ind_n(i2,j2,k2, nn);
          s->ajsurf->d[ind1] = s->nbsurf[0]->d[ind2];
        }
      } /* end zone loop */
    }
  }
  return;
}


/* set all ajsurf of a node */
void set_all_ajsurf(tNode *node)
{
  tDat *dat = node->dat;
  int face;

  if(!dat) return;

  for(face=0; face<6; face++)
    set_ajsurf_forall_vars(node, face);
}

/* get nbsurf for all nodes out of buffers and free the buffers */
void set_all_myln_ajsurf(tMesh *mesh)
{
  formylnodes(mesh)
  {
    tNode *node = MyLnode;
    set_all_ajsurf(node);
  }
}


/**********************************************************************/
/* some funcs to free and copy parts of the surfaces */
/**********************************************************************/

/* free all nbsurf in a tSurface, this can be useful to free up memory
   after ajsurf is constructed, and it does not touch ajsurf */
void free_nbsurf_only(tSurface *s)
{
  int i;

  if(!s) return;
  if(!s->nbsurf) return;

  /* free content of lists */
  /* free nbsurf[i] only if it is allocd  */
  for(i=0; i<s->nnbsurf; i++)
  {
    /* if it was allocated we free,
       unless someone else is still pointing to it */
    if(s->allocd_nbsurf[i])
    {
      /* if s->ajsurf is pointing to s->nbsurf[i] do not free the array */
      if(s->ajsurf == s->nbsurf[i])
      {
        if(s->allocd_ajsurf == 1)
          errorexit("if both s->ajsurf and s->nbsurf[i] were allocated,\n"
                    "they shouldn't point to the same array");
        s->allocd_ajsurf = 1;
      }
      else
        free_array(s->nbsurf[i]);
    }

    /* set it to NULL */
    s->nbsurf[i] = NULL;
  }

  /* make sure we know the s->nbsurf[i] are all gone */
  s->nnbsurf = 0;

  /* free lists */
  free(s->nbsurf);
  free(s->allocd_nbsurf);
  s->nbsurf = NULL;
  s->allocd_nbsurf = NULL;
}

/* free all nbsurf on face f of node */
void free_nbsurf_only_forall_vars(tNode *node, int f)
{
  tDat *dat = node->dat;
  int vi;

  if(!dat) return;

  /* free all nbsurf on face f */
  for(vi=0; vi<dat->nv; vi++)
    free_nbsurf_only(dat->s[f][vi]);
}

/* free all nbsurf surfaces on node */
void free_all_nbsurf_only(tNode *node)
{
  tDat *dat = node->dat;
  int vi,f;

  if(!dat) return;

  /* free all nbsurf */
  for(f=0; f<6; f++)
    for(vi=0; vi<dat->nv; vi++)
      free_nbsurf_only(dat->s[f][vi]);
}

/* free nbsurf on all nodes in the mesh */
void free_all_myln_nbsurf_only(tMesh *mesh)
{
  formylnodes(mesh)
  {
    tNode *node = MyLnode;
    free_all_nbsurf_only(node);
  }
}


/* free all ajsurf in a tSurface, this can be useful if n,pt_typ have
   changed and we want to recompute the ajsurf */
void free_ajsurf_only(tSurface *s)
{
  if(!s) return;
  if(!s->ajsurf) return;

  /* free ajsurf only if is allocd */
  if(s->allocd_ajsurf) free_array(s->ajsurf);

  /* make it NULL again, as it was before it was set */
  s->ajsurf = NULL;
}

/* free all ajsurf on face f of node */
void free_ajsurf_only_forall_vars(tNode *node, int f)
{
  tDat *dat = node->dat;
  int vi;

  if(!dat) return;

  /* free all ajsurf on face f */
  for(vi=0; vi<dat->nv; vi++)
    free_ajsurf_only(dat->s[f][vi]);
}

/* free all ajsurf surfaces on node */
void free_all_ajsurf_only(tNode *node)
{
  tDat *dat = node->dat;
  int vi,f;

  if(!dat) return;

  /* free all ajsurf */
  for(f=0; f<6; f++)
    for(vi=0; vi<dat->nv; vi++)
      free_ajsurf_only(dat->s[f][vi]);
}


/* Copy nbsurf data from s to sdest, this allocates sdest.
   It copies all nbsurf arrays, so that s could be freed after this. */
tSurface *surface_copy_with_nbsurf_only(tSurface *s)
{
  tSurface *sdest;
  int i;

  errorexit("use surface_point_to_same_nbsurf instead");

  if(!s) return NULL;
  if(!s->nbsurf) return NULL;

  /* make sdest and copy some info */
  sdest = alloc_empty_surface(s->nnbsurf);
  sdest->dat  = s->dat;
  sdest->face = s->face;
  sdest->vi   = s->vi;

  /* now copy all nbsurf arrays, and mark them as allocated */
  for(i=0; i<s->nnbsurf; i++)
  {
    sdest->nbsurf[i] = array_copy(s->nbsurf[i]);
    sdest->allocd_nbsurf[i] = 1;
  }
  return sdest;
}

/* copy all nbsurf from node_src to node_dest
   this allocates room for the surfaces but not their data. */
void surface_copy_all_nbsurf_only(tNode *node_src, tNode *node_dest)
{
  tDat *dat_src = node_src->dat;
  tDat *dat_dest = node_dest->dat;
  int vi,f;

  errorexit("use surface_copy_nbsurf_pointers instead");

  if(!dat_src || !dat_dest) return;

  /* loop over faces, vars and copy nbsurf of each surf */
  for(f=0; f<6; f++)
    for(vi=0; vi<dat_dest->nv; vi++)
      dat_dest->s[f][vi] = surface_copy_with_nbsurf_only(dat_src->s[f][vi]);
}

/* Copy mysurf,nbsurf,ajsurf pointers from s to sdest, this allocates sdest.
   mysurf,nbsurf,ajsurf will be marked as not allocated.
   If we free_surface the new surface, all the data that mysurf,nbsurf,ajsurf
   point to will not be freed. */
tSurface *surface_point_to_same_data(tSurface *s)
{
  tSurface *sdest;
  int i;

  if(!s) return NULL;
  if(!s->nbsurf) return NULL;

  /* make sdest and copy some info */
  sdest = alloc_empty_surface(s->nnbsurf);
  sdest->dat  = s->dat;
  sdest->face = s->face;
  sdest->vi   = s->vi;

  /* now copy all surface pointers, and mark them as not allocated */
  sdest->mysurf = s->mysurf;
  sdest->allocd_mysurf = 0;
  sdest->ajsurf = s->ajsurf;
  sdest->allocd_ajsurf = 0;
  for(i=0; i<s->nnbsurf; i++)
  {
    sdest->nbsurf[i] = s->nbsurf[i];
    sdest->allocd_nbsurf[i] = 0;
  }
  return sdest;
}

/* Copy mysurf,nbsurf,ajsurf pointers from node_src to node_dest,
   this allocates room for the surfaces but not their data. */
void surface_copy_all_pointers(tNode *node_src, tNode *node_dest)
{
  tDat *dat_src = node_src->dat;
  tDat *dat_dest = node_dest->dat;
  int vi,f;

  if(!dat_src || !dat_dest) return;

  /* loop over faces, vars and copy nbsurf of each surf */
  for(f=0; f<6; f++)
    for(vi=0; vi<dat_dest->nv; vi++)
      dat_dest->s[f][vi] = surface_point_to_same_data(dat_src->s[f][vi]);
}


/* Copy nbsurf pointers from s to sdest.
   nbsurf will be marked as not allocated in sdest.
   If we free_surface sdest, all the data that nbsurf points to will not
   be freed. */
tSurface *surface_point_to_same_nbsurf(tSurface *s, tSurface *sdest)
{
  int i;

  if(!s) return NULL;
  if(!s->nbsurf) return NULL;

  if(s->nnbsurf != sdest->nnbsurf)
    errorexit("the two nbsurf lists differ in size");

  /* now copy all nbsurf pointers, and mark them as not allocated */
  for(i=0; i<s->nnbsurf; i++)
  {
    sdest->nbsurf[i] = s->nbsurf[i];
    sdest->allocd_nbsurf[i] = 0;
  }
  return sdest;
}

/* Copy nbsurf pointers from node_src to node_dest */
void surface_copy_nbsurf_pointers(tNode *node_src, tNode *node_dest)
{
  tDat *dat_src = node_src->dat;
  tDat *dat_dest = node_dest->dat;
  int vi,f;

  if(!dat_src || !dat_dest) return;

  /* loop over faces, vars and copy nbsurf of each surf */
  for(f=0; f<6; f++)
    for(vi=0; vi<dat_dest->nv; vi++)
      surface_point_to_same_nbsurf(dat_src->s[f][vi], dat_dest->s[f][vi]);
}


/*************************************************************************/
/*************************************************************************/

/* NOTE:
   If we do non-uniform timestepping, each node needs to get surface info
   only from its neighbors. I.e. only the neighbors need to send data
   to the node (equiv. to MPI_Irecv without MPI_Isend on the node).
   Each node also only needs info for a VarList not all vars.
   So we should make functions that:
   1.: -create mysurf info for a VarList on one face of one node
   2.: -put these mysurf in a MPI window for one-sided communication
        (maybe use segmented array)
       -this needs to be done for all in VarList and on one face
       -we then use this func in a loop over all neighbor nodefaces
   3.: -create funcs that read from MPI windows (with mysurf) of all
        neighbor nodefaces to set nbsurf
*/

/*************************************************************************/
/*************************************************************************/

/* we need the functions below, but they are UNFINISHED!!! */


/* init surfaces for a VarList */
void init_all_vl_surfaces(tMesh *mesh, tVarList *vl)
{
  //FIXME: we should loop over vl only, but for now we just request
  //       surfaces for all vars and do this:
  init_all_myln_surfaces(mesh);
}

/* set mysurf for one node and a VarList */
void set_all_vl_mysurf(tNode *node, tVarList *vl)
{
  //FIXME: we should loop over vl only, but for now we just set
  //       surfaces for all vars and do this:
  set_all_mysurf(node);
}

/* request surface exchange for one node and a VarList */
void request_all_vl_surfaces(tNode *node, tVarList *vl)
{
  //FIXME: we should loop over vl only, but for now we just request
  //       surfaces for all vars and do this:
  request_all_surfaces_exchange(node);
  /* ^-We should really adapt this to using MPI Windows */

  // NOTE: request_all_surfaces_exchange has a Waitall in
  // free_dat_reqs_after_Waitall_com_send that leads to a deadlock when it
  // is called repeatedly from the same node, as I would like to do in RK4.
}

/* get surfaces on one node for a varlist */
void get_all_vl_surfaces(tNode *node, tVarList *vl)
{
  //FIXME: we should loop over vl only, but for now we just get
  //       surfaces for all vars and do this:
  get_all_surfaces(node);
  free_dat_reqs_after_Waitall_com_send(node); //not needed if we use MPI windows
}

/* free all surfaces on one node for a varlist */
void free_all_vl_surfaces(tNode *node, tVarList *vl)
{
  //FIXME: not needed if we use MPI windows
  free_dat_reqs_after_Waitall_com_send(node);

  //FIXME: we could free only the surfaces associated with vl, but for now:
  free_all_surfaces(node);
}
