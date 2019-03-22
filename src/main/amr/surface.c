/* surface.c */
/* Wolfgang Tichy, 2/2019 */

#include "nmesh.h"
#include "amr.h"

#define PR 0

/* functions to exchange surface data */

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
  int myid;
  formylnodes(mesh, myid)
  {
    tNode *node = MyNode(mesh, myid);
    free_all_surfaces(node);
  }
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
  alloc_mysurf = 1;
  if(n[dir] == zones) alloc_mysurf = 0;
  else                n[dir] = zones;

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
  int myid;
  formylnodes(mesh, myid)
  {
    tNode *node = MyNode(mesh, myid);
    init_all_surfaces(node);
  }
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
  int p = (node->n[dir] - 1) * (f%2); /* plane of surface */
  int vi = s->vi;

  /* do nothing if this surface is NULL */
  if(!s) return;

  if(s->allocd_mysurf)
    copy_array_plane(dat->v[vi], dir, p, s->mysurf, 0);
  else
    s->mysurf = dat->v[vi];
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
  int myid;
  formylnodes(mesh, myid)
  {
    tNode *node = MyNode(mesh, myid);
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
    printnode(node);
    printbfaces_on_f(node->pat, face);
    printf("neighbor node:\n");
    printnode(nb);
    printbfaces_on_f(nb->pat, nb_f);
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
    int rq, nb_rank, s_tag, r_tag;
    nMPI_Comm s_comm, r_comm;
    tCom *com = dat->com[face];
    int nb_n[3], nb_N;
    int my_n[3], my_N;
    double *sbuf, *rbuf; /* buffers for MPI */
    int zones, nvars, cnt;

    /* count number of vars that have surfaces to exchange and set my_N */
    find_nvars_vind_n_nbn(node,face, nb,nb_f, &nvars, &vi, my_n, nb_n);
    nb_N = nb_n[0] * nb_n[1] * nb_n[2];
    my_N = my_n[0] * my_n[1] * my_n[2];

    /* use MPI to recv nb->dat->s[nb_f][vi]->mysurf in s->nbsurf[ni],
       and also send s->mysurf to nb->dat->s[nb_f][vi]->nbsurf[nb_ni] */
    nb_rank = nb->datrank;
    r_tag = (node->nid)*6 + face;
    s_tag = (nb->nid)*6 + nb_f;
    r_comm = nb->comm;
    s_comm = node->comm;

    /* alloc send and recv buffers */
    sbuf = calloc(nvars * my_N, sizeof(double));
    rbuf = calloc(nvars * nb_N, sizeof(double));
    /* save buffers in com */
    rq = append_buffers_to_com(com, sbuf, nvars*my_N, rbuf,nvars*nb_N);
    //NOTE: it may be good to use a long segmented array as rbuf
    //      with dat->s[face][vi]->nbsurf[ni] pointing to the segments

    /* fill send buffer */
    for(cnt=0, vi=0; vi<node->dat->nv; vi++)
    {
      tSurface *s = dat->s[face][vi];

      zones = MeshVarSurfacezones(node->pat->mesh, vi);
      /* do nothing if var has no zones to exchange */
      if(s && zones && dat->v[vi])
      {
        /* allocate surface to later recv neighbor data */
        if(!s->nbsurf[ni])
        {
          s->nbsurf[ni] = alloc_array(nb_n);
          s->allocd_nbsurf[ni] = 1; //flag that we allocd
        }

        /* save MPI request number in the array */
        s->nbsurf[ni]->d[0] = rq;

        /* fill buffer for MPI exchange: sbuf[] = s->mysurf->d[] */
        memcpy(sbuf+cnt, s->mysurf->d, my_N * sizeof(sbuf[0]));
        cnt += my_N;
      }
    }
    /* now call MPI */
    nMPI_Isend_Irecv_double_com(com, rq, nb_rank, s_tag,r_tag, s_comm,r_comm);
  }
}


/* put nbsurf from all faces and variables for this node in buffers */
void request_all_surfaces_exchange(tNode *node)
{
  int face, ni;

  /* do nothing if this node is on other proc */
  if(!node->dat) return;

  /* free req, send/recv arrays */
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
}

/* request surface exchanges on all my nodes in the mesh
   Note: We need to call this! If we call request_all_surfaces_exchange(n1)
   for only node n1, MPI deadlocks because the other nodes are not sending
   to n1 or receiving from n1 */
void request_all_myln_surfaces_exchange(tMesh *mesh)
{
  int myid;

  TIMER_START;

  /* If we want threads in this loop, we need MPI_Init_thread with
     MPI_THREAD_MULTIPLE, instead of just MPI_Init in main. */
  formylnodes_noomp(mesh, myid)
  {
    tNode *node = MyNode(mesh, myid);
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
  int nb_n[3], nb_N;
  int my_n[3];
  int zones, nvars, rq, cnt;
  double *rbuf;

  /* do nothing if this node is on other proc */
  if(!dat) return;

  /* do nothing if com is empty */
  if(com->n_rq == 0) return;

  /* find face nb_f of nb that faces me */
  found = locate_facenb_in_fnbs(nb, node, &nb_f, &nb_ni);
  if(!found) errorexit("couldn't find nb face!!!");

  /* is nb is local we have already exchanged info  */
  if(nb->dat) return;

  /* nb is on other process, we have used MPI to exchange data */

  /* count number of vars that exchanged surfaces and set my_N */
  find_nvars_vind_n_nbn(node,face, nb,nb_f, &nvars, &vi, my_n, nb_n);
  nb_N = nb_n[0] * nb_n[1] * nb_n[2];
  /* do nothing if there are no vars that exchanged surfaces */
  if(!nvars) return;

  //PRF;printf(": nvars=%d\n", nvars);

  /* get MPI request number */
  rq = dat->s[face][vi]->nbsurf[ni]->d[0];
  /* find our recv buffer */
  rbuf = get_com_recv_buf(com, rq);

  /* wait for buffer */
  nMPI_Wait_com_recv(com, rq);

  /* get data out of recv buffer */
  for(cnt=0, vi=0; vi<node->dat->nv; vi++)
  {
    tSurface *s = dat->s[face][vi];

    zones = MeshVarSurfacezones(node->pat->mesh, vi);
    if(!zones) continue; /* do nothing if var has no zones to exchange */
    /* do nothing if var has no zones to exchange */
    if(zones && dat->v[vi])
    {
      /* get neighbor data from buffer */
      memcpy(s->nbsurf[ni]->d, rbuf+cnt, nb_N * sizeof(rbuf[0]));
      cnt += nb_N;
    }
  }
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
    nMPI_Waitall_com_recv(dat->com[face]);

    /* wait until all has been sent, then free all buffers for this face */
    nMPI_Waitall_com_send(dat->com[face]);
    realloc_dat_reqs(node->dat, 0, face); /* free req and send arrays */
  }
}

/* get nbsurf from all faces and variables for this node out of buffers */
void get_all_surfaces(tNode *node)
{
  int face, ni;
  tDat *dat = node->dat;

  if(!dat) return;

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
}

/* get nbsurf for all nodes out of buffers and free the buffers */
void get_all_myln_surfaces(tMesh *mesh)
{
  int myid;

  TIMER_START;

  /* If we want threads in this loop, we need MPI_Init_thread with
     MPI_THREAD_MULTIPLE, instead of just MPI_Init in main. */
  formylnodes_noomp(mesh, myid)
  {
    tNode *node = MyNode(mesh, myid);
    get_all_surfaces(node);
  }

  /* postpone Waitall until we have finished all nodefaces. This could have
     been already called in get_all_surfaces to free mem earlier.*/
  formylnodes_noomp(mesh, myid)
  {
    tNode *node = MyNode(mesh, myid);
    free_dat_reqs_after_Waitall_com_send(node);
  }

  TIMER_STOP;
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

/* do two nodes have same point number n, orthogonal to dir? */
int same_n_normal_to_dir(tNode *node1, tNode *node2, int dir)
{
  int d, samen = 1;

  for(d=0; d<3; d++)
    if(d!=dir)
    {
      if(node1->n[d] != node2->n[d]) { samen=0; break; }
    }
  return samen;
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
      if(nb->l == node->l)
      {
        int same_n = same_n_normal_to_dir(nb, node, dir);
        /* if number of points is the same we can copy or just point
           ajsurf to nbsurf[0] */
        if(same_n)
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
      found = locate_facenb_in_fnbs(nb, node, &nb_f, &nb_ni);
      if(!found) errorexit("couldn't find nb face!!!");

      /* if nb has only node as neighbor we may not need to interpolate */
      if(1 && nb->nfnb[nb_f] == 1)
      {
        tMesh *mesh = node->pat->mesh;
        int ix = Ind("x");
        double *px[] = { Vard(node,ix), Vard(node,ix+1), Vard(node,ix+2) };
        int *n = node->n;
        int pl = (n[dir]-1)*(f%2);
        int d1 = Dir1_norm(dir);
        int d2 = Dir2_norm(dir);
        int n1 = n[d1];
        int n2 = n[d2];
        int i,j,k, ind;
        int *nb_n = nb->n;
        int inb[3][3];
        double x0[3], dist[3], mx0;
        int odir = nb_f/2;
        int opl = (nb_n[odir]-1)*(nb_f%2);
        int od1 = Dir1_norm(odir);
        int od2 = Dir2_norm(odir);
        int on1 = nb_n[od1];
        int on2 = nb_n[od2];

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
            if(inb[1][od1] == 0 && inb[1][od2] == on2-1)
            {
              if(inb[2][od1] == on1-1 && inb[2][od2] == 0)
              {
                /* copy with two axis interchanged */
                if(PR) printf(" copy from nbsurf[0] with two axis interchanged\n");
                copy_ajsurf_from_nbsurf0(node,f,nb_f, 1,0,0);
                goto end_set_ajsurf_forall_vars;
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
                /* copy with axis1 reversed */
                if(PR) printf(" copy from nbsurf[0] with axis1 reversed\n");
                copy_ajsurf_from_nbsurf0(node,f,nb_f, 0,1,0);
                goto end_set_ajsurf_forall_vars;
              }
            }
            if(inb[1][od1] == on1-1 && inb[1][od2] == on2-1)
            {
              if(inb[2][od1] == 0 && inb[2][od2] == 0)
              {
                /* copy with two axis interchanged and axis1 reversed */
                if(PR) printf(" copy from nbsurf[0] with two axis interchanged and axis1 reversed\n");
                copy_ajsurf_from_nbsurf0(node,f,nb_f, 1,1,0);
                goto end_set_ajsurf_forall_vars;
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
                /* copy with axis2 reversed */
                if(PR) printf(" copy from nbsurf[0] with axis2 reversed\n");
                copy_ajsurf_from_nbsurf0(node,f,nb_f, 0,0,1);
                goto end_set_ajsurf_forall_vars;
              }
            }
            if(inb[1][od1] == 0 && inb[1][od2] == 0)
            {
              if(inb[2][od1] == on1-1 && inb[2][od2] == on2-1)
              {
                /* copy with two axis interchanged and axis2 reversed */
                if(PR) printf(" copy from nbsurf[0] with two axis interchanged and axis2 reversed\n");
                copy_ajsurf_from_nbsurf0(node,f,nb_f, 1,0,1);
                goto end_set_ajsurf_forall_vars;
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
                /* copy with both axis reversed */
                if(PR) printf(" copy from nbsurf[0] with both axis reversed\n");
                copy_ajsurf_from_nbsurf0(node,f,nb_f, 0,1,1);
                goto end_set_ajsurf_forall_vars;
              }
            }
            if(inb[1][od1] == on1-1 && inb[1][od2] == 0)
            {
              if(inb[2][od1] == 0 && inb[2][od2] == on2-1)
              {
                /* copy with two axis interchanged and both reversed */
                if(PR) printf(" copy from nbsurf[0] with both axis interchanged and both reversed\n");
                copy_ajsurf_from_nbsurf0(node,f,nb_f, 1,1,1);
                goto end_set_ajsurf_forall_vars;
              }
            }
          }
        }
        /* nothing matches if we get here and we need to interpolate */
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
      if(0 && f==1 && node->nid==7)
      {
        printf("ni=%d:\n", ni);
        printnode(node);
        printbfaces_on_f(node->pat, f);
        printnode(nb);
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
        nb = node->fnb[f][ni];
        found = locate_facenb_in_fnbs(nb, node, &nb_f, &nb_ni);
        if(!found) errorexit("couldn't find nb face!!!");
        nb_dir = nb_f/2;
        Lagrange_interpolate2d_toIpoints(nb, s->nbsurf[ni], nb_dir,0,
                                         Cb[ni],Ip[ni], Res[ni]);

if(0 && node->nid==17 && nb->nid==64 && vi==35)
{
tMesh *mesh = node->pat->mesh;
int ll;
printnode(node);
printnode(nb);
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
      int *n, *nn;

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
      forplaneN(dir, i,j,k, n, 0)
      {
        int ind1 = Ind_n(i,j,k, s->ajsurf->n);
        int I1 = i1_norm(i,j,k, dir);
        int I2 = i2_norm(i,j,k, dir);
        int i2,j2,k2, ind2;

        switch(odir)
        {
        case 0:  oi = 0;   oj = I1;  ok = I2;  break;
        case 1:  oi = I1;  oj = 0;   ok = I2;  break;
        case 2:  oi = I1;  oj = I2;  ok = 0;   break;
        default: errorexit("dir must be 0,1,2");
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
  int myid;
  formylnodes(mesh, myid)
  {
    tNode *node = MyNode(mesh, myid);
    set_all_ajsurf(node);
  }
}


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
  int myid;
  formylnodes(mesh, myid)
  {
    tNode *node = MyNode(mesh, myid);
    free_all_nbsurf_only(node);
  }
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
