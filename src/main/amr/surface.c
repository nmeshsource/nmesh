/* surface.c */
/* Wolfgang Tichy, 2/2019 */

#include "nmesh.h"
#include "amr.h"


/* functions to exchange surface data */

/**********************************************************************/
/* allocate and fill surfaces for vars that need it */
/**********************************************************************/
/* empty surface that we need to fill in */
tSurface *alloc_empty_surface(int nnb)
{
  tSurface *s = calloc(1, sizeof(*s));
  s->nbsurf = calloc(nnb, sizeof(s->nbsurf[0]));
  return s;
}

/* free all we need to in a surface */
void free_surface(tSurface *s)
{
  tNode *node;
  tDat *dat;
  int f, i;

  if(!s) return;
  dat = s->dat;
  f = s->face;
  node = dat->node;

  /* free content of lists */
  /* free mysurf only it has allocd */
  if(s->allocd_mysurf) free_array(s->mysurf);
  /* free nbsurf[i] only if it is not on this proc  */
  for(i=0; i<node->nfnb[f]; i++)
    if(!node->fnb[f][i]->dat) free_array(s->nbsurf[i]);

  /* free lists */
  free(s->nbsurf);
}


/* initialize a surface for var vi at face with nnb neighbors */
tSurface *init_surface(tNode *node, int vi, int face)
{
  int dir = face/2;
  int zones;
  tDat *dat;
  int i;
  tSurface *s;
  int n[3];
  int alloc_mysurf;

  /* do nothing if no data on this node */
  if(!node->dat) return NULL;
  dat = node->dat;

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
      dat->s[face][vi] = init_surface(node, vi, face);
      if(dat->s[face][vi]) cnt++;
    }
  }
  return cnt;
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
  if(s->allocd_mysurf)
    copy_array_plane(dat->v[vi], dir, p, s->mysurf, 0);
  else
    s->mysurf = dat->v[vi];
}


/* count number of vars that have surfaces to exchanged and set myN,
   input: node,my_f, nb,nb_f   output: nvars, vind, my_n, nb_n */
void find_nvars_vind_n_nbn(tNode *node, int my_f, tNode *nb, int nb_f,
                           int *nvars, int *vind, int my_n[3], int nb_n[3])
{
  int todo=1;
  int vi, i;
  int my_dir = my_f/2;
  int nb_dir = nb_f/2;

  /* count number of vars that have surfaces to exchanged and set myN */
  for(nvars=0, vi=0; vi<node->dat->nv; vi++)
  {
    int zones = MeshVarSurfacezones(node->pat->mesh, vi);
    if(zones)
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
  if(!found) errorexit("couldn't find nb face!!!");

  /* is nb local? */
  if(nb->dat)
  {
    /* nb is local so just point s->nbsurf[ni] to its data */
    for(vi=0; vi<node->dat->nv; vi++)
    {
      tSurface *s = dat->s[face][vi];
      tArray *nb_mysurf = nb->dat->s[nb_f][vi]->mysurf;
      s->nbsurf[ni] = nb_mysurf;
    }
  }
  else
  {
    /* nb is on other process so use MPI to exchange data */
    int n_rq, nb_rank, s_tag, r_tag;
    nMPI_Req *s_req, *r_req;
    int nb_n[3], nb_N;
    int my_n[3], my_N;
    double *sbuf, *rbuf; /* buffers for MPI */
    int zones, nvars, cnt;

    /* count number of vars that have surfaces to exchanged and set my_N */
    find_nvars_vind_n_nbn(node,face, nb,nb_f, &nvars, &vi, my_n, nb_n);
    nb_N = nb_n[0] * nb_n[1] * nb_n[2];
    my_N = my_n[0] * my_n[1] * my_n[2];

    /* use MPI to recv nb->dat->s[nb_f][vi]->mysurf in s->nbsurf[ni],
       and also send s->mysurf to nb->dat->s[nb_f][vi]->nbsurf[nb_ni] */
    nb_rank = nb->datrank;
    r_tag = ((node->nid)*256 + ni)*6 + face;
    s_tag = ((nb->nid)*256 + nb_ni)*6 + nb_f;
//FIXME: we need better tags!!!!

    n_rq = dat->n_rq[face]; /* number of MPI requests so far */
    /* make room for one more request */
    realloc_dat_reqs(dat, n_rq + 1, face);
    s_req = &(dat->send_rq[face][n_rq]);
    r_req = &(dat->recv_rq[face][n_rq]);
    /* alloc send and recv buffers */
    sbuf = calloc(nvars * my_N, sizeof(double));
    rbuf = calloc(nvars * nb_N, sizeof(double));
    dat->send_buf[face][n_rq] = sbuf; /* save buffers */
    dat->recv_buf[face][n_rq] = rbuf;
    //NOTE: it may be good to use a long segmented array as rbuf
    //      with dat->s[face][vi]->nbsurf[ni] pointing to the segments

    /* fill send buffer */
    for(cnt=0, vi=0; vi<node->dat->nv; vi++)
    {
      tSurface *s = dat->s[face][vi];

      zones = MeshVarSurfacezones(node->pat->mesh, vi);
      if(!zones) continue; /* do nothing if var has no zones to exchange */

      /* allocate surface to later recv neighbor data */
      s->nbsurf[ni] = alloc_array(nb_n);
      /* save MPI request number in the array */
      s->nbsurf[ni]->a[0] = n_rq;

      /* fill buffer for MPI exchange: sbuf[] = s->mysurf->a[] */
      memcpy(sbuf+cnt, s->mysurf->a, my_N);
      cnt += my_N;
    }
    /* no call MPI */
    nMPI_Isend_Irecv_double(sbuf, nvars*my_N,  rbuf, nvars*nb_N,
                            nb_rank, s_tag, r_tag, s_req, r_req);
  }
}


/* put nbsurf from all faces and variables for this node in buffers */
void request_all_surfaces_exchange(tNode *node)
{
  int face, ni;
  int ns = init_all_surfaces(node);

  /* do nothing if there are no surfaces */
  if(!ns) return;

  /* do nothing if this node is on other proc */
  if(!node->dat) return;

  for(face=0; face<6; face++)
  {
    realloc_dat_reqs(node->dat, 0, face); /* free req and send arrays */
    //nb_nid0 = node->fnb[face][0];
    for(ni=0; ni<node->nfnb[face]; ni++)
    {
      request_surfaces_exchange_for_all_vars(node, face, ni);
    }
  }
}


//request_surface_exchange_for_all_vars

/* get all surfaces from neighbor with index ni at face */
void get_surfaces_for_all_vars(tNode *node, int face, int ni)
{
  tNode *nb = node->fnb[face][ni];
  tDat *dat = node->dat;
  int nb_f, nb_ni;
  int found, vi;
  int nb_n[3], nb_N;
  int my_n[3];
  int zones, nvars, n_rq, cnt;
  nMPI_Stat stat;
  double *rbuf;

  /* do nothing if this node is on other proc */
  if(!dat) return;

  /* find face nb_f of nb that faces me */
  found = locate_facenb_in_fnbs(nb, node, &nb_f, &nb_ni);
  if(!found) errorexit("couldn't find nb face!!!");

  /* is nb is local we have already exchanged info  */
  if(nb->dat) return;

  /* nb is on other process, we have used MPI to exchange data */

  /* count number of vars that have surfaces to exchanged and set my_N */
  find_nvars_vind_n_nbn(node,face, nb,nb_f, &nvars, &vi, my_n, nb_n);
  nb_N = nb_n[0] * nb_n[1] * nb_n[2];

  /* get MPI request number */
  n_rq = dat->s[face][vi]->nbsurf[ni]->a[0];
  /* find our recv buffer */
  rbuf = dat->recv_buf[face][n_rq];
  /* wait for buffer */
  nMPI_Waitall(1, &(dat->recv_rq[face][n_rq]), &stat);

  /* get data out of recev buffer */
  for(cnt=0, vi=0; vi<node->dat->nv; vi++)
  {
    tSurface *s = dat->s[face][vi];

    zones = MeshVarSurfacezones(node->pat->mesh, vi);
    if(!zones) continue; /* do nothing if var has no zones to exchange */

    /* get neighbor data from buffer */
    memcpy(s->nbsurf[ni]->a, rbuf+cnt, nb_N);

    cnt += nb_N;
  }
}


/* get nbsurf from all faces and variables for this node out of buffers */
void get_all_surfaces(tNode *node)
{
  int face, ni;

  for(face=0; face<6; face++)
  {
    for(ni=0; ni<node->nfnb[face]; ni++)
    {
      get_surfaces_for_all_vars(node, face, ni);
    }
    realloc_dat_reqs(node->dat, 0, face); /* free req and send arrays */
  }
}





/*******************************************************************/
/* Everything below this line is untested and may not work
   Probably it should be removed !!!!! */
/*******************************************************************/

/* put the nbsurf from neighbor with index ni in s->face */
void get_nbsurf__old(tSurface *s, int ni, int zones)
{
  int vi = s->vi;
  int my_f = s->face;
  int nb_f, nb_ni;
  tNode *node = s->dat->node;
  tNode *nb = node->fnb[my_f][ni];

  /* find face nb_f of nb that faces me */
  locate_facenb_in_fnbs(nb, node, &nb_f, &nb_ni);

  /* is nb local? */
  if(nb->dat)
  {
    /* nb is local so just point s->nbsurf[ni] to its data */
    tArray *nb_mysurf = nb->dat->s[nb_f][vi]->mysurf;
    s->nbsurf[ni] = nb_mysurf;
  }
  else
  {
    /* nb is on other process so use MPI to exchange data */
    int nb_rank, s_tag, r_tag;
    nMPI_Req *s_req, *r_req;
    int nb_dir = nb_f/2;
    int nb_n[3];
    int i;

    /* set nb_n */
    for(i=0; i<3; i++) nb_n[i] = nb->n[i];
    nb_n[nb_dir] = zones;

    /* allocate surface to recv neighbor data */
    s->nbsurf[ni] = alloc_array(nb_n);

    /* use MPI to recv nb->dat->s[nb_f][vi]->mysurf in s->nbsurf[ni],
       and also send s->mysurf to nb->dat->s[nb_f][vi]->nbsurf[nb_ni] */
    nb_rank = nb->datrank;
    s_tag = node->nid;
    r_tag = nb->nid;
    //s_req = &(s->send_req[ni]);
    //r_req = &(s->recv_req[ni]);
    nMPI_Isend_Irecv_double(s->mysurf->a, s->mysurf->N,
                            s->nbsurf[ni]->a, s->nbsurf[ni]->N,
                            nb_rank, s_tag, r_tag, s_req, r_req);
  }
}

/* put nbsurf from all neighbors in s->face */
void get_all_nbsurf__old(tSurface *s)
{
  tNode *node = s->dat->node;
  int zones = MeshVarSurfacezones(node->pat->mesh, s->vi);
  int my_f = s->face;
  int nfnb = node->nfnb[my_f];
  int ni;
  for(ni=0; ni<nfnb; ni++) get_nbsurf__old(s, ni, zones);
}

/* get nbsurf from all faces and variables for this node */
void get_all_surfaces__old(tNode *node)
{
  tDat *dat = node->dat;
  tSurface *s;
  int face, vi;
  int ns = init_all_surfaces(node);

  if(!ns) return;

  for(face=0; face<6; face++)
  {
    for(vi=0; vi<dat->nv; vi++)
    {
      s = dat->s[face][vi];
      if(!s) continue;
      get_all_nbsurf__old(s);
    }
  }
}
