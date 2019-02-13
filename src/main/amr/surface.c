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
  s->recv_req = calloc(nnb, sizeof(s->recv_req[0]));
  s->send_req = calloc(nnb, sizeof(s->send_req[0]));
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
  free(s->recv_req);
  free(s->send_req);
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

/* put the nbsurf from neighbor with index ni in s->face */
void get_nbsurf(tSurface *s, int ni, int zones)
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
    s_req = &(s->send_req[ni]);
    r_req = &(s->recv_req[ni]);
    nMPI_Isend_Irecv_double(s->mysurf->a, s->mysurf->N,
                            s->nbsurf[ni]->a, s->nbsurf[ni]->N,
                            nb_rank, s_tag, r_tag, s_req, r_req);
  }
}

/* put nbsurf from all neighbors in s->face */
void get_all_nbsurf(tSurface *s)
{
  tNode *node = s->dat->node;
  int zones = MeshVarSurfacezones(node->pat->mesh, s->vi);
  int my_f = s->face;
  int nfnb = node->nfnb[my_f];
  int ni;
  for(ni=0; ni<nfnb; ni++) get_nbsurf(s, ni, zones);
}


/* get nbsurf from all faces and variables for this node */
void get_all_surfaces(tNode *node)
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
      get_all_nbsurf(s);
    }
  }
}
