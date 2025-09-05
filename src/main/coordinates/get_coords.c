/* get_coords.c */
/* Wolfgang Tichy, 2/2019 */

#include "nmesh.h"
#include "coordinates.h"


/* globals for coordinates */
extern tcoordinates coordinates[1];


/* get Xb from i,j,k */
void XbYbZb_of_ijk(tNode *node, int i, int j, int k, double Xb[3])
{
  tArray *A[3];
  int m[] = { i,j,k };
  int dir;

  node_Xb3(node, A);
  for(dir=0; dir<3; dir++) Xb[dir] = A[dir]->d[m[dir]];
}
/* get Xb from index ind */
void XbYbZb_of_ind(tNode *node, int ind, double Xb[3])
{
  int *n = node->n;
  int k = kOfInd_n(ind, n);
  int j = jOfInd_n_k(ind, n, k);
  int i = iOfInd_n_jk(ind, n, j,k);

  XbYbZb_of_ijk(node, i,j,k, Xb);
}

/* return i0 closest to Xb0 in direction dir */
int nearest_i0_of_Xb_indir(tNode *node, int dir, const double Xb0)
{
  int *n = node->n;
  double *Xb = node_Xb(node,dir)->d;
  int i,j,k;
  double a,b;
  int i0;

  i = 0;
  k = n[dir]-1;
  /* set a, b */
  Xb_of_X_indir(node, dir, &a, node->bbox[2*dir]);
  Xb_of_X_indir(node, dir, &b, node->bbox[2*dir + 1]);
  if(dgreater(Xb0, b)) return -k-1;
  if(dless(Xb0, a))    return -1;

  while(k-i>1)
  {
    j=(i+k)/2;

    a = Xb[i] - Xb0;
    b = Xb[j] - Xb0;
    if(a*b<=0.0) k=j;
    else         i=j;
  }
  if( fabs(Xb[i] - Xb0) < fabs(Xb[k] - Xb0) )
    i0 = i;
  else
    i0 = k;

  return i0;
}

/* find i,j,k closest to Xb0 */
void nearest_ijk_of_XbYbZb(tNode *node, int ijk[3], const double Xb0[3])
{
  int dir;
  for(dir=0; dir<3; dir++)
    ijk[dir] = nearest_i0_of_Xb_indir(node, dir, Xb0[dir]);
}

/* find i,j,k closest to X0 */
void nearest_ijk_of_XYZ(tNode *node, int ijk[3], const double X0[3])
{
  double Xb[3];
  XbYbZb_of_XYZ(node, Xb, X0);
  nearest_ijk_of_XbYbZb(node, ijk, Xb);
}

/* find i,j,k closest to X0+epsilon */
void nearest_ijk_of_XYZplus(tNode *node, int ijk[3], const double X0[3])
{
  double X0_plus[3];
  int dir;

  /* create a point X0_plus, that is moved in the positive dir (1,1,1) */
  for(dir=0; dir<3; dir++)
  {
    double X = X0[dir];
    X0_plus[dir] = X + (1. + fabs(X))*10.*dequaleps;
  }
  /* get ijk of X0_plus */
  nearest_ijk_of_XYZ(node, ijk, X0_plus);
}

/* find i,j,k closest to X0, but discard ijk[d]=0 in upper nodes */
void nearest_lowernode_ijk_of_XYZ(tNode *node, int ijk[3], const double X0[3])
{
  int nd_ijk = elm_get_ijk(node);
  tNode nb[1];
  int ijk2[3];

  nearest_ijk_of_XYZ(node, ijk, X0);

  /* invalidate ijk if node->ijk=nd_ijk is in the upper part of 0,1 or 2
     dirs, and neighbor in lower part has an ijk as well, i.e. if X0 is in
     the surface between both upper and lower nodes */
  if(ijk[0]==0 && nd_ijk%2)
  {
    //nb = node->parent->child[nd_ijk-1];
    amr_set_sibling_elm0(node, nd_ijk-1, nb);
    nearest_ijk_of_XYZ(nb, ijk2, X0);
    if(ijk2[0]==nb->n[0]-1)
      ijk[0] = -1;
  }
  if(ijk[1]==0 && (nd_ijk/2)%2)
  {
    //nb = node->parent->child[nd_ijk-2];
    amr_set_sibling_elm0(node, nd_ijk-2, nb);
    nearest_ijk_of_XYZ(nb, ijk2, X0);
    if(ijk2[1]==nb->n[1]-1)
      ijk[1] = -1;
  }
  if(ijk[2]==0 && nd_ijk>=4)
  {
    //nb = node->parent->child[nd_ijk-4];
    amr_set_sibling_elm0(node, nd_ijk-4, nb);
    nearest_ijk_of_XYZ(nb, ijk2, X0);
    if(ijk2[2]==nb->n[2]-1)
      ijk[2] = -1;
  }
}

/* find i,j,k closest to x0, return distance */
double nearest_ijk_of_xyz(tNode *node, int ijk[3], const double x0[3])
{
  tMesh *mesh = node->pat->mesh;
  int ix = Ind("x");
  double *px[] = { Vard(node,ix), Vard(node,ix+1), Vard(node,ix+2) };
  int *n = node->n;
  int i,j,k, dir, ind;
  double d1, d2, dist2 = DBL_MAX;

  /* if this node has no dat we need to compute x for each point from Xb */
  if(!px[0]) errorexit("we need x,y,z in this node");

  forijk(i,j,k, n)
  {
    ind = Ind_n(i,j,k, n);
    for(d2=0., dir=0; dir<3; dir++)
    {
      d1 = px[dir][ind] - x0[dir];
      d2 += d1*d1;
    }
    if(d2<dist2)
    {
      ijk[0] = i;
      ijk[1] = j;
      ijk[2] = k;
      dist2 = d2;
    }
  }
  return sqrt(dist2);
}

/* find i,j,k closest to x0 in plane pl normal to N, return distance */
double nearest_ijk_of_xyz_inplaneN(tNode *node, int N, int pl,
                                   int ijk[3], const double x0[3])
{
  tMesh *mesh = node->pat->mesh;
  int ix = Ind("x");
  double *px[] = { Vard(node,ix), Vard(node,ix+1), Vard(node,ix+2) };
  tArray *Cp[2];
  int *n = node->n;
  int n2[3];
  int dir1 = Dir1_norm(N);
  int dir2 = Dir2_norm(N);
  double Xpl=0, X[3], x[3];
  int dir, i,j,k, ind;
  double d1, d2, dist2 = DBL_MAX;

  /* if this node has no dat we need to compute x for each point from Xb */
  if(!px[0])
  {
    double Xb = node_Xb(node,N)->d[pl]; /* Xb coord in plane pl */

    /* convert Xb to X coords */
    X_of_Xb_indir(node, N, Xb, &Xpl);
    X[N]    = Xpl;

    /* make 2d n2 */
    for(dir=0; dir<3; dir++) n2[dir] = node->n[dir];
    n2[N] = 1;

    /* set X-coords in plane */
    Cp[0] = alloc_array(n2);
    Cp[1] = alloc_array(n2);
    fill_2arrays_with_nodepoints(node, N, Cp);
    /* convert Cp from Xb to X coords for node,
       these X are spread over the neighbor nodes */
    array_Xplane_of_Xb(node, N, Cp, Cp);
  }

  forplaneN(N, i,j,k, n, pl)
  {
    /* check if we have x on this node */
    if(!px[0])
    {
      ind = Ind_n_norm(i,j,k, n, N);
      /* set X and then get x */
      X[dir1] = Cp[0]->d[ind];
      X[dir2] = Cp[1]->d[ind];
      set_xyz(0, node,-1, X, x);
    }
    else
    {
      ind = Ind_n(i,j,k, n);
      for(dir=0; dir<3; dir++) x[dir] = px[dir][ind];
    }

    for(d2=0., dir=0; dir<3; dir++)
    {
      d1 = x[dir] - x0[dir];
      d2 += d1*d1;
    }
    if(d2<dist2)
    {
      ijk[0] = i;
      ijk[1] = j;
      ijk[2] = k;
      dist2 = d2;
    }
  }

  /* free mem */
  if(!px[0])
  {
    free_array(Cp[1]);
    free_array(Cp[0]);
  }

  return sqrt(dist2);
}

/* find i,j,k of corner closest to x0 in plane pl normal to N,
   return distance */
double nearest_corner_of_xyz_inplaneN(tNode *node, int N, int pl,
                                      int ijk[3], const double x0[3])
{
  tMesh *mesh = node->pat->mesh;
  int ix = Ind("x");
  double *px[] = { Vard(node,ix), Vard(node,ix+1), Vard(node,ix+2) };
  int *n = node->n;
  int dir1 = Dir1_norm(N);
  int dir2 = Dir2_norm(N);
  double Xpl=0, X[3], x[3];
  int dir, i,j,k, ind;
  double d1, d2, dist2 = DBL_MAX;

  /* if this node has no dat we need to compute x from Xb for each point */
  if(!px[0])
  {
    double Xb = node_Xb(node,N)->d[pl]; /* Xb coord in plane pl */

    /* convert Xb to X coords */
    X_of_Xb_indir(node, N, Xb, &Xpl);
    X[N] = Xpl;
  }

  forcornersN(N, i,j,k, n, pl)
  {
    /* check if we have x on this node */
    if(!px[0])
    {
      int i1 = i1_norm(i,j,k, N);
      int i2 = i2_norm(i,j,k, N);
      double Xb;

      /* set X */
      Xb = node_Xb(node,dir1)->d[i1];
      X_of_Xb_indir(node, dir1, Xb, &(X[dir1]));
      Xb = node_Xb(node,dir2)->d[i2];
      X_of_Xb_indir(node, dir2, Xb, &(X[dir2]));
      /* get x from X */
      set_xyz(0, node,-1, X, x);
    }
    else
    {
      ind = Ind_n(i,j,k, n);
      for(dir=0; dir<3; dir++) x[dir] = px[dir][ind];
    }

    for(d2=0., dir=0; dir<3; dir++)
    {
      d1 = x[dir] - x0[dir];
      d2 += d1*d1;
    }
    if(d2<dist2)
    {
      ijk[0] = i;
      ijk[1] = j;
      ijk[2] = k;
      dist2 = d2;
    }
  }

  return sqrt(dist2);
}

/* return an XYZ normal direction that is in a similar direction as the
   Cartesion normal direction cartN, returns -1 if not found */
int approxXYZnormal_of_xyznormal(tNode *node, int cartN)
{
  tPat *pat = node->pat;
  tCoordInfo *CI = pat->CI;

  switch(CI->type)
  {
  case PyramidFrustum:
  case innerCubedSphere:
  case outerCubedSphere:
  case CubedShell:
    if(cartN==2)
    {
      if(CI->dom<4) return 2;
      else          return -1;
    }
    else if(cartN==1)
    {
      if(CI->dom>=4 || CI->dom<=1) return 1;
      else                         return -1;
    }
    else if(cartN==0)
    {
      if(CI->dom>=4) return 2;
      if(CI->dom>=2) return 1;
      else           return -1;
    }
    return -1;
  default:
    return cartN; /* assume Cartesian coords */
  }
}

/* |\vec{x}|^2 */
double magnitude2_xyz(const double x[3])
{
  double d2;
  int dir;

  for(d2=0., dir=0; dir<3; dir++) d2 += x[dir] * x[dir];
  return d2;
}
/* |\vec{x}| */
double magnitude_xyz(const double x[3])
{
  return sqrt( magnitude2_xyz(x) );
}


/* get X from Xb */
void XYZ_of_XbYbZb(tNode *node, const double Xb[3], double X[3])
{
  double *nbb = node->bbox;
  int dir;

  for(dir=0; dir<3; dir++)
  {
    int f = dir*2;
    X[dir] = 0.5*( (nbb[f+1] - nbb[f]) * Xb[dir] + (nbb[f+1] + nbb[f]) );
  }
}

/* get dX/dXb */
void dXYZ_dXbYbZb(tNode *node, double dXdXb[3])
{
  double *nbb = node->bbox;
  int dir;

  for(dir=0; dir<3; dir++)
  {
    int f = dir*2;
    dXdXb[dir] = 0.5*( nbb[f+1] - nbb[f] );
  }
}
/* get dX/dXb */
void dXbYbZb_dXYZ(tNode *node, double dXbdX[3])
{
  double *nbb = node->bbox;
  int dir;

  for(dir=0; dir<3; dir++)
  {
    int f = dir*2;
    dXbdX[dir] = 2./( nbb[f+1] - nbb[f] );
  }
}

/* get one X from Xb in one direction */
void X_of_Xb_indir(tNode *node, int dir, double Xb, double *X)
{
  double *nbb = node->bbox;
  int f = dir*2;
  *X = 0.5*( (nbb[f+1] - nbb[f]) * Xb + (nbb[f+1] + nbb[f]) );
}

/* get X from Xb for entire arrays */
void array_XYZ_of_XbYbZb(tNode *node, tArray *aXb[3], tArray *aX[3])
{
  int dir, k;
  for(dir=0; dir<3; dir++)
  {
    int Nm = min2(aXb[dir]->N, aX[dir]->N);
    for(k=0; k<Nm; k++)
    {
      double Xb, X;
      Xb = aXb[dir]->d[k];
      X_of_Xb_indir(node, dir, Xb, &X);
      aX[dir]->d[k] = X;
    }
  }
}

/* get X from Xb for arrays in plane perp to dir */
void array_Xplane_of_Xb(tNode *node, int dir, tArray *aCb[2], tArray *aC[2])
{
  int d, d3, k, Nm;
  for(d=0; d<2; d++)
  {
    switch(dir)
    {
    case 0:
      d3 = d+1;
      break;
    case 1:
      d3 = d*2;
      break;
    case 2:
      d3 = d;
      break;
    default:
      errorexit("dir must be 0,1,2");
    }
    Nm = min2(aCb[d]->N, aC[d]->N);
    for(k=0; k<Nm; k++)
    {
      double Xb, X;
      Xb = aCb[d]->d[k];
      X_of_Xb_indir(node, d3, Xb, &X);
      aC[d]->d[k] = X;
    }
  }
}

/* get X from i,j,k */
void XYZ_of_ijk(tNode *node, int i, int j, int k, double X[3])
{
  XbYbZb_of_ijk(node, i,j,k, X);
  XYZ_of_XbYbZb(node, X, X);
}

/* get X from point index ind */
void XYZ_of_ind(tNode *node, int ind, double X[3])
{
  XbYbZb_of_ind(node, ind, X);
  XYZ_of_XbYbZb(node, X, X);
}

/* get x from point index ind */
void xyz_of_ind(tNode *node, int ind, double x[3])
{
  if(node->dat->coords_set)
  {
    int ix = coordinates->ix; // Ind("x");
    int dir;
    for(dir=0; dir<3; dir++) x[dir] = Vard_(node, ix+dir)[ind];
  }
  else
  {
    double X[3];
    XYZ_of_ind(node, ind, X);
    set_xyz(NULL, node, -1, X, x); //-1 index means don't read x,y,z vars
  }
}

/* get Xb from X */
void XbYbZb_of_XYZ(tNode *node, double Xb[3], const double X[3])
{
  double *nbb = node->bbox;
  int dir;

  for(dir=0; dir<3; dir++)
  {
    int f = dir*2;
    Xb[dir] = ( 2.*X[dir] - (nbb[f+1] + nbb[f]) )/(nbb[f+1] - nbb[f]);
  }
}

/* get one Xb from X in one direction */
void Xb_of_X_indir(tNode *node, int dir, double *Xb, const double X)
{
  double *nbb = node->bbox;
  int f = dir*2;
  *Xb = ( 2.*X - (nbb[f+1] + nbb[f]) )/(nbb[f+1] - nbb[f]);
}

/* get Xb from X for entire arrays */
void array_XbYbZb_of_XYZ(tNode *node, tArray *aXb[3], tArray *aX[3])
{
  int dir, k;
  for(dir=0; dir<3; dir++)
  {
    int Nm = min2(aXb[dir]->N, aX[dir]->N);
    for(k=0; k<Nm; k++)
    {
      double Xb, X;
      X = aX[dir]->d[k];
      Xb_of_X_indir(node, dir, &Xb, X);
      aXb[dir]->d[k] = Xb;
    }
  }
}

/* get Xb from X for arrays in plane perp to dir */
void array_Xbplane_of_X(tNode *node, int dir, tArray *aCb[2], tArray *aC[2])
{
  int d, d3, k, Nm;
  for(d=0; d<2; d++)
  {
    switch(dir)
    {
    case 0:
      d3 = d+1;
      break;
    case 1:
      d3 = d*2;
      break;
    case 2:
      d3 = d;
      break;
    default:
      errorexit("dir must be 0,1,2");
    }
    Nm = min2(aCb[d]->N, aC[d]->N);
    for(k=0; k<Nm; k++)
    {
      double Xb, X;
      X = aC[d]->d[k];
      Xb_of_X_indir(node, d3, &Xb, X);
      aCb[d]->d[k] = Xb;
    }
  }
}

/* is XYZ inside a node, also sets X to boundary value if it is very close
   to the boundary */
int XYZ_is_in_node(tNode *node, double X[3])
{
  double *nbb = node->bbox;
  int dir;

  /* return 0 if X is outside node */
  for(dir=0; dir<3; dir++)
  {
    int f = dir*2;
    if(dless(X[dir],nbb[f]))      return 0;
    if(dgreater(X[dir],nbb[f+1])) return 0;
  }

  /* set X to boundary value if it is very close to boundary */
  for(dir=0; dir<3; dir++)
  {
    int f = dir*2;
    if(X[dir] < nbb[f])   X[dir] = nbb[f];
    if(X[dir] > nbb[f+1]) X[dir] = nbb[f+1];
  }
  return 1;
}

/* Mark all points in aXP[0..2] within node by writing their index into
   aI. If a point is not in the node we write -1 into aI. I.e. we return
   a mask of points in node in Ip. */
void array_find_XYZ_in_node(tNode *node, tArray *aXP[3], tArray *aI)
{
  int k;
  forarray(aXP[0], k)
  {
    double X[] = { aXP[0]->d[k], aXP[1]->d[k], aXP[2]->d[k] };
    if(XYZ_is_in_node(node, X)) aI->i[k] = k;
    else                        aI->i[k] = -1;
  }
}


/* is a point C (in XYZ coords) in a plane (normal to dir) inside a node,
   also sets C to boundary value if it is very close to the boundary */
int Xplane_is_in_node(tNode *node, int dir, double C[2])
{
  double *nbb = node->bbox;
  int d, f;

  /* return 0 if C is outside node */
  for(d=0; d<2; d++)
  {
    switch(dir)
    {
    case 0:
      f = 2*d + 2;
      break;
    case 1:
      f = 4*d;
      break;
    case 2:
      f = 2*d;
      break;
    default:
      errorexit("dir has to be 0,1,2");
    }
    if(dless(C[d],nbb[f]))      return 0;
    if(dgreater(C[d],nbb[f+1])) return 0;
  }

  /* set X to boundary value if it is very close to boundary */
  for(d=0; d<2; d++)
  {
    switch(dir)
    {
    case 0:
      f = 2*d + 2;
      break;
    case 1:
      f = 4*d;
      break;
    case 2:
      f = 2*d;
    }
    if(C[d] < nbb[f])   C[d] = nbb[f];
    if(C[d] > nbb[f+1]) C[d] = nbb[f+1];
  }
  return 1;
}

/* Mark all points in aCP[0..1] in plane p normal to dir within node by
   writing their index into aI. If a point is not in the node we write -1
   into aI. I.e. we return a mask of points in node in Ip. */
void array_find_Xplane_in_node(tNode *node,int dir, tArray *aCP[2], tArray *aI)
{
  int k;
  forarray(aCP[0], k)
  {
    double C[]  = { aCP[0]->d[k], aCP[1]->d[k] };
    if(Xplane_is_in_node(node, dir, C)) aI->i[k] = k;
    else                                aI->i[k] = -1;
  }
}


/* find index of nb node in the node->fnb[f] list that contains
   the surface point C on face o_f of patch o_pat,
   return index of nb in fnb, or -1 if not found */
int fnb_containing_point(tNode *node, int f,
                         tPat *o_pat, int o_f, double C[2])
{
  int odir = o_f/2;
  int od1 = Dir1_norm(odir);
  int od2 = Dir2_norm(odir);
  int nfnb = node->nfnb[f];
  int i;

  for(i=0; i<nfnb; i++)
  {
    tNode *nb = node->fnb[f][i];
    double nbrct[] = { nb->bbox[2*od1], nb->bbox[2*od1+1],
                       nb->bbox[2*od2], nb->bbox[2*od2+1] };

    if(nb->pat != o_pat) continue;

    if(C_in_brct(nbrct , C)) return i;
  }
  return -1;
}


/* Mark all points in aC[0..1] of node that are in node nb by writing
   their index into aI. If a point is not in the node we write -1
   into aI. I.e. we return a mask of points in node in Ip.
   Note that aC and aoC contain the same points but in different coords.
   aC in X coords of the node and aoC in X coords of the neighbor nb. */
void mark_points_in_nb_f(tNode *node, int f, tArray *aC[2],
                         tNode *nb, int nb_f, tArray *aoC[2], tArray *aI)
{
  //tPat *pat = node->pat;
  //tPat *opat = nb->pat;
  int odir, od1, od2;
  double nbrct[4];
  int k;

  errorexit("mark_points_in_nb_f checks only 2 out of 3 coords!!! "
            "This is not good enough in many cases-> FIXME!!! "
            "Probably bound rect. should not be used for this...");

  /* get bounding rectangle of the node nb */
  odir = nb_f/2;
  od1 = Dir1_norm(odir);
  od2 = Dir2_norm(odir);
  nbrct[0] = nb->bbox[2*od1];
  nbrct[1] = nb->bbox[2*od1+1];
  nbrct[2] = nb->bbox[2*od2];
  nbrct[3] = nb->bbox[2*od2+1];

  forarray(aC[0], k)
  {
    //double C[]  = {  aC[0]->d[k],  aC[1]->d[k] };
    double oC[] = { aoC[0]->d[k], aoC[1]->d[k] };

    /* check if oC from opat is within bounding rectangle of nb */
    if( C_in_brct(nbrct, oC) )
      aI->i[k] = k;
    else
      aI->i[k] = -1;
  }
}


/* return patch index if x is inside this patch, if not return -1 */
int p_XYZ_of_xyz(tPat *pat, double X[3], const double x[3])
{
  int d, stat=0;

  /* get X */
  if(pat->XYZ_of_xyz)
    stat = pat->XYZ_of_xyz(pat, 0,-1, X, x);
  else
    for(d=0; d<3; d++) X[d] = x[d];

  //PRF;pr3v(": x",x);
  //printf("p=%d stat=%d ", pat->p, stat);pr3v("X",X);

  if(stat) return -1;

  for(d=0; d<3; d++)
    if(dless(X[d],pat->bbox[2*d]) || dless(pat->bbox[2*d+1],X[d]))
      return -1;

  //PRF;prbbox(pat->bbox,3);

  /* round X to inside box */
  for(d=0; d<3; d++)
  {
    if(X[d] < pat->bbox[2*d])   X[d] = pat->bbox[2*d];
    if(X[d] > pat->bbox[2*d+1]) X[d] = pat->bbox[2*d+1];
  }
  return pat->p;
}

/* go over pat list and find the one that contains x */
int p_XYZ_of_xyz_inpatlist(tMesh *mesh, intList *pl,
                           double X[3], const double x[3])
{
  tPat *pat;
  int i, p=-1;

  for(i=0; i<pl->n; i++)
  {
    p = pl->e[i];
    pat = mesh->pat[p];
    p = p_XYZ_of_xyz(pat, X, x);
    if(p>=0) break;
  }
  return p;
}

/* go over all patches and find the one that contains x */
int p_XYZ_of_xyz_mesh(tMesh *mesh, double X[3], const double x[3])
{
  int pi, p=-1;

  forpatches(mesh, pi)
  {
    tPat *pat = mesh->pat[pi];
    p = p_XYZ_of_xyz(pat, X, x);
    if(p>=0) break;
  }
  return p;
}

/* go over all patches and nodes and find the node name that contains x
   in:  mesh, namsiz, x
   out: name, X
   returns:  p */
int p_nodename_XYZ_of_xyz_mesh__OLD(tMesh *mesh, char *name, const int namsiz,
                               double X[3], const double x[3])
{
  int size = nMPI_size();
  int rank = nMPI_rank();
  int r, p;
  char *found; /* array with ranks that have x */
  char *found_local;

  /* set name="" as default */
  name[0] = 0;

  /* find patch p and set X */
  p = p_XYZ_of_xyz_mesh(mesh, X, x);
  //PRFs(": ");pr3v("X", X);printf(": p=%d\n", p);

  /* if x is not on mesh return -1 and leave name="" */
  if(p<0) return p;

  /* search among my leaf nodes (in patch p) for X */
  formylnodes_noomp(mesh)
  {
    tNode *node = MyLnode;
    if(node->pat->p == p)
      if(XYZ_is_in_node(node, X))
      {
        nodename(node, name, namsiz);
        break;
      }
  }

  /* mark found_local[rank] if we found X in one of my leaf nodes */
  found_local = calloc(size, sizeof(found_local[0]));
  found       = calloc(size, sizeof(found[0]));
  found_local[rank] = found[rank] = name[0];

  /* get global found */
  MCK( nMPI_Allreduce(found_local, found, size, nMPI_CHAR, nMPI_LOR) );

  /* find lowest rank r that has node with X */
  for(r=0; r<size; r++)
  {
    if(found[r]) break;
  }
  if(r>=size)
  {
    PRF;printf(": error: one rank must have the node with this x!\n");
    r = 0; /* to avoid failure in nMPI_Bcast */
  }

  /* broadcast node name from rank r to all MPI jobs */
  //PRF;printf(":|%s|r=%d\n", name, r);fflush(stdout);
  MCK( nMPI_Bcast(name, namsiz, nMPI_CHAR, r) );

  free(found);
  free(found_local);
  return p;
}

/* return node and set X to where x is located */
tNode *node_XYZ_of_xyz_mesh__OLD(tMesh *mesh, double X[3], const double x[3])
{
  char name[99];

  if(p_nodename_XYZ_of_xyz_mesh__OLD(mesh, name,99, X, x) < 0)
    return NULL;

  if(name[0]==0)
    return NULL;

  return node_from_nodename(mesh, name);
}

// replaces p_nodename_XYZ_of_xyz_mesh:
/* go over all elms on all ranks and find the elm eid name that contains x
   In:  mesh, x   Out: eid, X
   Returns:  p */
int p_eid_XYZ_of_xyz_mesh(tMesh *mesh, ulong *eid,
                          double X[3], const double x[3])
{
  int size = nMPI_size();
  int rank = nMPI_rank();
  int r, p;
  char *found; /* array with ranks that have x */
  char *found_local;
  char flag;

  /* if eid stays at this val we didn't find it */
  *eid = EID_INVALID;

  /* find patch p and set X */
  p = p_XYZ_of_xyz_mesh(mesh, X, x);
  //PRFs(": ");pr3v("X", X);printf(": p=%d\n", p);

  /* if x is not on mesh return -1 and leave name="" */
  if(p<0) return p;

  /* search among my leaf nodes (in patch p) for X */
  flag = 0;
  formyelms_noomp(mesh)
  {
    tElm *elm = MyElm;
    if(elm->pat->p == p)
      if(XYZ_is_in_node(elm, X))
      {
        *eid = Elm_eid(elm);
        flag = 1;
        break;
      }
  }

  /* mark found_local[rank] if we found X in one of my leaf nodes */
  found_local = checked_calloc(size, sizeof(found_local[0]));
  found       = checked_calloc(size, sizeof(found[0]));
  found_local[rank] = found[rank] = flag;

  /* get global found */
  MCK( nMPI_Allreduce(found_local, found, size, nMPI_CHAR, nMPI_LOR) );

  /* find lowest rank r that has node with X */
  for(r=0; r<size; r++)
  {
    if(found[r]) break;
  }
  if(r>=size)
  {
    PRF;printf(": error: one rank must have the node with this x!\n");
    r = 0; /* to avoid failure in nMPI_Bcast */
  }

  /* broadcast eid from rank r to all MPI jobs */
  MCK( nMPI_Bcast(eid,1, nMPI_UNSIGNED_LONG, r) );
  //PRF;printf(": *eid=%lu r=%d\n", *eid, r);fflush(stdout);

  free(found);
  free(found_local);
  return p;
}

// replaces node_XYZ_of_xyz_mesh
/* return elm and set eid,elmindex,elmrank,X to where x is located */
tElm *elm_XYZ_of_xyz_mesh(tMesh *mesh,
                          ulong *eid, ulong *elmindex, int *elmrank,
                          double X[3], const double x[3])
{
  tElm *elm;

  if(p_eid_XYZ_of_xyz_mesh(mesh, eid, X, x) < 0)
    return NULL;

  if(*eid == EID_INVALID)
    return NULL;

  elm = elm_from_eid(mesh, *eid, elmindex, elmrank);
  //PRF;printf(": *eid=%lu\n", *eid);
  //printf("*elmindex=%lu *elmrank=%d\n", *elmindex, *elmrank);
  return elm;
}

/* return node and set X to where x is located */
tNode *node_XYZ_of_xyz_mesh(tMesh *mesh, double X[3], const double x[3])
{
  ulong eid, elmindex;
  int elmrank;
  return elm_XYZ_of_xyz_mesh(mesh, &eid,&elmindex,&elmrank, X, x);
}

/* set elm0,eid,elmindex,elmrank,X to where x is located
   In: mesh,x   Out: elm0,eid,elmindex,elmrank,X
   NOTE: elm0,eid,elmindex,elmrank,X are invalid if eid==EID_INVALID
   returns: elm if elm is on this rank, NULL otherwise */
tElm *set_elm0_XYZ_of_xyz_mesh(tMesh *mesh, tElm0 elm0[1],
                               ulong *eid, ulong *elmindex, int *elmrank,
                               double X[3], const double x[3])
{
  tElm *elm = elm_XYZ_of_xyz_mesh(mesh, eid,elmindex,elmrank, X, x);

  /* if we have elm set elm0 from it */
  if(elm) memcpy(elm0, elm, sizeof(elm0[1]));

  /* if we have a valid eid, we should have valid data, and we thus
     Bcast elm0 from elmrank to all ranks */
  if(*eid != EID_INVALID)
    MCK( nMPI_Bcast(elm0, sizeof(tElm0), nMPI_CHAR, *elmrank) );

  return elm;
}


// replaces l_XYZ_of_xyz__old:
/* set X and return 1 if x is inside this elm, otherwise return 0 */
int elmXYZ_of_xyz(tElm *elm, int ind, double X[3], const double x[3])
{
  tPat *pat = elm->pat;
  int d, stat=0;

  /* get X */
  if(pat->XYZ_of_xyz)
    //stat = pat->XYZ_of_xyz(pat, (tNode *)elm,ind, X, x);
    stat = pat->XYZ_of_xyz(pat, elm,ind, X, x);
  else
    for(d=0; d<3; d++) X[d] = x[d];

  if(stat) return 0;

  for(d=0; d<3; d++)
    if(dless(X[d],elm->bbox[2*d]) || dless(elm->bbox[2*d+1],X[d]))
      return 0;

  /* round X to inside box */
  for(d=0; d<3; d++)
  {
    if(X[d] < elm->bbox[2*d])   X[d] = elm->bbox[2*d];
    if(X[d] > elm->bbox[2*d+1]) X[d] = elm->bbox[2*d+1];
  }

  return 1;
}

/* find the faces a point X is on within tol, face[2]=1 if X is on face2  */
int XYZ_on_face_tol(tPat *pat, int *face, const double X[3], double tol)
{
  double *bb=pat->bbox;
  double diag = pat->bbdiag;
  int f;
  int nf;

  /* find all faces we are on */
  for(nf=0, f=0; f<6; f++)
  {
    int d=f/2;
    if(dequal_tol(X[d], bb[f], tol*diag)) { face[f]=1; nf++; }
    else                                  { face[f]=0; }
  }
  return nf; /* number of faces point is on */
}

/* find the faces a point X is on, face[2]=1 if X is on face2  */
int XYZ_on_face(tPat *pat, int *face, const double X[3])
{
  return XYZ_on_face_tol(pat, face, X, 1e-10);
}

/* find the elm-faces a point X is on within tol:
   e.g. face[2]=1 if X is on face2
   Returns: number of faces point is on  */
int XYZ_on_elmface_tol(tElm *elm, int *face, const double X[3], double tol)
{
  tPat *pat   = elm->pat;
  double *bb  = elm->bbox;
  double diag = pat->bbdiag;
  //int l       = Elm_l(elm);
  //ulong fac   = 1<<l;
  //double diag = pat->bbdiag / fac;
  int f;
  int nf;

  /* find all faces we are on */
  for(nf=0, f=0; f<6; f++)
  {
    int d=f/2;
    if(dequal_tol(X[d], bb[f], tol*diag)) { face[f]=1; nf++; }
    else                                  { face[f]=0; }
  }
  return nf; /* number of faces point is on */
}

/* find the elm-faces a point X is on:  e.g. face[2]=1 if X is on face2
   Returns: number of faces point is on  */
int XYZ_on_elmface(tElm *elm, int *face, const double X[3])
{
  return XYZ_on_elmface_tol(elm, face, X, 1e-10);
}

/* check if i,j,k is on any elm face */
int ijk_on_elmface(tElm *elm, int i, int j, int k, int *face)
{
  int *n = elm->n;
  int In[] = { i,j,k };
  int f, nf;

  for(nf=0, f=0; f<6; f++)
  {
    int d = f/2;
    int pl = (n[d]-1)*(f%2);
    if(In[d]==pl) { face[f] = 1; nf++; }
    else          { face[f] = 0; }
  }
  if(0)
  {
    printf("%d %d %d  ", In[0],In[1],In[2]);
    for(f=0; f<6; f++) printf("%d ", face[f]);
    printf(" -> nf=%d\n ", nf);
  }
  return nf;
}

/* check if ind is on any elm face */
int ind_on_elmface(tElm *elm, int ind, int *face)
{
  int *n = elm->n;
  int k  = kOfInd_n(ind, n);
  int j  = jOfInd_n_k(ind, n, k);
  int i  = iOfInd_n_jk(ind, n, j,k);
  return ijk_on_elmface(elm, i,j,k, face);
}

/* check if ind is on a node face */
/*
int ind_on_nodeface(tNode *node, int ind, int *face)
{
  int *n = node->n;
  int In[3];
  int f, nf;

  In[2] = kOfInd_n(ind, n);
  In[1] = jOfInd_n_k(ind, n, In[2]);
  In[0] = iOfInd_n_jk(ind, n, In[1],In[2]);

  for(nf=0, f=0; f<6; f++)
  {
    int d = f/2;
    int pl = (n[d]-1)*(f%2);
    if(In[d]==pl) { face[f] = 1; nf++; }
    else          { face[f] = 0; }
  }
  if(0)
  {
    printf("%d: %d %d %d  ", ind, In[0],In[1],In[2]);
    for(f=0; f<6; f++) printf("%d ", face[f]);
    printf(" -> nf=%d\n ", nf);
  }
  return nf;
}
*/

/* check if ind is on outer boundary */
int ind_on_outerbound(tNode *node, int ind)
{
  int face[6];

  if(ind_on_elmface(node, ind, face))
  {
    int f;
    for(f=0; f<6; f++)
      if(face[f])
      {
        //tBface *bfaces = node->pat->bfaces[f];
        //if(node->patface[f] && bfaces && bfaces->boundary==OUTERBOUND)
        if(Elm_on_OUTERBOUND(node,f))
          return 1;
      }
  }
  return 0;
}


/* set x at X */
int set_xyz(tPat *pat, tNode *node, int ind, const double X[3], double x[3])
{
  if(!pat) pat = node ? node->pat : NULL;
  if(!pat) errorexit("pat and node must not both be NULL");

  /* now set x, dXb/dx */
  if(pat->xyz_of_XYZ)
  {
    return pat->xyz_of_XYZ(pat, node, ind, X, x);
  }
  else /* assume X,Y,Z are Cartesian*/
  {
    int d;

    for(d=0; d<3; d++)
      x[d] = X[d];
    return 0;
  }
}

/* set X at x */
int set_XYZ(tPat *pat, tNode *node, int ind, double X[3], const double x[3])
{
  if(!pat) pat = node ? node->pat : NULL;

  /* now set x, dXb/dx */
  if(pat->XYZ_of_xyz)
  {
    return pat->XYZ_of_xyz(pat, node, ind, X, x);
  }
  else /* assume X,Y,Z are Cartesian*/
  {
    int d;

    for(d=0; d<3; d++)
      X[d] = x[d];
    return 0;
  }
}

/* set x and dXYZ/dxyz at X */
int set_xyz_dXYZdxyz(tPat *pat, tNode *node, int ind,
                     const double X[3], double x[3], double dXYZdxyz[3][3])
{
  /* now set x, dXdx, det(dXb/dx) */
  if(pat->dXYZ_dxyz)
  {
    return pat->dXYZ_dxyz(pat, node, ind, X, x, dXYZdxyz);
  }
  else /* assume X,Y,Z are Cartesian */
  {
    int d, e;

    for(d=0; d<3; d++)
    {
      x[d] = X[d];
      for(e=0; e<3; e++)
        dXYZdxyz[d][e] = (d==e);
    }
    return 0;
  }
}

/* get the bounding rectangle of a nodeface, norm = face/2 */
void brct_nodeface(tNode *node, int norm, double brct[4])
{
  int i;
  switch(norm)
  {
  case 0:
    for(i=0; i<4; i++) brct[i] = node->bbox[i+2];
    break;
  case 1:
    for(i=0; i<2; i++) brct[i] = node->bbox[i];
    for(i=2; i<4; i++) brct[i] = node->bbox[i+2];
    break;
  case 2:
    for(i=0; i<4; i++) brct[i] = node->bbox[i];
    break;
  default:
    errorexit("norm must be 0,1,2");
  }
}

/* get the bounding rectangle of a patchface, norm = face/2 */
void brct_patface(tPat *pat, int norm, double brct[4])
{
  int i;
  switch(norm)
  {
  case 0:
    for(i=0; i<4; i++) brct[i] = pat->bbox[i+2];
    break;
  case 1:
    for(i=0; i<2; i++) brct[i] = pat->bbox[i];
    for(i=2; i<4; i++) brct[i] = pat->bbox[i+2];
    break;
  case 2:
    for(i=0; i<4; i++) brct[i] = pat->bbox[i];
    break;
  default:
    errorexit("norm must be 0,1,2");
  }
}

/* change brct by an amount eps */
void resize_brct(double brct[4], double eps)
{
  double L0 = brct[1] - brct[0];
  double L1 = brct[3] - brct[2];
  double dL0 = L0 * eps * 0.5;
  double dL1 = L1 * eps * 0.5;

  brct[0] = brct[0] + dL0;
  brct[1] = brct[1] - dL0;
  brct[2] = brct[2] + dL1;
  brct[3] = brct[3] - dL1;
}

/* make bounding rectangle large enough to fit the point X[3],
   if expand=0, include just this one point */
void expand_brct_to_include_X(double brct[4], int norm,
                              const double X[3], int expand)
{
  int d;
  double C[2]; /* point coords in face */

  switch(norm)
  {
  case 0:
    C[0] = X[1];  C[1] = X[2];
    break;
  case 1:
    C[0] = X[0];  C[1] = X[2];
    break;
  case 2:
    C[0] = X[0];  C[1] = X[1];
    break;
  default:
    errorexit("norm must be 0,1,2");
  }
  if(expand)
  {
    /* expand rectangle */
    for(d=0; d<2; d++)
    {
      if(C[d] < brct[2*d])   brct[2*d]   = C[d];
      if(C[d] > brct[2*d+1]) brct[2*d+1] = C[d];
    }
  }
  else
  {
    brct[1] = brct[0] = C[0];
    brct[3] = brct[2] = C[1];
   }
}

/* put intersection of 2 bounding rectangles into brct, if intersection
   is empty return 0 */
int intersection_brct1_brct2(const double brct1[4], const double brct2[4],
                             double brct[4])
{
  int d, isec=1;

  for(d=0; d<2; d++)
  {
    if((brct1[2*d] >= brct2[2*d]) && (brct1[2*d] <= brct2[2*d+1]))
    {
      brct[2*d] = brct1[2*d];
      if(brct1[2*d+1] < brct2[2*d+1])
        brct[2*d+1] = brct1[2*d+1];
      else
        brct[2*d+1] = brct2[2*d+1];
    }
    else if((brct2[2*d] >= brct1[2*d]) && (brct2[2*d] <= brct1[2*d+1]))
    {
      brct[2*d] = brct2[2*d];
      if(brct2[2*d+1] < brct1[2*d+1])
        brct[2*d+1] = brct2[2*d+1];
      else
        brct[2*d+1] = brct1[2*d+1];
    }
    else
    {
      isec = 0;
      break;
    }
  }

  /* test if the intersection is empty */
  if(isec)
    for(d=0; d<2; d++)
      if(dequal(brct[2*d], brct[2*d+1]))
      {
        isec = 0;
        break;
      }

  return isec;
}

/* put intersection of 2 bounding boxes into bb, if intersection
   is empty return 0 */
int intersection_bb1_bb2(const double bb1[6], const double bb2[6],
                         double bb[6])
{
  int d, isec=1;

  for(d=0; d<3; d++)
  {
    if((bb1[2*d] >= bb2[2*d]) && (bb1[2*d] <= bb2[2*d+1]))
    {
      bb[2*d] = bb1[2*d];
      if(bb1[2*d+1] < bb2[2*d+1])
        bb[2*d+1] = bb1[2*d+1];
      else
        bb[2*d+1] = bb2[2*d+1];
    }
    else if((bb2[2*d] >= bb1[2*d]) && (bb2[2*d] <= bb1[2*d+1]))
    {
      bb[2*d] = bb2[2*d];
      if(bb2[2*d+1] < bb1[2*d+1])
        bb[2*d+1] = bb2[2*d+1];
      else
        bb[2*d+1] = bb1[2*d+1];
    }
    else
    {
      isec = 0;
      break;
    }
  }

  /* test if the intersection is empty */
  if(isec)
    for(d=0; d<3; d++)
      if(dequal(bb[2*d], bb[2*d+1]))
      {
        isec = 0;
        break;
      }

  return isec;
}


/* put intersection of 2 bounding boxes into bb, if intersection
   is empty return 0. But return 1 if they touch approximately. */
int touch_or_intersect_bb1_bb2(const double bb1[6], const double bb2[6],
                               double bb[6])
{
  int d, isec=1;

  for(d=0; d<3; d++)
  {
    //if((bb1[2*d] >= bb2[2*d]) && (bb1[2*d] <= bb2[2*d+1]))
    if(dgreatereq(bb1[2*d], bb2[2*d]) && dlesseq(bb1[2*d], bb2[2*d+1]))
    {
      bb[2*d] = bb1[2*d];
      if(bb1[2*d+1] < bb2[2*d+1])
        bb[2*d+1] = bb1[2*d+1];
      else
        bb[2*d+1] = bb2[2*d+1];
    }
    //else if((bb2[2*d] >= bb1[2*d]) && (bb2[2*d] <= bb1[2*d+1]))
    else if(dgreatereq(bb2[2*d], bb1[2*d]) && dlesseq(bb2[2*d], bb1[2*d+1]))
    {
      bb[2*d] = bb2[2*d];
      if(bb2[2*d+1] < bb1[2*d+1])
        bb[2*d+1] = bb2[2*d+1];
      else
        bb[2*d+1] = bb1[2*d+1];
    }
    else
    {
      isec = 0;
      break;
    }
  }

  return isec;
}


/* transform XYZ from patch1 to XYZ of other patch2 */
int XYZpat2_of_XYZpat1(tPat *pat1, const double X1[3],
                        tPat *pat2, double X2[3])
{
  double x[3];
  int p2;

  set_xyz(pat1, NULL, -1, X1, x);
  p2 = p_XYZ_of_xyz(pat2, X2, x);
  return p2;
}

/* get Coords on face from X */
void C_from_X_on_face(const double X[3], int face, double C[2])
{
  switch(face/2)
  {
  case 0:
    C[0] = X[1];
    C[1] = X[2];
    break;
  case 1:
    C[0] = X[0];
    C[1] = X[2];
    break;
  case 2:
    C[0] = X[0];
    C[1] = X[1];
    break;
  default:
    errorexit("face/2 must be 0,1,2");
  }
}

/* get X from Coords on face */
void X_from_C_on_face(tPat *pat, int face, const double C[2], double X[3])
{
  switch(face/2)
  {
  case 0:
    X[0] = pat->bbox[face];
    X[1] = C[0];
    X[2] = C[1];
    break;
  case 1:
    X[0] = C[0];
    X[1] = pat->bbox[face];
    X[2] = C[1];
    break;
  case 2:
    X[0] = C[0];
    X[1] = C[1];
    X[2] = pat->bbox[face];
    break;
  default:
    errorexit("face/2 must be 0,1,2");
  }
}

/* transform point C1 on face f1 of pat1 to point C2 on face2 of pat2 */
int Cpat2_of_Cpat1(tPat *pat1, int f1, const double C1[2],
                    tPat *pat2, int f2, double C2[2])
{
  double X1[3], X2[3];
  int p2;

  X_from_C_on_face(pat1,f1, C1, X1);
  p2 = XYZpat2_of_XYZpat1(pat1, X1, pat2, X2);
  /* if(p2<0) X2 is not within pat2 or NAN */

  C_from_X_on_face(X2,f2, C2);

  return p2;
}

/* transform brct from face f1 of patch1 to brct of face f2 on patch2 */
int brctpat2_of_brctpat1__old(tPat *pat1, int f1, const double brct1[4],
                         tPat *pat2, int f2, double brct2[4])
{
  double C1[2], C2[2], sw;
  int p2_1, p2_2, problem = 0;

  C1[0] = brct1[0];
  C1[1] = brct1[2];
  p2_1 = Cpat2_of_Cpat1(pat1,f1,C1,   pat2,f2,C2);
  /* set lower bounds of brct2 */
  brct2[0] = C2[0];
  brct2[2] = C2[1];


  C1[0] = brct1[1];
  C1[1] = brct1[3];
  p2_2 = Cpat2_of_Cpat1(pat1,f1,C1, pat2,f2,C2);
  /* set upper bounds of brct2 */
  brct2[1] = C2[0];
  brct2[3] = C2[1];

  /* check if C2 is not NAN */
  if(p2_1<0 || p2_2<0)
  {
    int j;

    /* check if all is ok, or if there is a NAN */
    for(j=0; j<4; j++) if(isnan(brct2[j])) problem = 1;
  }

  /* swap brct2 entries to have min in brct2[0], brct2[2] */
  if(brct2[1] < brct2[0])
  {
    sw = brct2[1];
    brct2[1] = brct2[0];
    brct2[0] = sw;
  }
  if(brct2[3] < brct2[2])
  {
    sw = brct2[3];
    brct2[3] = brct2[2];
    brct2[2] = sw;
  }
  return problem;
}

/* transform brct from face f1 of patch1 to brct of face f2 on patch2 */
int brctpat2_of_brctpat1(tPat *pat1, int f1, const double brct1[4],
                         tPat *pat2, int f2, double brct2[4])
{
  double C1[5][2]; /* 5 points, 3rd one is (C1[2][0],C1[2][1]) */
  double C2[5][2]; /* 5 points, 3rd one is (C2[2][0],C2[2][1]) */
  double f[5];
  int i,j, im[5], p2[5], prob[5], problem = 0;

  C1[0][0] = brct1[0];
  C1[0][1] = brct1[2];
  C1[1][0] = brct1[1];
  C1[1][1] = brct1[2];
  C1[2][0] = brct1[0];
  C1[2][1] = brct1[3];
  C1[3][0] = brct1[1];
  C1[3][1] = brct1[3];
  C1[4][0] = 0.5*(brct1[0] + brct1[1]);
  C1[4][1] = 0.5*(brct1[2] + brct1[3]);

  for(i=0; i<5; i++)
    p2[i] = Cpat2_of_Cpat1(pat1,f1,C1[i],   pat2,f2,C2[i]);

  /* check if C2 is not NAN */
  if(p2[0]<0 || p2[1]<0 || p2[2]<0 || p2[3]<0 ||  p2[4]<0)
  {
    /* check if all is ok, or if there is a NAN */
    for(i=0; i<5; i++)
    {
      prob[i] = 0;

      for(j=0; j<2; j++)
      {
        if(isinf(C2[i][j]))      prob[i] = problem = 1;
        else if(isnan(C2[i][j])) prob[i] = problem = 10;
      }
    }
  }

  /* find min and max in coords */
  for(j=0, i=0; i<5; i++) if(!prob[i]) f[j++] = C2[i][0];
  brct2[0] = min_in_1d_array(f,j, im);

  for(j=0, i=0; i<5; i++) if(!prob[i]) f[j++] = C2[i][0];
  brct2[1] = max_in_1d_array(f,j, im);

  for(j=0, i=0; i<5; i++) if(!prob[i]) f[j++] = C2[i][1];
  brct2[2] = min_in_1d_array(f,j, im);

  for(j=0, i=0; i<5; i++) if(!prob[i]) f[j++] = C2[i][1];
  brct2[3] = max_in_1d_array(f,j, im);

  return problem;
}

/* check if a point is in brct,
   also sets C to boundary value if it is very close to the boundary */
int C_in_brct(const double brct[4], double C[2])
{
  int d;

  /* return 0 if C is clearly outside */
  for(d=0; d<2; d++)
  {
    if(dless(C[d],    brct[2*d]))   return 0;
    if(dgreater(C[d], brct[2*d+1])) return 0;
  }

  /* set C to boundary value if it is very close to boundary */
  for(d=0; d<2; d++)
  {
    if(C[d] < brct[2*d])   C[d] = brct[2*d];
    if(C[d] > brct[2*d+1]) C[d] = brct[2*d+1];
  }
  return 1;
}


/* set x,y,z in face-plane */
void set_xyz_in_face(tNode *node, int face, int i, int j, int k, double x[3])
{
  int dir = face/2;
  double Xb[3], X[3];

  /* get Xb and then X of i,j,k */
  XbYbZb_of_ijk(node, i,j,k, Xb);
  XYZ_of_XbYbZb(node, Xb, X);

  /* set one X to value on face and then get x */
  X[dir] = node->bbox[face];
  set_xyz(NULL, node, -1, X, x);
}


/* convert X-coords on face f of node to X-coords of neighboring node,
   write them into nbC, and mark points outside nb with -1 in nb I */
void array_find_nbXface_of_Xface(tNode *node, int f, tNode *nb, int nb_f,
                                 tArray *nbC[2], tArray *nbI)
{
  tPat *pat = node->pat;
  tMesh *mesh = pat->mesh;
  int *n = node->n;
  int dir = f/2;
  int n_dir = n[dir];
  int pl = (n_dir-1)*(f%2);
  int ix = Ind("x");
  double *px[] = { Vard(node,ix), Vard(node,ix+1), Vard(node,ix+2) };
  double *oC[] = { Arrd(nbC[0]), Arrd(nbC[1]) };
  int *oI = Arri(nbI);
  int i,j,k;
  long loc;
  int odir = nb_f/2;
  int od1 = Dir1_norm(odir);
  int od2 = Dir2_norm(odir);

  forplaneN(dir, i,j,k, n, pl)
  {
    int ijk = Ind_n(i,j,k, n);
    int ind = Ind_n_norm(i,j,k, n, dir);
    double x[] = { px[0][ijk], px[1][ijk], px[2][ijk] };
    double oX[3];

    /* if n_dir=1 the x above is not on the face,
       so we set it the a value on face f */
    if(n_dir<2) set_xyz_in_face(node, f, i,j,k, x);

    /* find point x in nb */
    //loc = l_XYZ_of_xyz(nb,-1, oX, x);
    //if(loc>=0) /* point was found inside nb */
    /* find point x in nb */
    loc = elmXYZ_of_xyz(nb,-1, oX, x);
    if(loc) /* point was found inside nb */
    {
      oC[0][ind] = oX[od1];
      oC[1][ind] = oX[od2];
      oI[ind] = ind;
    }
    else /* point is not inside opat */
      oI[ind] = -1;
  }
}


/* write node mid point at i,j,k in dir into Xbm in Xb-coords,
   return value is 1 if mid point is found
   return value is 0 otherwise, then Xbm is not set! */
int set_nodemidpoint_XbYbZb(tNode *node, int i, int j, int k, int dir,
                            double Xbm[3])
{
  int *n = node->n;
  int nm1 = n[dir]-1;

  /* return zero if there is only 1 regular grid point in dir */
  if(nm1<=0) return 0; /* we do not set Xbm in this case */

  switch(dir)
  {
  case 0:
    if(i<nm1)
    {
      Xbm[dir] = 0.5*(node_Xb(node,dir)->d[i] + node_Xb(node,dir)->d[i+1]);
      Xbm[1] = node_Xb(node,1)->d[j];
      Xbm[2] = node_Xb(node,2)->d[k];
      return 1;
    }
    return 0;
  case 1:
    if(j<nm1)
    {
      Xbm[0] = node_Xb(node,0)->d[i];
      Xbm[dir] = 0.5*(node_Xb(node,dir)->d[j] + node_Xb(node,dir)->d[j+1]);
      Xbm[2] = node_Xb(node,2)->d[k];
      return 1;
    }
    return 0;
  case 2:
    if(k<nm1)
    {
      Xbm[0] = node_Xb(node,0)->d[i];
      Xbm[1] = node_Xb(node,1)->d[j];
      Xbm[dir] = 0.5*(node_Xb(node,dir)->d[k] + node_Xb(node,dir)->d[k+1]);
      return 1;
    }
    return 0;
  default:
    errorexit("dir must be 0,1,2");
  }
}

/* Write nm node midpoints (starting at id) in dir into Xbmd in Xb-coords.
   If set_nodemidpoint_XbYbZb fails it returns 0. Then Xbm is never set, and
   thus Xbmd may contain invalid values. We signal this by returning the 0
   from set_nodemidpoint_XbYbZb in this case. */
int set_nm_nodemidpoints_Xb_dir(tNode *node, int nm, int id, int dir,
                                double *Xbmd)
{
  double ret = 0;
  double Xbm[3] = {0}; //unnecessary init to fix gcc warning
  int i;

  switch(dir)
  {
  case 0:
    for(i=0; i<nm; i++)
    {
      ret = set_nodemidpoint_XbYbZb(node, id+i,0,0, dir, Xbm);
      Xbmd[i] = Xbm[dir];
    }
    return ret;
  case 1:
    for(i=0; i<nm; i++)
    {
      ret = set_nodemidpoint_XbYbZb(node, 0,id+i,0, dir, Xbm);
      Xbmd[i] = Xbm[dir];
    }
    return ret;
  case 2:
    for(i=0; i<nm; i++)
    {
      ret = set_nodemidpoint_XbYbZb(node, 0,0,id+i, dir, Xbm);
      Xbmd[i] = Xbm[dir];
    }
    return ret;
  default:
    errorexit("dir must be 0,1,2");
  }
}

/* write node mid point at id in dir into Xbmd in Xb-coords */
int set_nodemidpoint_Xb_dir(tNode *node, int id, int dir, double *Xbmd)
{
  return set_nm_nodemidpoints_Xb_dir(node,1, id, dir, Xbmd);
}

/* write all node midpoints plus the points on the 2 faces in dir into
   Xbmdf in Xb-coords. These are n+1 values if the node has n gridpoints! */
int set_nodemidpoints_2facepoints_Xb_dir(tNode *node, int dir, double *Xbmdf)
{
  int n = node->n[dir];

  /* first set all n-1 interior midpoints */
  set_nm_nodemidpoints_Xb_dir(node, n-1,0, dir, Xbmdf+1);

  /* add face points as first and last */
  Xbmdf[0] = -1.;
  Xbmdf[n] = +1.;
  return n+1;
}

/* write min distance in Xb-coords of nodemidpoints into Xbdist,
   if there is no midpoint write -2 in distXb for this face */
void set_nodemidpoints_to_face_distXb(tNode *node, double distXb[6])
{
  int *n = node->n;
  double *bb = node->bbox;
  int f;

  for(f=0; f<6; f++)
  {
    double Xbm, Xb[3];
    double X[3] = {0}; /* init to avoid NANs in it */
    int dir = f/2;
    int right = f%2;
    int nm1 = n[dir]-1;
    int id = (nm1-1)*(right);

    if(nm1>0)
    {
      set_nodemidpoint_Xb_dir(node, id, dir, &Xbm);
      X[dir] = bb[f];
      XbYbZb_of_XYZ(node, Xb, X);
      distXb[f] = fabs(Xbm - Xb[dir]);
    }
    else /* if we have only 1 regular grid point return -(node width) */
    {
      distXb[f] = -2.;
    }
  }
}

/* this computes all distances between the n-1 midpoints Xbmd in direc. dir
   and the distances to the nodeface on both ends and puts them into
   dXb, this results in n different distances */
int set_nm_nodemidpoint_distsXb_dir(tNode *node, int dir,
                                    const double *Xbmid, double *dXb)
{
  double distXb[6];
  int n = node->n[dir];
  int i;

  if(n<=1)
  {
    dXb[0] = 2.;
    return 0;
  }

  /* dist of first and last midpoint to nodeface */
  set_nodemidpoints_to_face_distXb(node, distXb);
  dXb[0]   = distXb[dir*2];
  dXb[n-1] = distXb[dir*2+1];

  /* distances between two midpoints */
  for(i=1; i<n-1; i++)
  {
    dXb[i] = Xbmid[i] - Xbmid[i-1];
  }
  //printf("dXb[0]=%g dXb[2]=%g dXb[n-1]=%g\n", dXb[0], dXb[2], dXb[n-1]);

  return 1;
}

/* with fin. vol. we sometimes want to interpret the standard grid points
   that are on the faces as grid points at the cell center between the face
   and the nearest midpoint.
   This func shifts them in this way. */
void shift_Xb0_XbN_toward_Xbm0_XbmN(const double *Xbm, int n, double *Xb)
{
  int nm = n-1;
  Xb[0]  = 0.5*(Xb[0]   + Xbm[0]);
  Xb[nm] = 0.5*(Xb[nm] + Xbm[nm-1]);
}
