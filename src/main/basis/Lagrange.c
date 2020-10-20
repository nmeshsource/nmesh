/* Lagrange.c */
/* Wolfgang Tichy, 2/2018
   some functions to deal with Lagrange interpolation polynomials */

#include "nmesh.h"
#include "basis.h"


/* ************************************************************************ */
/* various functions needed for Lagrange interpolation                      */
/* ************************************************************************ */

/* get Lagrange interpolation weights w_interp from the n points in x,
   this is coming from the denominator in Lagrange interpolation only */
void Lagrange_winterp(int n, const double *x, double *w_interp)
{
  int m, i;
  double denom;

  for(i = 0; i < n; i++)
  {
    denom = 1.;
    for(m = 0; m < n; m++)
      if(m != i) denom *= (x[i] - x[m]);

    w_interp[i] = 1./denom;
  }
}


/* find matrix D for taking derivatives: D_{ij} = \partial_x l_j(x_i),
   this sets the transpose D^T if DT is interpreted as stored in
   column-major form */
void Lagrange_DT(int n, const double *x, const double *w_interp, double *DT)
{
  int i, j;
  double Dii, Dij;

  for(i = 0; i < n; i++)
  {
    Dii = 0;
    for(j = 0; j < n; j++)
    {
      if(i != j)
      {
        Dij  = (w_interp[j] / w_interp[i]) / (x[i] - x[j]);
        Dii -= Dij;
        DT[j + i*n] = Dij;
        /* NOTE: this DT is D_{ij} in row-major form or its transpose
                 D_{ji} in column-major form */
      }
    }
    DT[i*n + i] = Dii;
  }
}

/* get Lagrange basis function l_k(x),
   here x_p are the grid points, w_interp the interp. weights */
double Lagrange_of_x(int k, double x, int np,
                     const double *x_p, const double *w_interp)
{
  int m;
  double prod = 1.;

  for(m=0; m<np; m++) if(m!=k) prod *= (x - x_p[m]);
  return prod * w_interp[k];
}

/***********************************************************************/
/* interpolate */
/***********************************************************************/

/* 3d interpolation:
   interpolate to the point (Xb[0],Xb[1],Xb[2]) for variable in array var
   Note: for Lagrange interpolation the coeffs are simply the function
         values at grid points
   Note2: the only info we really retrieve from the node (in Xb3_n and WL3_n)
          is node->pt_typ, but not node->n */
double Lagrange_array_interpolate(tNode *node, tArray *var, double Xb[3])
{
  int *n = var->n;
  tArray *Xb_n[3];
  tArray *WL_n[3];
  double *xp[3]; /* points */
  double *w[3];  /* weights */
  double *restrict B0 = dmalloc(n[0]);   /* basis */
  double *restrict B1 = dmalloc(n[1]);
  double *restrict B2 = dmalloc(n[2]);
  int k;
  double sum;

  /* get arrays with points and weights for n = var->n */
  Xb3_n(node, n, Xb_n);
  WL3_n(node, n, WL_n);
  /* now set data pointers to points and weights */
  for(k=0; k<3; k++)
  {
    xp[k] = Xb_n[k]->d;
    w[k]  = WL_n[k]->d;
  }

  /* save basis func values at (Xb[0],Xb[1],Xb[2]) in B0,... */
  for(k=0; k<n[0]; k++) B0[k] = Lagrange_of_x(k, Xb[0], n[0], xp[0], w[0]);
  for(k=0; k<n[1]; k++) B1[k] = Lagrange_of_x(k, Xb[1], n[1], xp[1], w[1]);
  for(k=0; k<n[2]; k++) B2[k] = Lagrange_of_x(k, Xb[2], n[2], xp[2], w[2]);

  /* interpolate to (Xb[0],Xb[1],Xb[2]) */
  sum = 0.;
  //SGRID_LEVEL3_Pragma(omp parallel for reduction(+:sum))
  for(k=0; k<n[2]; k++)
  {
    int j,i;
    for(j=0; j<n[1]; j++)
    for(i=0; i<n[0]; i++)
      sum += var->d[Ind_n(i,j,k, n)] * B0[i] * B1[j] * B2[k];
  }

  free(B2);
  free(B1);
  free(B0);
  return sum;
}

/* 2d interpolation:
   interpolate to the point (Cb1, Cb2) for variable in array var
   in plane p orthogonal to direction dir
   NOTE: We can set node=neighbor when we call this, even if the var is not
         on neighbor. We can use this to interpolate on a surface that was
         copied from a neighbor node!
   Note2: the only info we really retrieve from the node (in Xb3_n and WL3_n)
          is node->pt_typ, but not node->n */
double Lagrange_array_interpolate2d(tNode *node, tArray *var, int dir, int p,
                                    double Cb[2])
{
  int *n = var->n;
  tArray *Xb_n[3];
  tArray *WL_n[3];
  double *xp[3]; /* points */
  double *w[3];  /* weights */
  double *restrict B1 = dmalloc(max3(n[0],n[1],n[2]));
  double *restrict B2 = dmalloc(max3(n[0],n[1],n[2]));
  int i,j,k;
  double sum;

  /* get arrays with points and weights for n = var->n */
  Xb3_n(node, n, Xb_n);
  WL3_n(node, n, WL_n);
  /* now set data pointers to points and weights */
  for(k=0; k<3; k++)
  {
    xp[k] = Xb_n[k]->d;
    w[k]  = WL_n[k]->d;
  }

  switch(dir)
  {
  case 0:
    /* save basis func values */
    for(k=0; k<n[1]; k++) B1[k] = Lagrange_of_x(k, Cb[0], n[1], xp[1], w[1]);
    for(k=0; k<n[2]; k++) B2[k] = Lagrange_of_x(k, Cb[1], n[2], xp[2], w[2]);

    /* interpolate */
    sum = 0.;
    for(k=0; k<n[2]; k++)
    for(j=0; j<n[1]; j++)
      sum += var->d[Ind_n(p,j,k, n)] * B1[j] * B2[k];
    break;
  case 1:
    /* save basis func values */
    for(k=0; k<n[0]; k++) B1[k] = Lagrange_of_x(k, Cb[0], n[0], xp[0], w[0]);
    for(k=0; k<n[2]; k++) B2[k] = Lagrange_of_x(k, Cb[1], n[2], xp[2], w[2]);

    /* interpolate */
    sum = 0.;
    for(k=0; k<n[2]; k++)
    for(i=0; i<n[0]; i++)
      sum += var->d[Ind_n(i,p,k, n)] * B1[i] * B2[k];
    break;
  case 2:
    /* save basis func values */
    for(k=0; k<n[0]; k++) B1[k] = Lagrange_of_x(k, Cb[0], n[0], xp[0], w[0]);
    for(k=0; k<n[1]; k++) B2[k] = Lagrange_of_x(k, Cb[1], n[1], xp[1], w[1]);

    /* interpolate */
    sum = 0.;
    for(j=0; j<n[1]; j++)
    for(i=0; i<n[0]; i++)
      sum += var->d[Ind_n(i,j,p, n)] * B1[i] * B2[j];
    break;
  default:
    errorexit("dir must be 0,1,2");
  }
  free(B2);
  free(B1);
  return sum;
}


/* make 3 arrays Xp[0..2] that contain all points of a node. The
   arrays Xp[0..2] are in Xb coords. */
void fill_3arrays_with_nodepoints(tNode *node, tArray *Xp[3])
{
  int i,j,k, dir, *n = node->n;
  forijk(i,j,k, n)
  {
    double Xb[] = { node_Xb(node,0)->d[i], node_Xb(node,1)->d[j],
                    node_Xb(node,2)->d[k] };
    for(dir=0; dir<3; dir++)
      Xp[dir]->d[Ind_n(i,j,k, n)] = Xb[dir];
  }
}

/* make 2 arrays Cp[0..1] that contain all points of a node, in plane p
   orthogonal to direction dir. The arrays Cp[0..1] are in Xb coords. */
void fill_2arrays_with_nodepoints(tNode *node, int dir, tArray *Cp[2])
{
  int i,j,k, d0,d1, *m0,*m1, ai, c;

  switch(dir)
  {
  case 0:
    d0 = 1;  /* Yb */
    d1 = 2;  /* Zb */
    m0 = &j;
    m1 = &k;
    break;
  case 1:
    d0 = 0;  /* Xb */
    d1 = 2;  /* Zb */
    m0 = &i;
    m1 = &k;
    break;
  case 2:
    d0 = 0;
    d1 = 1;
    m0 = &i;
    m1 = &j;
    break;
  default:
    d0=d1=0; m0=m1=NULL;
    errorexit("dir must be 0,1,2");
  }

  ai = 0;
  forplaneN(dir, i,j,k, node->n, 0)
  {
    double Cb[2] = { node_Xb(node,d0)->d[*m0], node_Xb(node,d1)->d[*m1] };
    for(c=0; c<2; c++) Cp[c]->d[ai] = Cb[c];
    ai++;
  }
}

/* 3d interpolation from array var in node to a set of points given in
   arrays Xp[0..2]. The arrays Xp[0..2] are in Xb coords. The result will
   be written into array interp */
void Lagrange_interpolate_topoints(tNode *node, tArray *var,
                                   tArray *Xp[3], tArray *interp)
{
  int k;
  forarray(Xp[0], k)
  {
    double Xb[]  = { Xp[0]->d[k], Xp[1]->d[k], Xp[2]->d[k] };
    interp->d[k] = Lagrange_array_interpolate(node, var, Xb);
  }
}

/* 3d interpolation from array var in node to a set of points indicated by
   the arrays Xp[0..2] and Ip. Xp[0..2] has the point coords in Xb coords
   and Ip has the index where the interpolation result is written to in
   interp. For points where Ip<0 nothing will be written into interp. */
void Lagrange_interpolate_toIpoints(tNode *node, tArray *var,
                                    tArray *Xp[3], tArray *Ip, tArray *interp)
{
  int k;
  forarray(Xp[0], k)
  {
    double Xb[]  = { Xp[0]->d[k], Xp[1]->d[k], Xp[2]->d[k] };
    int idx = Ip->i[k];
    if(idx>=0)
      interp->d[idx] = Lagrange_array_interpolate(node, var, Xb);
  }
}

/* 2d interpolation from array var in node to a set of points given in
   arrays Cp[0..1], in plane p orthogonal to direction dir. The arrays
   Cp[0..1] are in Xb coords. The result will be written into array interp */
void Lagrange_interpolate2d_topoints(tNode *node, tArray *var, int dir, int p,
                                     tArray *Cp[2], tArray *interp)
{
  int k;
  forarray(Cp[0], k)
  {
    double Cb[]  = { Cp[0]->d[k], Cp[1]->d[k] };
    interp->d[k] = Lagrange_array_interpolate2d(node, var, dir,p, Cb);
  }
}

/* 2d interpolation from array var in node to a set of points indicated by
   the arrays Cp[0..1] and Ip. Cp[0..1] has the point coords in Xb coords
   and Ip has the index where the interpolation result is written to in
   interp. For points where Ip<0 nothing will be written into interp. */
void Lagrange_interpolate2d_toIpoints(tNode *node, tArray *var, int dir,int p,
                                      tArray *Cp[2], tArray *Ip,
                                      tArray *interp)
{
  int k;
  forarray(Cp[0], k)
  {
    double Cb[]  = { Cp[0]->d[k], Cp[1]->d[k] };
    int idx = Ip->i[k];
    if(idx>=0)
      interp->d[k] = Lagrange_array_interpolate2d(node, var, dir,p, Cb);
  }
}


/* insert the values in interp2d into plane p of var */
void insert_array_inplane(tArray *var, int dir, int p, tArray *interp2d)
{
  int i,j,k, *n = var->n;
  int ai = 0;
  forplaneN(dir, i,j,k, n, p)
    var->d[Ind_n(i,j,k, n)] = interp2d->d[ai++];
}
