/* interpolate.c */
/* Wolfgang Tichy, 12/2020
   some functions to deal with interpolation polynomials */

#include "nmesh.h"
#include "basis.h"



/* ************************************************************************ */
/* various functions needed for piecewise const or linear interpolation     */
/* ************************************************************************ */

/* piecewise const basis functions */
double basis_pw_const(int k, double x, int np,
                      const double *x_p, const double *w_interp)
{
  int n;
  double xml, xmr;

  /* special case for just 1 point */
  if(np<=1) return 1.;

  if(k<=0) /* lowest k */
  {
    n = 0;
    xmr = 0.5*(x_p[n] + x_p[n+1]);
    xml = -DBL_MAX;
  }
  else if(k<np-1) /* k in the middle */
  {
    n = k;
    xml = 0.5*(x_p[k] + x_p[k-1]);
    xmr = 0.5*(x_p[k] + x_p[k+1]);
  }
  else /* highest k */
  {
    n = np-1;
    xml = 0.5*(x_p[n] + x_p[n-1]);
    xmr = DBL_MAX;
  }

  if(x >= xml && x < xmr)
    return 1.;
  else
    return 0.;
}

/* piecewise linear basis functions */
double basis_pw_linear(int k, double x, int np,
                       const double *x_p, const double *w_interp)
{
  /* special case for just 1 point: do same as in basis_pw_const */
  if(np<=1) return 1.;

  /* x is left */
  if(x < x_p[1])
  {
    switch(k)
    {
    case 0:
      return (x - x_p[1])/(x_p[0] - x_p[1]);
    case 1:
      return (x - x_p[0])/(x_p[1] - x_p[0]);
    default:
      return 0.;
    }
  }

  /* x is right */
  if(x >= x_p[np-2])
  {
    int n0 = np-2;
    int n1 = np-1;
    int l = k - n0;
    switch(l)
    {
    case 0:
      return (x - x_p[n1])/(x_p[n0] - x_p[n1]);
    case 1:
      return (x - x_p[n0])/(x_p[n1] - x_p[n0]);
    default:
      return 0.;
    }
  }

  /* x in middle: */
  if(k < np-1)
  {
    if(x >= x_p[k] && x < x_p[k+1])
      return (x - x_p[k+1])/(x_p[k] - x_p[k+1]);
  }
  if(k > 0)
  {
    if(x >= x_p[k-1] && x < x_p[k])
      return (x - x_p[k-1]) / (x_p[k] - x_p[k-1]);
  }

  return 0.;
}

/* piecewise parabolic basis functions */
double basis_pw_parab(int k, double x, int np,
                      const double *x_p, const double *w_interp)
{
  int m, m1, m2;
  double xml, xmr, xmll, xmrr;

  /* special case for less than 3 points */
  if(np<3) return basis_pw_linear(k, x, np, x_p, w_interp);

  /* special case for only 3 points*/
  if(np==3)
  {
    switch(k)
    {
    case 0:
      return  (     x - x_p[1])*(     x - x_p[2]) /
             ((x_p[0] - x_p[1])*(x_p[0] - x_p[2]));
    case 1:
      return  (     x - x_p[0])*(     x - x_p[2]) /
             ((x_p[1] - x_p[0])*(x_p[1] - x_p[2]));
    case 2:
      return  (     x - x_p[0])*(     x - x_p[1]) /
             ((x_p[2] - x_p[0])*(x_p[2] - x_p[1]));
    default:
      return 0.;
    }
  }

  /* x is left */
  xmr = 0.5*(x_p[0] + x_p[1]);
  if(x < xmr)
  {
    switch(k)
    {
    case 0:
      return  (     x - x_p[1])*(     x - x_p[2]) /
             ((x_p[0] - x_p[1])*(x_p[0] - x_p[2]));
    case 1:
      return  (     x - x_p[0])*(     x - x_p[2]) /
             ((x_p[1] - x_p[0])*(x_p[1] - x_p[2]));
    case 2:
      return  (     x - x_p[0])*(     x - x_p[1]) /
             ((x_p[2] - x_p[0])*(x_p[2] - x_p[1]));
    default:
      return 0.;
    }
  }

  /* x is right */
  xml = 0.5*(x_p[np-2] + x_p[np-1]);
  if(x >= xml)
  {
    int n0 = np-3;
    int n1 = np-2;
    int n2 = np-1;
    int l = k - n0;

    switch(l)
    {
    case 0:
      return  (      x - x_p[n1])*(      x - x_p[n2]) /
             ((x_p[n0] - x_p[n1])*(x_p[n0] - x_p[n2]));
    case 1:
      return  (      x - x_p[n0])*(      x - x_p[n2]) /
             ((x_p[n1] - x_p[n0])*(x_p[n1] - x_p[n2]));
    case 2:
      return  (      x - x_p[n0])*(      x - x_p[n1]) /
             ((x_p[n2] - x_p[n0])*(x_p[n2] - x_p[n1]));
    default:
      return 0.;
    }
  }

  /* x is in middle: */
  if(k<=0) /* lowest k */
  {
    m = 0;
    xmr = 0.5*(x_p[m] + x_p[m+1]);
    xml = x_p[m] - (xmr - x_p[m]);
    xmll = xml - (xmr - xml);
    xmrr = 0.5*(x_p[m+1] + x_p[m+2]);
    m1 = 1;
    m2 = 2;
  }
  else if(k<np-1) /* k in the middle */
  {
    m = k;
    xml = 0.5*(x_p[m] + x_p[m-1]);
    xmr = 0.5*(x_p[m] + x_p[m+1]);
    if(k<np-2) xmrr = 0.5*(x_p[m+1] + x_p[m+2]);
    else       xmrr = xmr + (xmr - xml);
    if(k>=2) xmll = 0.5*(x_p[m-1] + x_p[m-2]);
    else     xmll = xml - (xmr - xml);
    m1 = k-1;
    m2 = k+1;
  }
  else /* highest k */
  {
    m = np-1;
    xml = 0.5*(x_p[m-1] + x_p[m]);
    xmr = x_p[m] + (x_p[m] - xml);
    xmll = 0.5*(x_p[m-2] + x_p[m-1]);
    xmrr = xmr + (xmr - xml);
    m1 = m-2;
    m2 = m-1;
  }

  if(x >= xml && x < xmr)
    return  (     x - x_p[m1])*(     x - x_p[m2]) /
           ((x_p[m] - x_p[m1])*(x_p[m] - x_p[m2]));

  if(x >= xmll && x < xml)
    return  (     x - x_p[m-2])*(     x - x_p[m-1]) /
           ((x_p[m] - x_p[m-2])*(x_p[m] - x_p[m-1]));

  if(x >= xmr && x < xmrr)
    return  (     x - x_p[m+1])*(     x - x_p[m+2]) /
           ((x_p[m] - x_p[m+1])*(x_p[m] - x_p[m+2]));

  return 0.;
}


/***********************************************************************/
/* interpolate using some basis */
/***********************************************************************/

/* 3d interpolation:
   interpolate to the point (Xb[0],Xb[1],Xb[2]) for variable in array var
   using the basis function basis we pass in.
   Note: for Lagrange interpolation the coeffs are simply the function
         values at grid points
   Note2: the only info we really retrieve from the node (in Xb3_n and WL3_n)
          is node->pt_typ, but not node->n */
double basis_array_interp(tNode *node, tArray *var, double Xb[3],
                          double basis(int k, double x, int np,
                                       const double *x_p,
                                       const double *w_interp))
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
  for(k=0; k<n[0]; k++) B0[k] = basis(k, Xb[0], n[0], xp[0], w[0]);
  for(k=0; k<n[1]; k++) B1[k] = basis(k, Xb[1], n[1], xp[1], w[1]);
  for(k=0; k<n[2]; k++) B2[k] = basis(k, Xb[2], n[2], xp[2], w[2]);

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
double basis_array_interp2d(tNode *node, tArray *var, int dir, int p,
                            double Cb[2],
                            double basis(int k, double x, int np,
                                         const double *x_p,
                                         const double *w_interp))
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
    for(k=0; k<n[1]; k++) B1[k] = basis(k, Cb[0], n[1], xp[1], w[1]);
    for(k=0; k<n[2]; k++) B2[k] = basis(k, Cb[1], n[2], xp[2], w[2]);

    /* interpolate */
    sum = 0.;
    for(k=0; k<n[2]; k++)
    for(j=0; j<n[1]; j++)
      sum += var->d[Ind_n(p,j,k, n)] * B1[j] * B2[k];
    break;
  case 1:
    /* save basis func values */
    for(k=0; k<n[0]; k++) B1[k] = basis(k, Cb[0], n[0], xp[0], w[0]);
    for(k=0; k<n[2]; k++) B2[k] = basis(k, Cb[1], n[2], xp[2], w[2]);

    /* interpolate */
    sum = 0.;
    for(k=0; k<n[2]; k++)
    for(i=0; i<n[0]; i++)
      sum += var->d[Ind_n(i,p,k, n)] * B1[i] * B2[k];
    break;
  case 2:
    /* save basis func values */
    for(k=0; k<n[0]; k++) B1[k] = basis(k, Cb[0], n[0], xp[0], w[0]);
    for(k=0; k<n[1]; k++) B2[k] = basis(k, Cb[1], n[1], xp[1], w[1]);

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
void basis_interp_topoints(tNode *node, tArray *var,
                           tArray *Xp[3], tArray *interp,
                           double basis(int k, double x, int np,
                                        const double *x_p,
                                        const double *w_interp))
{
  int k;
  forarray(Xp[0], k)
  {
    double Xb[]  = { Xp[0]->d[k], Xp[1]->d[k], Xp[2]->d[k] };
    interp->d[k] = basis_array_interp(node, var, Xb, basis);
  }
}

/* 3d interpolation from array var in node to a set of points indicated by
   the arrays Xp[0..2] and Ip. Xp[0..2] has the point coords in Xb coords
   and Ip has the index where the interpolation result is written to in
   interp. For points where Ip<0 nothing will be written into interp. */
void basis_interp_toIpoints(tNode *node, tArray *var,
                            tArray *Xp[3], tArray *Ip, tArray *interp,
                            double basis(int k, double x, int np,
                                         const double *x_p,
                                         const double *w_interp))
{
  int k;
  forarray(Xp[0], k)
  {
    double Xb[]  = { Xp[0]->d[k], Xp[1]->d[k], Xp[2]->d[k] };
    int idx = Ip->i[k];
    if(idx>=0)
      interp->d[idx] = basis_array_interp(node, var, Xb, basis);
  }
}

/* 2d interpolation from array var in node to a set of points given in
   arrays Cp[0..1], in plane p orthogonal to direction dir. The arrays
   Cp[0..1] are in Xb coords. The result will be written into array interp */
void basis_interp2d_topoints(tNode *node, tArray *var, int dir, int p,
                             tArray *Cp[2], tArray *interp,
                             double basis(int k, double x, int np,
                                          const double *x_p,
                                          const double *w_interp))
{
  int k;
  forarray(Cp[0], k)
  {
    double Cb[]  = { Cp[0]->d[k], Cp[1]->d[k] };
    interp->d[k] = basis_array_interp2d(node, var, dir,p, Cb, basis);
  }
}

/* 2d interpolation from array var in node to a set of points indicated by
   the arrays Cp[0..1] and Ip. Cp[0..1] has the point coords in Xb coords
   and Ip has the index where the interpolation result is written to in
   interp. For points where Ip<0 nothing will be written into interp. */
void basis_interp2d_toIpoints(tNode *node, tArray *var, int dir,int p,
                              tArray *Cp[2], tArray *Ip,
                              tArray *interp,
                              double basis(int k, double x, int np,
                                           const double *x_p,
                                           const double *w_interp))
{
  int k;
  forarray(Cp[0], k)
  {
    double Cb[]  = { Cp[0]->d[k], Cp[1]->d[k] };
    int idx = Ip->i[k];
    if(idx>=0)
      interp->d[k] = basis_array_interp2d(node, var, dir,p, Cb, basis);
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


/***********************************************************************/
/* Lagrange interpolation */
/***********************************************************************/

/* 3d interpolation:
   interpolate to the point (Xb[0],Xb[1],Xb[2]) for field in array vals
   the pts[3] arrays contain the point coords in the 3 dirs of vals.
   vals and pts[3] are set with extract_vals_pts_around_Xb */
double Lagrange_array_interp(tArray *vals, double *pts[3], double Xb[3])
{
  int *n = vals->n;
  double *w[3];  /* weights */
  double *B[3];  /* basis */
  int d, k;
  double sum;

errorexit("Lagrange_array_interp is untested!!!");

  /* get arrays for basis and weights for n = vals->n */
  for(d=0; d<3; d++)
  {
    w[d] = dmalloc(n[d]);
    B[d] = dmalloc(n[d]);
  }

  /* set weights */
  for(d=0; d<3; d++) Lagrange_winterp(n[d], pts[d], w[d]);

  /* save basis func values at (Xb[0],Xb[1],Xb[2]) in B[0],... */
  for(d=0; d<3; d++)
    for(k=0; k<n[0]; k++)
      B[d][k] = Lagrange_of_x(k, Xb[d], n[d], pts[d], w[d]);

  /* interpolate to (Xb[0],Xb[1],Xb[2]) */
  sum = 0.;
  //SGRID_LEVEL3_Pragma(omp parallel for reduction(+:sum))
  for(k=0; k<n[2]; k++)
  {
    int j,i;
    for(j=0; j<n[1]; j++)
    for(i=0; i<n[0]; i++)
      sum += vals->d[Ind_n(i,j,k, n)] * B[0][i] * B[1][j] * B[2][k];
  }

  /* free saved stuff */
  for(d=2; d>=0; d--)
  {
    free(B[d]);
    free(w[d]);
  }
  return sum;
}

/* Note1: the size of how much we extract from var depends on dimensions of
   vals array
   Note2: the only info we really retrieve from the node (in Xb3_n)
          is node->pt_typ, but not node->n */
void extract_vals_pts_around_Xb(tNode *node, tArray *var, double Xb[3],
                                tArray *vals, double *pts[3])
{
  int *n = var->n;
//  int *ne = vals->n;
  tArray *Xb_n[3];   /* points */
  int d, k;

errorexit("extract_vals_pts_around_Xb is unfinished!!!");

  /* get arrays with points and weights for n = var->n */
  Xb3_n(node, n, Xb_n);

  /* find 2 closest points in each dir */
  for(d=0; d<3; d++)
  {
    double d1=DBL_MAX, d2=DBL_MAX;
    int k1=-1, k2=-1;

    for(k=0; k<n[d]; k++)
    {
      double x = Xb_n[d]->d[k];
      double dist = fabs(x - Xb[d]);
      if(dist < d1)
      {
        d2 = d1;
        k2 = k1;
        d1 = dist;
        k1 = k;
      }
      else if(dist < d2)
      {
        d2 = dist;
        k2 = k;
      }
    }
    if(k2<0) exit(9);
  }
  //... set points in pts and then vals

}
