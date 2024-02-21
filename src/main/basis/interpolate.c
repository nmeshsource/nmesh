/* interpolate.c */
/* Wolfgang Tichy, 12/2020
   some functions to deal with interpolation polynomials */

#include "nmesh.h"
#include "basis.h"


/* frequently used global vars */
extern tbasis basis[1];

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
                          double Basis(int k, double x, int np,
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
  for(k=0; k<n[0]; k++) B0[k] = Basis(k, Xb[0], n[0], xp[0], w[0]);
  for(k=0; k<n[1]; k++) B1[k] = Basis(k, Xb[1], n[1], xp[1], w[1]);
  for(k=0; k<n[2]; k++) B2[k] = Basis(k, Xb[2], n[2], xp[2], w[2]);

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
                            double Basis(int k, double x, int np,
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
    for(k=0; k<n[1]; k++) B1[k] = Basis(k, Cb[0], n[1], xp[1], w[1]);
    for(k=0; k<n[2]; k++) B2[k] = Basis(k, Cb[1], n[2], xp[2], w[2]);

    /* interpolate */
    sum = 0.;
    for(k=0; k<n[2]; k++)
    for(j=0; j<n[1]; j++)
      sum += var->d[Ind_n(p,j,k, n)] * B1[j] * B2[k];
    break;
  case 1:
    /* save basis func values */
    for(k=0; k<n[0]; k++) B1[k] = Basis(k, Cb[0], n[0], xp[0], w[0]);
    for(k=0; k<n[2]; k++) B2[k] = Basis(k, Cb[1], n[2], xp[2], w[2]);

    /* interpolate */
    sum = 0.;
    for(k=0; k<n[2]; k++)
    for(i=0; i<n[0]; i++)
      sum += var->d[Ind_n(i,p,k, n)] * B1[i] * B2[k];
    break;
  case 2:
    /* save basis func values */
    for(k=0; k<n[0]; k++) B1[k] = Basis(k, Cb[0], n[0], xp[0], w[0]);
    for(k=0; k<n[1]; k++) B2[k] = Basis(k, Cb[1], n[1], xp[1], w[1]);

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


/* Take 3 1d arrays and make them into 3 3d arrays. This is useful to
   convert the 3 1d Xb3[3] into real 3d point arrays Xp[3] with coords. */
void array_1d1d1d_coords_to_3d_coords(tArray *X1d[3], tArray *Xp[3])
{
  int i,j,k, d;

  for(d = 0; d < 3; d++)
  {
    double *X1dd = X1d[d]->d;
    double *Xpd = Xp[d]->d;
    int *n = Xp[d]->n;
    int d0 = (d==0);
    int d1 = (d==1);
    int d2 = (d==2);

    for(k = 0; k < n[2]; k++)
    for(j = 0; j < n[1]; j++)
    for(i = 0; i < n[0]; i++)
    {
      int ind = Ind_n(i,j,k, n);
      int l = i*d0 + j*d1 + k*d2;
      Xpd[ind] = X1dd[l];
    }
  }
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
                           double Basis(int k, double x, int np,
                                        const double *x_p,
                                        const double *w_interp))
{
  int k;
  forarray(Xp[0], k)
  {
    double Xb[]  = { Xp[0]->d[k], Xp[1]->d[k], Xp[2]->d[k] };
    interp->d[k] = basis_array_interp(node, var, Xb, Basis);
  }
}

/* 3d interpolation from array var in node to a set of points indicated by
   the arrays Xp[0..2] and Ip. Xp[0..2] has the point coords in Xb coords
   and Ip has the index where the interpolation result is written to in
   interp. For points where Ip<0 nothing will be written into interp. */
void basis_interp_toIpoints(tNode *node, tArray *var,
                            tArray *Xp[3], tArray *Ip, tArray *interp,
                            double Basis(int k, double x, int np,
                                         const double *x_p,
                                         const double *w_interp))
{
  int k;
  forarray(Xp[0], k)
  {
    double Xb[]  = { Xp[0]->d[k], Xp[1]->d[k], Xp[2]->d[k] };
    int idx = Ip->i[k];
    if(idx>=0)
      interp->d[idx] = basis_array_interp(node, var, Xb, Basis);
  }
}

/* 2d interpolation from array var in node to a set of points given in
   arrays Cp[0..1], in plane p orthogonal to direction dir. The arrays
   Cp[0..1] are in Xb coords. The result will be written into array interp */
void basis_interp2d_topoints(tNode *node, tArray *var, int dir, int p,
                             tArray *Cp[2], tArray *interp,
                             double Basis(int k, double x, int np,
                                          const double *x_p,
                                          const double *w_interp))
{
  int k;
  forarray(Cp[0], k)
  {
    double Cb[]  = { Cp[0]->d[k], Cp[1]->d[k] };
    interp->d[k] = basis_array_interp2d(node, var, dir,p, Cb, Basis);
  }
}

/* 2d interpolation from array var in node to a set of points indicated by
   the arrays Cp[0..1] and Ip. Cp[0..1] has the point coords in Xb coords
   and Ip has the index where the interpolation result is written to in
   interp. For points where Ip<0 nothing will be written into interp. */
void basis_interp2d_toIpoints(tNode *node, tArray *var, int dir,int p,
                              tArray *Cp[2], tArray *Ip,
                              tArray *interp,
                              double Basis(int k, double x, int np,
                                           const double *x_p,
                                           const double *w_interp))
{
  int k;
  forarray(Cp[0], k)
  {
    double Cb[]  = { Cp[0]->d[k], Cp[1]->d[k] };
    int idx = Ip->i[k];
    if(idx>=0)
      interp->d[k] = basis_array_interp2d(node, var, dir,p, Cb, Basis);
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


/* 3d interpolation from var iu in node to a number of points (set by
   Arrn(interp) = interp->n) with point type pt_typ. The result will be
   written into array interp. */
void basis_interp_to_pt_typ(tNode *node, int iu, int pt_typ[3],
                            tArray *interp)
{
  int *n_interp = Arrn(interp);
  tArray *Xb[3]; /* 3 1d arrays */
  tArray *Xp[3]; /* 3 3d arrays */

  /* set Xb */
  Xb3_pt_typ_n(pt_typ, n_interp, Xb);

  /* set Xp */
  Xp[0] = alloc_array(n_interp);
  Xp[1] = alloc_array(n_interp);
  Xp[2] = alloc_array(n_interp);
  array_1d1d1d_coords_to_3d_coords(Xb, Xp);

  /* interpolate to points Xp */
  basis_interp_topoints(node, VarA(node,iu), Xp, interp, Lagrange_of_x);

  free_3_arrays(Xp);
}

/***********************************************************************/
/* Interpolate using a particular interpolation scheme */
/***********************************************************************/

/* call Lagrange, or WENO depending on scheme */
double interpolate1d_ds(double x, int n, const double *x_p,
                        int scheme, const double *w_interp,
                        const double *f, int ds, double fscal)
{
  switch(scheme)
  {
  case INTERP_LAGRANGE:
    return Lagrange_interp_barycentric2_ds(x, n,x_p, w_interp, f, ds, fscal);
  case INTERP_WENO:
    return interpolate_WENO_n_ds(x, n,x_p, w_interp, f, ds, fscal);
  default:
    errorexit("unknown scheme");
  }
}

/* 3d interp of var onto point Xb using np points around Xb */
double interp_to_Xb0(tElm *elm, tArray *var, double Xb0[3], int np[3],
                     int scheme, double vscal)
{
  int *nn = elm->n;
  int *pt_typ = elm->pt_typ;
  double *Xb[3];
  double *x_p[3];
  double *w[3];
  int CenterOnXb0, b0[3], nb[3];
  double *vd;
  double *r2;
  double *r1;
  double interp;
  int d, j,k;

  /* center interp box for WENO only */
  CenterOnXb0 = (scheme==INTERP_WENO);

  for(d=0; d<3; d++)
  {
    /* WENO does currently not work for non-uniform grids */
    if( (scheme==INTERP_WENO) && (pt_typ[d]!=P_UNIFORM) )
      errorexit("INTERP_WENO works only on uniform grids!!!");

    /* reset np from elm->n, if it is zero */
    if(np[d]<1) np[d] = nn[d];

    /* get point coords in elm */
    Xb[d] = node_Xb(elm,d)->d;

    /* get index range into b0[3], nb[3] */
    IndexRange_Xb0_get(elm, d, Xb0[d], np[d], CenterOnXb0,
                       &(b0[d]), &(nb[d]));
    /* get coords in box */
    x_p[d] = Xb[d] + b0[d];

    /* get interpolation weights if needed */
    if(scheme==INTERP_LAGRANGE)
    {
      w[d] = dmalloc(nb[d]);    /* alloc w */
      Lagrange_winterp(nb[d], x_p[d], w[d]);
    }
    else
    {
      w[d] = NULL;
    }
    /*
    PRF;printf(": d=%d  b0[d]=%d nb[d]=%d\n", d, b0[d], nb[d]);
    printf("  x_p[d] =");
    for(k=0; k<nb[d]; k++) printf(" %g", x_p[d][k]);
    printf("\n");
    if(w[d])
    {
      printf("  w[d]   =");
      for(k=0; k<nb[d]; k++) printf(" %g", w[d][k]);
      printf("\n");
      //printf("  WL[d]  =");
      //for(k=0; k<nb[d]; k++) printf(" %g", node_WL(elm,d)->d[k]);
      //printf("\n");
    }
    */
  }
  //d=0;
  //printf("  x_p[d] =");
  //for(k=0; k<nb[d]; k++) printf(" %g", x_p[d][k]);
  //printf("\n");

  /* get pointer vd to start of var data */
  vd = Arrd(var);

  /* interp vd along X for all Y,Z */
  r2 = dmalloc(nb[1]*nb[2]);
  for(k=0; k<nb[2]; k++)
  for(j=0; j<nb[1]; j++)
    r2[j + nb[1]*k] = interpolate1d_ds(Xb0[0], nb[0], x_p[0], scheme, w[0],
                                       vd + Ind_n(b0[0],b0[1]+j,b0[2]+k, nn),
                                       1, vscal);
  //printf("  r2 =");
  //for(k=0; k<nb[1]*nb[2]; k++) printf(" %g", r2[k]);
  //printf("\n");

  /* interp r2 along Y for all Z */
  r1 = dmalloc(nb[2]);
  for(k=0; k<nb[2]; k++)
    r1[k] = interpolate1d_ds(Xb0[1], nb[1], x_p[1], scheme, w[1],
                             r2 + nb[1]*k, 1, vscal);
  //printf("  r1 =");
  //for(k=0; k<nb[2]; k++) printf(" %g", r1[k]);
  //printf("\n");

  /* interp r1 along Z */
  interp = interpolate1d_ds(Xb0[2], nb[2], x_p[2], scheme, w[2],
                            r1, 1, vscal);
  free(r2);
  free(r1);
  for(d=2; d>=0; d--) free(w[d]);

  return interp;
}

/* 3d interpolation from array var in elm onto a set of points given in
   arrays Xp[0..2]. The arrays Xp[0..2] are in Xb coords. The result will
   be written into array interp.
   The interpolation will use np[3] points around Xp. scheme describes
   the interpolation scheme, and vscal is the scale use in WENO, usually 1 */
void interpolate_topoints(tElm *elm, tArray *var, tArray *Xp[3],
                          int np[3], int scheme, double vscal,
                          tArray *interp)
{
  int k;
  forarray(Xp[0], k)
  {
    double Xb[]  = { Xp[0]->d[k], Xp[1]->d[k], Xp[2]->d[k] };
    interp->d[k] = interp_to_Xb0(elm, var, Xb, np, scheme, vscal);
  }
}

/* 3d interpolation from array var in node to a set of points indicated by
   the arrays Xp[0..2] and Ip. Xp[0..2] has the point coords in Xb coords
   and Ip has the index where the interpolation result is written to in
   interp. For points where Ip<0 nothing will be written into interp. */
void interpolate_toIpoints(tElm *elm, tArray *var, tArray *Xp[3],
                           int np[3], int scheme, double vscal,
                           tArray *Ip, tArray *interp)
{
  int k;
  forarray(Xp[0], k)
  {
    double Xb[]  = { Xp[0]->d[k], Xp[1]->d[k], Xp[2]->d[k] };
    int idx = Ip->i[k];
    if(idx>=0)
      interp->d[idx] = interp_to_Xb0(elm, var, Xb, np, scheme, vscal);
  }
}

/* Return scheme and set np and vscal from basis pars.
   scheme_pref is a preference that can be overridden on some grids */
int interpolate_scheme_get(tElm *elm, int vi,
                           int np[3], int scheme_pref, double *vscal)
{
  tMesh *mesh = Elm_mesh(elm);
  int *pt_typ = elm->pt_typ;
  int is_UNI, scheme, npts;
  int d;

  errorexit("this func is not tested!");

  /* set scheme to preferred value */
  scheme = scheme_pref;

  /* check grid */
  is_UNI=1;
  for(d=0; d<3; d++)
    if(pt_typ[d]!=P_UNIFORM) { is_UNI=0; break; }

  /* switch from INTERP_WENO tp INTERP_LAGRANGE on non-unif. grids */
  if( (scheme==INTERP_WENO) && (!is_UNI) ) scheme = INTERP_LAGRANGE;;

  /* set np */
  switch(scheme)
  {
  case INTERP_LAGRANGE:
    npts = Geti(basis->Lagrange_interp_np);
    for(d=0; d<3; d++)
    {
      if(npts<1) np[d] = elm->n[d];
      else       np[d] = npts;
    }
    break;
  case INTERP_WENO:
    npts = Geti(basis->WENO_interp_np);
    for(d=0; d<3; d++) np[d] = npts;
    break;
  default:
    errorexit("unknown scheme");
  }

  /* set vscal to 1 for now */
  *vscal = 1.;
  return scheme;
}


/***********************************************************************/
/* interpolation across MPI procs if based on Xb*/
/***********************************************************************/

/* 3d interpolation:
   interpolate var vi to the point (Xb[0],Xb[1],Xb[2]) locally */
double interpolate_var_local(tElm *elm, int vi, double Xb[3],
                             int npts, int scheme)
{
  tArray *v;
  int np[] = {npts,npts,npts};
  double vscal = 1.; /* set vscal to 1 for now */
  double val;

  v = VarA(elm, vi);
  if(!v) return 0.; /* return 0 as interp value if var vi has no storage */

  /* interp var to Xb in elm */
  val = interp_to_Xb0(elm, v, Xb, np, scheme, vscal);
  return val;
}

/* 3d interpolation:
   call basis_var_interpolate_local and then send interp. val around
   Returns: 1 if success
            0 if if all MPI procs have node=NULL */
int interpolate_var_ok(tNode *node, int vi, double Xb[3],
                       int npts, int scheme, double *vinterp)
{
  double Val, val=0.;
  int Haveval, haveval=0;

  if(node) if(node->dat)
  {
    val = interpolate_var_local(node, vi, Xb, npts, scheme);
    haveval = 1;
  }
  Val = val;
  Haveval = haveval;

  /* find out how many have a value, and add all of them */
  nMPI_Allreduce(&haveval, &Haveval, 1, nMPI_INT, nMPI_SUM);
  nMPI_Allreduce(&val, &Val, 1, nMPI_DOUBLE, nMPI_SUM);
  if(!Haveval) return 0; /* could not get interp on any node */
  Val = Val/Haveval;

  //PRF;printf(": Val=%g Haveval=%d\n", Val, Haveval);
  *vinterp = Val;
  return 1; /* got interp value */
}

/* 3d interpolation: call basis_var_interpolate_ok */
double interpolate_var(tNode *node, int vi, double Xb[3],
                       int npts, int scheme)
{
  double Val;
  int Haveval = interpolate_var_ok(node, vi, Xb, npts, scheme, &Val);
  if(!Haveval) errorexit("one MPI proc should have this node");

  //PRF;printf(": Val=%g Haveval=%d\n", Val, Haveval);
  return Val;
}

/* 3d interpolation:
   interpolate var vi to the point (x[0],x[1],x[2]).
   out: val
   returns: 1 if success, or 0 if failure to find x */
int interpolate_var_mesh(tMesh *mesh, int vi, const double x[3],
                         int npts, int scheme, double *val)
{
  int Haveval;
  double X[3], Xb[3];
  tNode *node = node_XYZ_of_xyz_mesh(mesh, X, x);

  /* set Xb in node */
  if(node) XbYbZb_of_XYZ(node, Xb, X);

  /* interp var vi to Xb in node */
  Haveval = interpolate_var_ok(node, vi, Xb, npts, scheme, val);
  return Haveval;
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


/***********************************************************************/
/* interpolate to a given x,y,z */
/***********************************************************************/

/* use basis_array_interp to interpolate var ivar onto xyz[3],
   IN: mesh, ivar, xyz, basis   OUT: XYZ, value   RETURN: p (patchnumber) */
int basis_var_interpolate_xyz(tMesh *mesh, int ivar, const double xyz[3],
                              double Basis(int k, double x, int np,
                                           const double *x_p,
                                           const double *w_interp),
                              double XYZ[3], double *value)
{
  double Val, val;
  int Haveval, haveval;
  int p, npts;

  /* find patch p and set XYZ */
  p = p_XYZ_of_xyz_mesh(mesh, XYZ, xyz);
  //PRFs(": ");pr3v("XYZ", XYZ);printf(": p=%d\n", p);

  /* if xyz is not on mesh return -1 */
  if(p<0) return p;

  /* search among my leaf nodes (in patch p) for XYZ */
  val = 0.;
  npts = 0;
  haveval = 0;
  formyelms_noomp(mesh)
  {
    tElm *elm = MyElm;
    if(elm->pat->p == p)
      if(XYZ_is_in_node(elm, XYZ))
      {
        double XbYbZb[3];

        npts++; /* count how often we found XYZ */
        XbYbZb_of_XYZ(elm, XbYbZb, XYZ);
        val += basis_array_interp(elm, VarA(elm, ivar), XbYbZb, Basis);
        break; //FIXME: remove this!!!!
      }
  }
  /* if we found points with XYZ, set val to average value */
  if(npts)
  {
    val = val/npts; /* average val */
    haveval = 1;
  }

  Val = val;
  Haveval = haveval;
  //PRF;printf(": Val=%g Haveval=%d\n", Val, Haveval);

  /* find out how many procs have a value, and add all of them */
  nMPI_Allreduce(&haveval, &Haveval, 1, nMPI_INT, nMPI_SUM);
  nMPI_Allreduce(&val, &Val, 1, nMPI_DOUBLE, nMPI_SUM);
  if(!Haveval) errorexit("one MPI proc should have this value");
  Val = Val/Haveval;
  *value = Val;
  //PRF;printf(": Val=%g Haveval=%d\n", Val, Haveval);
  return p;
}

/* find node and Xb of xyz[3] and then use use basis_array_interp
   to interpolate var ivar onto xyz[3] */
double basis_var_interp_xyz(tMesh *mesh, int ivar, double xyz[3],
                            double Basis(int k, double x, int np,
                                         const double *x_p,
                                         const double *w_interp))
{
  double XYZ[3];
  double value;
  int p = basis_var_interpolate_xyz(mesh, ivar, xyz, Basis, XYZ, &value);
  if(p<0)
  {
    pr3v("xyz",xyz);
    printf("p=%d\n", p);
    errorexit("cannot find xyz on mesh");
  }
  return value;
}

/* same as basis_var_interp_xyz but with differnt interface */
double basis_var_interp_x_y_z(tMesh *mesh, int ivar,
                              double x,double y,double z,
                              double Basis(int k, double x, int np,
                                         const double *x_p,
                                         const double *w_interp))
{
  double xyz[] = {x,y,z};
  return basis_var_interp_xyz(mesh, ivar, xyz, Basis);
}
