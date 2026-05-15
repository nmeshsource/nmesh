/* interpolate.c */
/* Wolfgang Tichy, 12/2020
   some functions to deal with interpolation polynomials */

#include "nmesh.h"
#include "basis.h"

/* frequently used global vars */
extern tbasis basis[1];
extern tGridPoints gridpoints[1];


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


/***********************************************************************/
/* Interpolate arrays using a particular interpolation scheme */
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

/* 3d interp of array var onto point Xb using np points around Xb */
double interp_to_Xb0(tElm *elm, tArray *var, double Xb0[3], int np[3],
                     int scheme, double vscal)
{
  tMesh *mesh = Elm_mesh(elm);
  int *nn = elm->n;
  int *pt_typ = elm->pt_typ;
  int schm[] = {scheme,scheme,scheme}; //defaults
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
    int bou;

    /* WENO does currently not work for non-uniform grids */
    if( (scheme==INTERP_WENO) && (pt_typ[d]!=P_UNIFORM) )
      errorexit("INTERP_WENO works only on uniform grids!!!");

    /* reset np from elm->n, if it is zero */
    if(np[d]<1) np[d] = nn[d];

    /* get point coords in elm */
    Xb[d] = node_Xb(elm,d)->d;

    /* get index range into b0[3], nb[3] */
    bou = IndexRange_Xb0_get(elm, d, Xb0[d], np[d], CenterOnXb0,
                             &(b0[d]), &(nb[d]));
    /* get coords in box */
    x_p[d] = Xb[d] + b0[d];

    /* if we want to force Lagrange at boundary */
    if(bou) //if at boundary
      if(Getb(basis->boundary_interp_Lagrange))
        schm[d] = INTERP_LAGRANGE;

    /* get interpolation weights if needed */
    if(schm[d]==INTERP_LAGRANGE)
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
  //for(d=1; d<3; d++)
  //{
  //  printf("d=%d: b0=%d nb=%d ", d, b0[d], nb[d]);
  //  printf("  x_p[%d] =", d);
  //  for(k=0; k<nb[d]; k++) printf(" %g", x_p[d][k]);
  //  printf("\n");
  //}

  /* get pointer vd to start of var data */
  vd = Arrd(var);

  /* interp vd along X for all Y,Z */
  r2 = dmalloc(nb[1]*nb[2]);
  for(k=0; k<nb[2]; k++)
  for(j=0; j<nb[1]; j++)
    r2[j + nb[1]*k] = interpolate1d_ds(Xb0[0], nb[0], x_p[0], schm[0], w[0],
                                       vd + Ind_n(b0[0],b0[1]+j,b0[2]+k, nn),
                                       1, vscal);
  //printf("  r2 =");
  //for(k=0; k<nb[1]*nb[2]; k++) printf(" %g", r2[k]);
  //printf("\n");

  /* interp r2 along Y for all Z */
  r1 = dmalloc(nb[2]);
  for(k=0; k<nb[2]; k++)
    r1[k] = interpolate1d_ds(Xb0[1], nb[1], x_p[1], schm[1], w[1],
                             r2 + nb[1]*k, 1, vscal);
  //printf("  r1 =");
  //for(k=0; k<nb[2]; k++) printf(" %g", r1[k]);
  //printf("\n");

  /* interp r1 along Z */
  interp = interpolate1d_ds(Xb0[2], nb[2], x_p[2], schm[2], w[2],
                            r1, 1, vscal);
  free(r1);
  free(r2);
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
void interpolate_toIpoints(tElm *elm, tArray *var, tArray *Xp[3], tArray *Ip,
                           int np[3], int scheme, double vscal,
                           tArray *interp)
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

/* 3d interpolation from array var in elm onto a set of points given in
   arrays Xp[0..2]. The arrays Xp[0..2] are in Xb coords. The result will
   be written into array interp.
   For INTERP_LAGRANGE this func calls either:
     basis_interp_topoints  OR  interpolate_topoints
   For INTERP_WENO it always uses interpolate_topoints */
void interp_topoints(tElm *elm, tArray *var, tArray *Xp[3],
                     int npts, int scheme, double vscal,
                     tArray *interp)
{
  int *n = elm->n;

  /* if npts is either very small or very large and we want LAGRANGE we use
     the more optimized basis_interp_topoints */
  if( (scheme == INTERP_LAGRANGE) &&
      ( (npts<1) || ((npts>=n[0]) && (npts>=n[1]) && (npts>=n[2])) ) )
  {
    basis_interp_topoints(elm, var, Xp, interp, Lagrange_of_x);
  }
  else
  {
    int np[] = {npts, npts, npts};
    interpolate_topoints(elm, var, Xp, np,scheme,vscal, interp);
  }
}

/* like interp_topoints, but switch the array_MatrixInterp3_scheme
   for certain scheme values */
void interp_topoints_scheme(tElm *elm, tArray *var, tArray *Xp[3],
                            int npts, int scheme, double vscal,
                            tArray *interp)
{
  switch(scheme)
  {
  case INTERP_UNIFORM_TO_n_LGL:
  case INTERP_UNIFORM_TO_nO2_LGL:
    array_MatrixInterp3_scheme(scheme, var, interp);
    break;
  default:
    interp_topoints(elm, var, Xp, npts,scheme,vscal, interp);
  }
}


/* 3d interpolation from array var in node to a set of points indicated by
   the arrays Xp[0..2] and Ip. Xp[0..2] has the point coords in Xb coords
   and Ip has the index where the interpolation result is written to in
   interp. For points where Ip<0 nothing will be written into interp.
   For INTERP_LAGRANGE this func calls either:
     basis_interp_toIpoints  OR  interpolate_toIpoints
   For INTERP_WENO it always uses interpolate_toIpoints */
void interp_toIpoints(tElm *elm, tArray *var, tArray *Xp[3], tArray *Ip,
                      int npts, int scheme, double vscal, tArray *interp)
{
  int *n = elm->n;

  /* if npts is either very small or very large and we want LAGRANGE we use
     the more optimized basis_interp_toIpoints */
  if( (scheme == INTERP_LAGRANGE) &&
      ( (npts<1) || ((npts>=n[0]) && (npts>=n[1]) && (npts>=n[2])) ) )
  {
    basis_interp_toIpoints(elm, var, Xp, Ip, interp, Lagrange_of_x);
  }
  else
  {
    int np[] = {npts, npts, npts};
    interpolate_toIpoints(elm, var, Xp,Ip, np,scheme,vscal, interp);
  }
}


/* 2d interp of array var onto point Xb using np points around Xb.
   interp2d_to_Xb0 is very similar to the 3d interp in interp_to_Xb0 */
double interp2d_to_Xb0(tElm *elm, tArray *var, int dir, int p, double Xb0[2],
                       int np[2], int scheme, double vscal)
{
  tMesh *mesh = Elm_mesh(elm);
  int *nn = elm->n;
  int *an = var->n;
  int *pt_typ = elm->pt_typ;
  int schm[] = {scheme,scheme}; //defaults
  double *Xb[2];
  double *x_p[2];
  double *w[2];
  int CenterOnXb0, b0[2], nb[2];
  double *vd;
  double *r1;
  double interp;
  int d,d2, j,k;

  /* center interp box for WENO only */
  CenterOnXb0 = (scheme==INTERP_WENO);

  d2 = 0;
  for(d=0; d<3; d++)
    if(d != dir)
    {
      int bou;

      /* WENO does currently not work for non-uniform grids */
      if( (scheme==INTERP_WENO) && (pt_typ[d]!=P_UNIFORM) )
        errorexit("INTERP_WENO works only on uniform grids!!!");

      if(an[d] != nn[d])
        errorexiti("elm and var need same dimension in dir%d", d2);

      /* reset np from elm->n, if it is zero */
      if(np[d2]<1) np[d2] = an[d];

      /* get point coords in elm */
      Xb[d2] = node_Xb(elm,d)->d;

      /* get index range into b0[2], nb[2] */
      bou = IndexRange_Xb0_get(elm, d, Xb0[d2], np[d2], CenterOnXb0,
                               &(b0[d2]), &(nb[d2]));
      /* get coords in box */
      x_p[d2] = Xb[d2] + b0[d2];

      /* if we want to force Lagrange at boundary */
      if(bou) //if at boundary
        if(Getb(basis->boundary_interp_Lagrange))
          schm[d2] = INTERP_LAGRANGE;

      /* get interpolation weights if needed */
      if(schm[d2]==INTERP_LAGRANGE)
      {
        w[d2] = dmalloc(nb[d2]);    /* alloc w */
        Lagrange_winterp(nb[d2], x_p[d2], w[d2]);
      }
      else
      {
        w[d2] = NULL;
      }
      /* increment d2 */
      d2++;
    }
  //for(d2=0; d2<2; d2++)
  //{
  //  printf("d2=%d: b0=%d nb=%d ", d2, b0[d2], nb[d2]);
  //  printf("  x_p[%d] =", d2);
  //  for(k=0; k<nb[d2]; k++) printf(" %g", x_p[d2][k]);
  //  printf("\n");
  //}

  /* get pointer vd to start of var data */
  vd = Arrd(var);

  switch(dir)
  {
  case 0:
    /* interp vd along Y for all Z. Note Y has stride an[0]. */
    r1 = dmalloc(nb[1]);
    for(k=0; k<nb[1]; k++)
      r1[k] = interpolate1d_ds(Xb0[0], nb[0], x_p[0], schm[0], w[0],
                               vd + Ind_n(p,b0[0],b0[1]+k, an), an[0], vscal);
    /* interp r1 along Z */
    interp = interpolate1d_ds(Xb0[1], nb[1], x_p[1], schm[1], w[1],
                              r1, 1, vscal);
    free(r1);
    break;
  case 1:
    /* interp vd along X for all Z */
    r1 = dmalloc(nb[1]);
    for(k=0; k<nb[1]; k++)
      r1[k] = interpolate1d_ds(Xb0[0], nb[0], x_p[0], schm[0], w[0],
                               vd + Ind_n(b0[0],p,b0[1]+k, an), 1, vscal);
    /* interp r1 along Z */
    interp = interpolate1d_ds(Xb0[1], nb[1], x_p[1], schm[1], w[1],
                              r1, 1, vscal);
    free(r1);
    break;
  case 2:
    /* interp vd along X for all Y */
    r1 = dmalloc(nb[1]);
    for(j=0; j<nb[1]; j++)
      r1[j] = interpolate1d_ds(Xb0[0], nb[0], x_p[0], schm[0], w[0],
                               vd + Ind_n(b0[0],b0[1]+j,p, an), 1, vscal);
    /* interp r1 along Y */
    interp = interpolate1d_ds(Xb0[1], nb[1], x_p[1], schm[1], w[1],
                              r1, 1, vscal);
    free(r1);
    break;
  default:
    errorexit("dir must be 0,1,2");
  }

  for(d2=1; d2>=0; d2--) free(w[d2]);

  return interp;
}

/* 2d interpolation from array var in elm onto a set of points given in
   arrays Xp[0..1]. The arrays Xp[0..1] are in Xb coords. The result will
   be written into array interp.
   The interpolation will use np[2] points around Xp. scheme describes
   the interpolation scheme, and vscal is the scale used in WENO, usually 1 */
void interpolate2d_topoints(tElm *elm, tArray *var, int dir, int p,
                            tArray *Xp[2], int np[2], int scheme,
                            double vscal, tArray *interp)
{
  int k;
  forarray(Xp[0], k)
  {
    double Xb[]  = { Xp[0]->d[k], Xp[1]->d[k] };
    interp->d[k] = interp2d_to_Xb0(elm, var, dir,p, Xb, np, scheme, vscal);
  }
}

/* 2d interpolation from array var in node to a set of points indicated by
   the arrays Xp[0..1] and Ip. Xp[0..1] has the point coords in Xb coords
   and Ip has the index where the interpolation result is written to in
   interp. For points where Ip<0 nothing will be written into interp. */
void interpolate2d_toIpoints(tElm *elm, tArray *var, int dir, int p,
                             tArray *Xp[2], tArray *Ip, int np[2],
                             int scheme, double vscal, tArray *interp)
{
  int k;
  forarray(Xp[0], k)
  {
    double Xb[]  = { Xp[0]->d[k], Xp[1]->d[k] };
    int idx = Ip->i[k];
    if(idx>=0)
      interp->d[idx] = interp2d_to_Xb0(elm, var, dir,p, Xb, np,scheme,vscal);
  }
}


/* 2d interpolation from array var in elm onto a set of points given in
   arrays Xp[0..1]. The arrays Xp[0..1] are in Xb coords. The result will
   be written into array interp.
   FIXME: For INTERP_LAGRANGE this func should call either:
     basis_interp2d_topoints  OR  interpolate2d_topoints
   For INTERP_WENO it always uses interpolate2d_topoints */
void interp2d_topoints(tElm *elm, tArray *var, int dir, int p,
                       tArray *Xp[2], int npts, int scheme,
                       double vscal, tArray *interp)
{
  int np[] = {npts, npts, npts};
  interpolate2d_topoints(elm, var, dir,p, Xp, np,scheme,vscal, interp);
}

/* 2d interpolation from array var in node to a set of points indicated by
   the arrays Xp[0..1] and Ip. Xp[0..1] has the point coords in Xb coords
   and Ip has the index where the interpolation result is written to in
   interp. For points where Ip<0 nothing will be written into interp.
   FIXME: For INTERP_LAGRANGE this func should call either:
     basis_interp2d_toIpoints  OR  interpolate2d_toIpoints
     (see func interp_toIpoints to see how this is done)
   For INTERP_WENO it always uses interpolate2d_toIpoints */
void interp2d_toIpoints(tElm *elm, tArray *var, int dir, int p,
                        tArray *Xp[2], tArray *Ip, int npts,
                        int scheme, double vscal, tArray *interp)
{
  int np[] = {npts, npts};
  interpolate2d_toIpoints(elm, var, dir,p, Xp,Ip, np,scheme,vscal, interp);
}


/***********************************************************************/
/* Interpolate mesh vars using a particular interpolation scheme */
/***********************************************************************/

/* 3d interpolation from var iu in node to a number of points (set by
   Arrn(interp) = interp->n) with point type pt_typ. The result will be
   written into array interp. */
//once this was:
//void basis_interp_to_pt_typ(tNode *node, int iu, int pt_typ[3], tArray *interp)
void interp_to_pt_typ(tNode *node, int iu, int pt_typ[3],
                      int npts, int scheme, double vscal, tArray *interp)
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
  interp_topoints(node, VarA(node,iu), Xp, npts,scheme,vscal, interp);
  //was: basis_interp_topoints(node, VarA(node,iu), Xp, interp, Lagrange_of_x);

  free_3_arrays(Xp);
}



/* 3d interpolation:
   interpolate var vi to the point (Xb[0],Xb[1],Xb[2]) locally */
double interpolate_var_local(tElm *elm, int vi, double Xb[3],
                             int npts, int scheme, double vscal)
{
  tArray *v;
  int np[] = {npts,npts,npts};
  double val;

  v = VarA(elm, vi);
  if(!v) return 0.; /* return 0 as interp value if var vi has no storage */

  /* interp var to Xb in elm */
  val = interp_to_Xb0(elm, v, Xb, np, scheme, vscal);
  return val;
}

/* same as interpolate_var_local, but more optimized for speed */
double interp_var_local(tElm *elm, int vi, double Xb[3],
                        int npts, int scheme, double vscal)
{
  int *n = elm->n;

  /* if npts is either very small or very large and we want LAGRANGE we use
     the more optimized basis_array_interp */
  if( (scheme == INTERP_LAGRANGE) &&
      ( (npts<1) || ((npts>=n[0]) && (npts>=n[1]) && (npts>=n[2])) ) )
  {
    tArray *v;
    v = VarA(elm, vi);
    if(!v) return 0.; /* return 0 as interp value if var vi has no storage */

    return basis_array_interp(elm, v, Xb, Lagrange_of_x);
  }
  return interpolate_var_local(elm, vi, Xb, npts, scheme, vscal);
}

/* like interp_var_local, but specify point in x,y,z-coords AND interpolate
   for a list of nvars variables
   IN:  elm, nvars,vi, xyz, npts,scheme,vscal, ds<--stride in value-array
   OUT: XYZ, value <-- C-array with values (one for each varindex in vi)
   Returns: pat number of elm, OR <0 if xyz is not in elm */
int interp_vars_xyz_local(tElm *elm, int nvars, int *vi, double xyz[3],
                          int npts, int scheme, double vscal,
                          double XYZ[3], double *value, int ds)
{
  tPat *pat = elm->pat;
  int p, l;
  double XbYbZb[3];

  /* check if xyz is in pat, and set XYZ */
  p = p_XYZ_of_xyz(pat, XYZ, xyz);

  /* if xyz is not on mesh return -1 */
  if(p<0) return -1;

  /* if XYZ is not in elm return a more neg. number */
  if(!XYZ_is_in_node(elm, XYZ)) return -p-2;

  /* set XbYbZb */
  XbYbZb_of_XYZ(elm, XbYbZb, XYZ);

  /* set values for all vars at XYZ */
  for(l=0; l<nvars; l++)
  {
    //double val = interp_var_local(elm, vi[l], XbYbZb, npts, scheme, vscal);
    //GEN_Pragma(omp atomic write)
    //value[l*ds] = val;
    GEN_Pragma(omp atomic write)
    value[l*ds] = interp_var_local(elm, vi[l], XbYbZb, npts, scheme, vscal);
  }

  /* return patch where xyz is in */
  return p;
}
/* same as interp_vars_xyz_local with only 1 var
   IN:  elm, vi, xyz, npts, scheme, vscal     OUT: XYZ, value
   Returns: pat number of elm, OR <0 if xyz is not in elm */
int interp_var_xyz_local(tElm *elm, int vi, double xyz[3],
                         int npts, int scheme, double vscal,
                         double XYZ[3], double *value)
{
  return interp_vars_xyz_local(elm, 1,&vi, xyz, npts,scheme,vscal,
                               XYZ, value,1);
}

/* use interp_vars_xyz_local to interpolate a varlist onto a set of points
   contained in xp[3]
   IN: mesh, vl, xp, npts,scheme,vscal
   OUT: Value <--array with interp values
   value of var at point is here: Arrd(Value)[pt_index + np*vl_index]
   Returns: +1 if all is ok
            -(point index) of 1st point that was not found on mesh */
int interp_VL_xp(tMesh *mesh, tVarList *vl, tArray *xp[3],
                 int npts, int scheme, double vscal, tArray *Value)
{
  int nvars = VLn(vl);
  int *vi   = &(Vind(vl, 0));
  int np = ArrN(xp[0]); /* number of points in xp */
  tArray *value = alloc_array1d(np*nvars);
  double *Val = Arrd(Value);
  double *val = Arrd(value);
  int myrank = nMPI_rank();
  int size = nMPI_size();
  int ret = +1; /* default return value */
  int *rank_pt = imalloc(np);
  int *Rank_pt = imalloc(np);
  int pt;

  if(ArrN(Value) < np*nvars) errorexit("Array named Value is too small!");

  /* init Rank_pt,rank_pt */
  for(pt=0; pt<np; pt++) Rank_pt[pt] = rank_pt[pt] = INT_MAX;

  formyelms(mesh)
  {
    tElm *elm = MyElm;
    int ind;

    for(ind=0; ind<np; ind++)
    {
      double xyz[3], XYZ[3], *valpt;
      int p;

      xyz[0] = Arrd(xp[0])[ind];
      xyz[1] = Arrd(xp[1])[ind];
      xyz[2] = Arrd(xp[2])[ind];
      //PRF;printf(": ind=%d  ", ind);pr3v("xyz", xyz);printf("\n");

      valpt = val + ind; /* pointer to var0 at point ind */
      p = interp_vars_xyz_local(elm, nvars,vi, xyz, npts, scheme, vscal,
                                XYZ, valpt, np);
      /* NOTE: if several threads have the point, we get the interp result from
               whichever thread writes last! */
      /* check if elm has this point */
      if(p>=0)
      {
        /* signal that we have this point */
        GEN_Pragma(omp atomic write)
        rank_pt[ind] = myrank;
      }
    } /* end for i */
  }

  /* copy rank_pt into Rank_pt */
  for(pt=0; pt<np; pt++) Rank_pt[pt] = rank_pt[pt];

  /* get min rank that has point pt into Rank_pt */
  MCK( nMPI_Allreduce(rank_pt, Rank_pt, np, nMPI_INT, nMPI_MIN) );
  /* the val at point pt is later taken from Rank_pt, i.e. we take val from
     lowest rank */

  /* zero val if point pt is on another rank, and copy val into Val */
  for(pt=0; pt<np; pt++)
  {
    int l;
    if((Rank_pt[pt] >= size) && (ret > 0))
    {
      //errorexiti("could not find point %d", pt);
      PRF;printf(": could not find point %d\n", pt);
      ret = -pt;
    }

    /* zero my vals if I don't own them */
    if(Rank_pt[pt] != myrank)
      for(l=0; l<nvars; l++) val[pt + np*l] = 0.;

    /* copy val into Val */
    for(l=0; l<nvars; l++) Val[pt + np*l] = val[pt + np*l];
  }

  /* sum val from all ranks and put result into Val,
     this should give the interp val from the highest rank on each point */
  MCK( nMPI_Allreduce(val, Val, np*nvars, nMPI_DOUBLE, nMPI_SUM) );

  free(Rank_pt);
  free(rank_pt);
  free_array(value);
  return ret;
}


/***********************************************************************/
/* interpolation across MPI procs if based on Xb*/
/***********************************************************************/

/* 3d interpolation:
   call interp_var_local and then send interp. val around
   Returns: 1 if success
            0 if if all MPI procs have node=NULL */
int interpolate_var_ok(tNode *node, int vi, double Xb[3],
                       int npts, int scheme, double *vinterp)
{
  double vscal = 1.; /* set vscal to 1 for now */
  double Val, val=0.;
  int Haveval, haveval=0;

  /* FIXME: use interp_var_xyz instead */
  /* In output0d.c this func is used for xyz interpolation.
     output0d.c should be changed to use interp_var_xyz. */
  //errorexit("this func can only be used for non-critical things "
  //          "like output");

  if(node) if(node->dat)
  {
    val = interp_var_local(node, vi, Xb, npts, scheme, vscal);
    haveval = 1;
  }
  Val = val;
  Haveval = haveval;

  /* find out how many have a value, and add all of them */
  MCK( nMPI_Allreduce(&haveval, &Haveval, 1, nMPI_INT, nMPI_SUM) );
  MCK( nMPI_Allreduce(&val, &Val, 1, nMPI_DOUBLE, nMPI_SUM) );
  if(!Haveval) return 0; /* could not get interp on any node */
  Val = Val/Haveval;

  //PRF;printf(": Val=%g Haveval=%d\n", Val, Haveval);
  *vinterp = Val;
  return 1; /* got interp value */
}

/* 3d interpolation: call interpolate_var_ok */
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
  /* NOTE: node_XYZ_of_xyz_mesh leads to call of p_eid_XYZ_of_xyz_mesh
           BUT node_XYZ_of_xyz_mesh finds only one elm, EVEN IF there are
           several that have x. This is not good! */
  errorexit("this func can only be used for non-critical things "
            "like output");

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
/* interpolate to a given x,y,z using basis_array_interp, i.e.
   using a func: double Basis(...).
   This is not the best choice!!!
   Use interpolate_var_xyz instead. */
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

  /* NOTE: this func is deprecated!!! */
  printf("WARNING: basis_var_interp*xyz funcs are deprecated.\n");
  printf("         Use interp_var_xyz or interp_var_x_y_z instead!\n");

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
      }
  }
  /* if we found points with XYZ, set val to average value */
  if(npts)
    haveval = npts;

  Val = val;
  Haveval = haveval;
  //PRF;printf(": Val=%g Haveval=%d\n", Val, Haveval);

  /* find out how many procs have a value, and add all of them */
  MCK( nMPI_Allreduce(&haveval, &Haveval, 1, nMPI_INT, nMPI_SUM) );
  MCK( nMPI_Allreduce(&val, &Val, 1, nMPI_DOUBLE, nMPI_SUM) );
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
  errorexit("use interp_var_x_y_z instead");
  return basis_var_interp_xyz(mesh, ivar, xyz, Basis);
}


/***********************************************************************/
/* interpolate to a given x,y,z using
   double interp_var_local(tElm *elm, int vi, double Xb[3],
                           int npts, int scheme, double vscal) */
/***********************************************************************/

/* Use interp_var_local to interpolate var ivar onto xyz[3].
   First find node and Xb of xyz[3] and then use use interp_var_local
   to interpolate var ivar onto xyz[3]
   IN: mesh, ivar, xyz, npts, scheme, vscal
   OUT: XYZ, value   RETURN: p (patchnumber) */
int interpolate_var_xyz(tMesh *mesh, int ivar, const double xyz[3],
                        int np, int scheme, double vscal,
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
        val += interp_var_local(elm, ivar, XbYbZb, np, scheme, vscal);
      }
  }
  /* if we found points with XYZ, record how many we found */
  if(npts)
    haveval = npts;

  Val = val;
  Haveval = haveval;
  //PRF;printf(": Val=%g Haveval=%d\n", Val, Haveval);

  /* find out how many values we found, and add all of them */
  MCK( nMPI_Allreduce(&haveval, &Haveval, 1, nMPI_INT, nMPI_SUM) );
  MCK( nMPI_Allreduce(&val, &Val, 1, nMPI_DOUBLE, nMPI_SUM) );
  if(!Haveval) errorexit("at least one MPI proc should have this value");
  Val = Val/Haveval;
  *value = Val;
  //PRF;printf(": Val=%g Haveval=%d\n", Val, Haveval);
  return p;
}

/* interpolate var ivar onto xyz[3] */
double interp_var_xyz(tMesh *mesh, int ivar, const double xyz[3],
                      int np, int scheme, double vscal)
{
  double XYZ[3];
  double value;
  int p = interpolate_var_xyz(mesh, ivar, xyz, np, scheme, vscal, XYZ, &value);
  if(p<0)
  {
    pr3v("xyz",xyz);
    printf("p=%d\n", p);
    errorexit("cannot find xyz on mesh");
  }
  return value;
}

/* same as basis_var_interp_xyz but with differnt interface */
double interp_var_x_y_z(tMesh *mesh, int ivar, double x,double y,double z,
                        int np, int scheme, double vscal)
{
  double xyz[] = {x,y,z};
  return interp_var_xyz(mesh, ivar, xyz, np, scheme, vscal);
}


/***********************************************************************/
/* set interpolation matrices */
/***********************************************************************/

/* write transpose of Lagrange interp. matrix into array Mt */
void Lagrange_InterpMatT(tArray *Xb, tArray *WL, tArray *Yb, tArray *Mt)
{
  if(Xb && Yb)
  {
    int nx = Xb->n[0];
    int ny = Yb->n[0];
    if(nx>0 && ny>0) //do it only if we have grid points
    {
      double *x  = Xb->d;
      double *w  = WL->d;
      double *y  = Yb->d;
      double *MT = Mt->d;
      Lagrange_InterpMatrixT(nx,x, w, ny,y, MT);
    }
  }
}

/* Use transposed Lagrange interp matrix Pt (to interp from LGL to UNIFORM),
   the LGL quad. weights wq, and the UNIFORM quad. weights rq to calculate
   the inverse matrix Rt that allows to transform back from UNIFORM to LGL.
   Both Pt and Rt are the transposes of the interpolation matrices.
   Here we compute Rt using the specified UNIFORM quad. weights rq:
   R = Phi P^T - c Phi (w.w^T) Phi P^T + c Phi (w.r^T)  */
void Inverse_InterpMatT_rq(tArray *Pt, tArray *wq, tArray *rq, tArray *Rt)
{
  if(wq)
  {
    int nu = Pt->n[0]; //dim of u in:     ubar = P u
    int nb = Rt->n[0]; //dim of ubar in:  u    = R ubar
    tArray *Phi   = alloc_array2d(nu, nu);
    tArray *PhiPt = alloc_array2d(nu, nb);
    tArray *tmp_v = alloc_array2d(nu, 1);
    tArray *wwt   = alloc_array2d(nu, nu);
    tArray *R      = alloc_array2d(nu, nb);
    tArray *tmp_R1 = alloc_array2d(nu, nb);
    tArray *tmp_R2 = alloc_array2d(nu, nb);
    double c;

    PRFs(": Pt");printarray_matrix0(Pt);

    /* use tmp_R1 to temporarily store P: tmp_R1 = Pt^T */
    array_swap_dim01(tmp_R1);
    /* set Phi = (P^T P)^{-1}. Note: Phi is symmetric  */
    array_transpose01(Pt, tmp_R1);
    PRFs(": tmp_R1 = P");printarray_matrix0(tmp_R1);
    mm_array_indir(tmp_R1, tmp_R1, 0, Phi); // here tmp_R1 = P
    array_swap_dim01(tmp_R1);
    PRFs(": Phi^{-1}");printarray_matrix0(Phi);
    array_inverse01_inplace(Phi);
    PRFs(": Phi");printarray_matrix0(Phi);

    /* set Phi Pt */
    mm_array_indir(Phi, Pt, 0, PhiPt);
    PRFs(": PhiPt");printarray_matrix0(PhiPt);

    /* set tmp_v = Phi w, and then c = 1/(w^T Phi w) = 1/(wt tmp_v) */
    mm_array_indir(Phi, wq, 0, tmp_v);
    c = 1./array1d_inner_vectorproduct(wq, tmp_v);
    PRFs(": w");printarray_matrix0(wq);
    PRFs(": tmp_v = Phi w");printarray_matrix0(tmp_v);
    printf("c=%g\n", c);

    /* set wwt = w \otimes w^T */
    array1d_outer_vectorproduct(wwt, wq, wq); // wwt = w \otimes w^T
    PRFs(": wwt");printarray_matrix0(wwt);

    /* Set tmp_R1 = (w \otimes w^T) Phi Pt */
    mm_array_indir(wwt, PhiPt, 0, tmp_R1);
    /* Set tmp_R2 = Phi (w \otimes w^T) Phi Pt */
    mm_array_indir(Phi, tmp_R1, 0, tmp_R2);
    PRFs(": tmp_R1 = wwt PhiPt");printarray_matrix0(tmp_R1);
    PRFs(": tmp_R2 = Phi wwt PhiPt");printarray_matrix0(tmp_R2);

    /* Get 1st two terms into R:  R = Phi Pt - c Phi (w \otimes w^T) Phi Pt */
    array_add(R, 1.,PhiPt, -c,tmp_R2);

    /* set tmp_R1 = w \otimes r^T */
    array1d_outer_vectorproduct(tmp_R1, wq, rq); // tmp_R1 = w \otimes r^T
    /* Set tmp_R2 = Phi (w \otimes r^T) */
    mm_array_indir(Phi, tmp_R1, 0, tmp_R2);
    PRFs(": r");printarray_matrix0(rq);
    PRFs(": tmp_R1 = w rt");printarray_matrix0(tmp_R1);
    PRFs(": tmp_R2 = Phi w rt");printarray_matrix0(tmp_R2);

    /* Add last term to R += c Phi (w \otimes r^T) */
    array_addto(R, c,tmp_R2);

    /* finally set Rt = R^T */
    array_transpose01(R, Rt);

    /* check how good the R is. I.e. is  R P = 1  ??? */
    array_transpose01(Pt, tmp_R1);
    mm_array_indir(Rt, tmp_R1, 0, wwt);
    PRFs(": R");printarray_matrix0(R);
    PRFs(": P");printarray_matrix0(tmp_R1);
    PRFs(": RP");printarray_matrix0(wwt);

    /* free temp arrays */
    free_array(tmp_R2);
    free_array(tmp_R1);
    free_array(R);
    free_array(wwt);
    free_array(tmp_v);
    free_array(PhiPt);
    free_array(Phi);
  }
}

/* Use transposed Lagrange interp matrix Pt (to interp from LGL to UNIFORM),
   the LGL quad. weights wq, and the UNIFORM quad. weights rq to calculate
   the inverse matrix Rt that allows to transform back from UNIFORM to LGL.
   Both Pt and Rt are the transposes of the interpolation matrices.
   Here we compute Rt using UNIFORM quad. weights rq coming from
   uniform_x_wGaussquad(ni, Xb, rq);
   Then the formula for R simplifies to: R = Phi P^T */
void Inverse_InterpMatT_best_rq(tArray *Pt, tArray *Rt)
{
  if(Pt)
  {
    int nu = Pt->n[0]; //dim of u in:     ubar = P u
    int nb = Rt->n[0]; //dim of ubar in:  u    = R ubar
    tArray *P   = alloc_array2d(nb, nu);
    tArray *Phi = alloc_array2d(nu, nu);
    tArray *R   = alloc_array2d(nu, nb);

    //PRFs(": Pt");printarray_matrix0(Pt);

    /* set P = P^T  */
    array_transpose01(Pt, P);
    //PRFs(": P");printarray_matrix0(P);

    /* set Phi = (P^T P)^{-1}. Note: Phi is symmetric  */
    mm_array_indir(P, P, 0, Phi); // here tmp_R1 = P
    //PRFs(": Phi^{-1}");printarray_matrix0(Phi);
    array_inverse01_inplace(Phi);
    //PRFs(": Phi");printarray_matrix0(Phi);

    /* set R = Phi Pt */
    mm_array_indir(Phi, Pt, 0, R);
    //PRFs(": R");printarray_matrix0(R);

    /* finally set Rt = R^T */
    array_transpose01(R, Rt);

    /* check how good the R is. I.e. is  R P = 1  ??? */
    //tArray *RP = alloc_array2d(nu, nu);
    //mm_array_indir(Rt, P, 0, RP);
    //PRFs(": P");printarray_matrix0(P);
    //PRFs(": RP");printarray_matrix0(RP);
    //free_array(RP);

    /* free temp arrays */
    free_array(R);
    free_array(Phi);
    free_array(P);
  }
}


/***********************************************************************/
/* functions to interpolate with interpolation matrices */
/***********************************************************************/

/* interp from var into Ivar using the 3 1d interp matrices in Mt[3] */
void array_MatrixInterp3(tArray *Mt[3], tArray *var, tArray *Ivar)
{
  int allocd;

  mm_array_indir(Mt[0], var, 0, Ivar);
  allocd = redim_array(Ivar, Mt[0]->n[1], var->n[1], var->n[2]);
  if(allocd) errorexit("Ivar is to small");

  mm_array_indir(Mt[1], Ivar, 1, Ivar);
  allocd = redim_array(Ivar, -1, Mt[1]->n[1], -1);
  if(allocd) errorexit("Ivar is to small");

  mm_array_indir(Mt[2], Ivar, 2, Ivar);
  allocd = redim_array(Ivar, -1, -1, Mt[2]->n[1]);
  if(allocd) errorexit("Ivar is to small");
}

/* choose inpterp matrix based on scheme, and then interp */
void array_MatrixInterp3_scheme(int scheme, tArray *var, tArray *Ivar)
{
  int n0 = var->n[0];
  int n1 = var->n[1];
  int n2 = var->n[2];

  switch(scheme)
  {
  case INTERP_UNIFORM_TO_n_LGL:
    {
      tArray *Mt[] = { gridpoints->UNI_to_nLGLt[n0],
                       gridpoints->UNI_to_nLGLt[n1],
                       gridpoints->UNI_to_nLGLt[n2] };
      array_MatrixInterp3(Mt, var, Ivar);
    }
    break;
  case INTERP_UNIFORM_TO_nO2_LGL:
    {
      tArray *Mt[] = { gridpoints->UNI_to_no2LGLt[n0],
                       gridpoints->UNI_to_no2LGLt[n1],
                       gridpoints->UNI_to_no2LGLt[n2] };
      tArray *Iv2 = alloc_array(var->n); //temp array of correct size
      array_MatrixInterp3(Mt, var, Iv2);
      copy_array_data(Iv2, Ivar);
      free_array(Iv2);
    }
    break;
  default:
    errorexiti("unknown interpolation scheme %d", scheme);
  }
}
