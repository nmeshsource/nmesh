/* rec1d.c */
/* Wolfgang Tichy, October 2020 */

/* several reconstruction methods in 1d */

#include "nmesh.h"
#include "dg.h"




/* Interpolate a field u to midpoint with index im.
   Here we interpolate in the positive direction (p) from the left of the
   midpoint to the midpoint to obtain umid_p. */
double rec1d_p_1(int n, const double *u, int im, double u_scale)
{
  return u[im]; // one sided 0-th order interpolation
}

/* Interpolate a field u to midpoint with index im.
   Here we interpolate in the negative direction (m) from the right of the
   midpoint to the midpoint to obtain umid_m. */
double rec1d_m_1(int n, const double *u, int im, double u_scale)
{
  return u[im+1]; // one sided 0-th order interpolation
}



/* ideal weights*3 and epsilon for WENO3,
   taken from https://math.la.asu.edu/~gardner/weno.pdf */
#define WENO3_3id_gamma1 1.
#define WENO3_3id_gamma2 2.
#define WENO3_epsilon 1e-6

/* Interpolate a field u to midpoint at i+1/2 with index im.
   Here we interpolate in the positive direction (p) from the left of the
   midpoint to the midpoint to obtain umid_p.
   We use u at the points i-1, i, i+1 */
double rec1d_p_WENO3_uniform__old(int n, const double *u, int im, double u_scale)
{
  /* u at 3 grid points around point i=im */
  double u_im1 = u[im-1];
  double u_i   = u[im];
  double u_ip1 = u[im+1];

  /* diffs */
  double d1 = u_i - u_im1;
  double d2 = u_ip1 - u_i;

  /* smoothness indicators */
  double beta1 = d1*d1;
  double beta2 = d2*d2;
  double us2 = u_scale*u_scale;
  double beta1_p_eps = beta1 + WENO3_epsilon*us2;
  double beta2_p_eps = beta2 + WENO3_epsilon*us2;

  /* non-normalized weights */
  double omegab1 = WENO3_3id_gamma1/(beta1_p_eps*beta1_p_eps);
  double omegab2 = WENO3_3id_gamma2/(beta2_p_eps*beta2_p_eps);
  double omegab_sum = omegab1 + omegab2;

  /* normalized weights */
  double omega1 = omegab1/(omegab_sum);
  double omega2 = omegab2/(omegab_sum);

  /* 2 linear reconstructions */
  double u1 = -0.5*u_im1 + 1.5*u_i;
  double u2 = 0.5*(u_i + u_ip1);

  /* final reconstruction */
  return omega1*u1 + omega2*u2;
}

/* Interpolate a field u to midpoint i+1/2 with index im.
   Here we interpolate in the negative direction (m) from the right of the
   midpoint to the midpoint to obtain umid_m.
   We use u at the points i, i+1, i+2 */
double rec1d_m_WENO3_uniform__old(int n, const double *u, int im, double u_scale)
{
  /* u at 3 grid points around point i=im */
  double u_i   = u[im];
  double u_ip1 = u[im+1];
  double u_ip2 = u[im+2];

  /* diffs */
  double d1 = u_ip1 - u_i;
  double d2 = u_ip2 - u_ip1;

  /* smoothness indicators */
  double beta1 = d1*d1;
  double beta2 = d2*d2;
  double us2 = u_scale*u_scale;
  double beta1_p_eps = beta1 + WENO3_epsilon*us2;
  double beta2_p_eps = beta2 + WENO3_epsilon*us2;

  /* non-normalized weights */
  double omegab1 = WENO3_3id_gamma2/(beta1_p_eps*beta1_p_eps);
  double omegab2 = WENO3_3id_gamma1/(beta2_p_eps*beta2_p_eps);
  double omegab_sum = omegab1 + omegab2;

  /* normalized weights */
  double omega1 = omegab1/(omegab_sum);
  double omega2 = omegab2/(omegab_sum);

  /* 2 linear reconstructions */
  double u1 = 0.5*(u_i + u_ip1);
  double u2 = 1.5*u_ip1 - 0.5*u_ip2;

  /* final reconstruction */
  return omega1*u1 + omega2*u2;
}


/* Interpolate a field u to midpoint at i+1/2 with index im.
   Here we interpolate in the positive direction (p) from the left of the
   midpoint to the midpoint to obtain umid_p.
   We use u at the points i-1, i, i+1 */
/* more general WENO that uses weights in global WENOweights var */
double rec1d_p_WENO3(int n, const double *u, int im, double u_scale,
                     tWENO3weight *W3)
{
  /* u at 3 grid points around point i=im */
  double u_im1 = u[im-1];
  double u_i   = u[im];
  double u_ip1 = u[im+1];

  /* diffs */
  double d1 = u_i - u_im1; // maybe we should divide here by grid spacing???
  double d2 = u_ip1 - u_i;

  /* smoothness indicators */
  double beta1 = d1*d1;
  double beta2 = d2*d2;
  double us2 = u_scale*u_scale;
  double beta1_p_eps = beta1 + WENO3_epsilon*us2;
  double beta2_p_eps = beta2 + WENO3_epsilon*us2;

  /* non-normalized weights */
  double omegab1 = W3->optw[0]/(beta1_p_eps*beta1_p_eps);
  double omegab2 = W3->optw[1]/(beta2_p_eps*beta2_p_eps);
  double omegab_sum = omegab1 + omegab2;

  /* normalized weights */
  double omega1 = omegab1/(omegab_sum);
  double omega2 = omegab2/(omegab_sum);

  /* 2 linear reconstructions */
  double u1 = W3->lw[0][0]*u_im1 + W3->lw[0][1]*u_i;
  double u2 = W3->lw[1][0]*u_i   + W3->lw[1][1]*u_ip1;

  /* final reconstruction */
  return omega1*u1 + omega2*u2;
}

/* more general WENO that uses weights in global WENOweights var */
double rec1d_p_WENO3_at_im(int n, const double *u, int im, double u_scale)
{
  tWENO3weight *W3 = WENOweights_global_p_WENO3_at_(im);
  return rec1d_p_WENO3(n,u, im, u_scale, W3);
}
/* more general WENO that uses weights in global WENOweights var */
double rec1d_p_WENO3_at_last_minus_l(int n, const double *u, int l,
                                     double u_scale)
{
  int nm = n-1;
  int im = nm-1 - l;
  tWENO3weight *W3 = WENOweights_global_p_WENO3_at_last_minus_(l);
  return rec1d_p_WENO3(n,u, im, u_scale, W3);
}

/* use rec1d_p_WENO3 and set weights for uniform grid */
double rec1d_p_WENO3_uniform(int n, const double *u, int im, double u_scale)
{
  tWENO3weight W3[1];
  W3->lw[0][0] = -0.5;
  W3->lw[0][1] = 1.5;
  W3->lw[1][0] = 0.5;
  W3->lw[1][1] = 0.5;
  W3->optw[0] = WENO3_3id_gamma1;
  W3->optw[1] = WENO3_3id_gamma2;
  return rec1d_p_WENO3(n,u, im, u_scale, W3);
}


/* Interpolate a field u to midpoint i+1/2 with index im.
   Here we interpolate in the negative direction (m) from the right of the
   midpoint to the midpoint to obtain umid_m.
   We use u at the points i, i+1, i+2 */
/* more general WENO that uses weights in global WENOweights var */
double rec1d_m_WENO3(int n, const double *u, int im, double u_scale,
                     tWENO3weight *W3)
{
  /* u at 3 grid points around point i=im */
  double u_i   = u[im];
  double u_ip1 = u[im+1];
  double u_ip2 = u[im+2];

  /* diffs */
  double d1 = u_ip1 - u_i; // maybe we should divide here by grid spacing???
  double d2 = u_ip2 - u_ip1;

  /* smoothness indicators */
  double beta1 = d1*d1;
  double beta2 = d2*d2;
  double us2 = u_scale*u_scale;
  double beta1_p_eps = beta1 + WENO3_epsilon*us2;
  double beta2_p_eps = beta2 + WENO3_epsilon*us2;

  /* non-normalized weights */
  double omegab1 = W3->optw[0]/(beta1_p_eps*beta1_p_eps);
  double omegab2 = W3->optw[1]/(beta2_p_eps*beta2_p_eps);
  double omegab_sum = omegab1 + omegab2;

  /* normalized weights */
  double omega1 = omegab1/(omegab_sum);
  double omega2 = omegab2/(omegab_sum);

  /* 2 linear reconstructions */
  double u1 = W3->lw[0][0]*u_i   + W3->lw[0][1]*u_ip1;
  double u2 = W3->lw[1][0]*u_ip1 + W3->lw[1][1]*u_ip2;

  /* final reconstruction */
  return omega1*u1 + omega2*u2;
}

/* more general WENO that uses weights in global WENOweights var */
double rec1d_m_WENO3_at_im(int n, const double *u, int im, double u_scale)
{
  tWENO3weight *W3 = WENOweights_global_m_WENO3_at_(im);
  return rec1d_m_WENO3(n,u, im, u_scale, W3);
}
/* more general WENO that uses weights in global WENOweights var */
double rec1d_m_WENO3_at_last_minus_l(int n, const double *u, int l,
                                     double u_scale)
{
  int nm = n-1;
  int im = nm-1 - l;
  tWENO3weight *W3 = WENOweights_global_m_WENO3_at_last_minus_(l);
  return rec1d_m_WENO3(n,u, im, u_scale, W3);
}

/* use rec1d_m_WENO3 and set weights for uniform grid */
double rec1d_m_WENO3_uniform(int n, const double *u, int im, double u_scale)
{
  tWENO3weight W3[1];
  W3->lw[0][0] = 0.5;
  W3->lw[0][1] = 0.5;
  W3->lw[1][0] = 1.5;
  W3->lw[1][1] = -0.5;
  W3->optw[0] = WENO3_3id_gamma2;
  W3->optw[1] = WENO3_3id_gamma1;
  return rec1d_m_WENO3(n,u, im, u_scale, W3);
}



/* Use WENO3 inside and WENO1 (i.e. just copying) near the boundary.
   Even if the stencil would fit we stay one more away from the boundary
   because the grid points at i=0 and i=n-1 are considerd to be moved in
   by h/4. */
double rec1d_p_WENO3_if2away(int n, const double *u, int im, double u_scale)
{
  if(im>1 && im<n-2) return rec1d_p_WENO3_uniform(n, u, im, u_scale);
  else               return rec1d_p_1(n, u, im, u_scale);
}
double rec1d_m_WENO3_if2away(int n, const double *u, int im, double u_scale)
{
  if(im>0 && im<n-3) return rec1d_m_WENO3_uniform(n, u, im, u_scale);
  else               return rec1d_m_1(n, u, im, u_scale);
}

/* Use WENO3 inside and WENO1 (i.e. just copying) at the boundary.
   The grid points at i=0 and i=n-1 are considerd to be moved in by h/4. */
double rec1d_p_WENO3_if1away(int n, const double *u, int im, double u_scale)
{
  if(im>1 && im<n-2)
  {
    return rec1d_p_WENO3_uniform(n, u, im, u_scale);
  }
  else
  {
    if(im==1)         return rec1d_p_WENO3_at_im(n, u, im, u_scale);
    else if(im==n-2)  return rec1d_p_WENO3_at_last_minus_l(n, u, 0, u_scale);
    else              return rec1d_p_1(n, u, im, u_scale);
  }
}
double rec1d_m_WENO3_if1away(int n, const double *u, int im, double u_scale)
{
  if(im>0 && im<n-3)
  {
    return rec1d_m_WENO3_uniform(n, u, im, u_scale);
  }
  else
  {
    if(im==0)        return rec1d_m_WENO3_at_im(n, u, im, u_scale);
    else if(im==n-3) return rec1d_m_WENO3_at_last_minus_l(n, u, 1, u_scale);
    else             return rec1d_m_1(n, u, im, u_scale);
  }
}


/* Use WENO3 inside and some 2nd order thing at the boundary.
   The grid points at i=0 and i=n-1 are considerd to be moved in by h/4.
   The n-2 midpoints are at im=0,...,n-2
   The 2 face points are at im=-1 & im = n-1 */
double rec1d_p_WENO3_2(int n, const double *u, int im, double u_scale)
{
  int lm;

  /* inside */
  if(im>1 && im<n-2)
    return rec1d_p_WENO3_uniform(n, u, im, u_scale);

  /* on right end: lm=0 is last midpoint lm=-1 is facepoint */
  lm = (n-2) - im;
  switch(lm)
  {
  case 0:
    if(n>2) return rec1d_p_WENO3_at_last_minus_l(n, u, 0, u_scale);
    else    return rec1d_p_1(n, u, im, u_scale);
  case -1:
    {
      //errorexit("rt???");

      double lw0 = -0.33333333333333333333333333333;
      double lw1 = 1. - lw0;
      return lw0*u[n-2] + lw1*u[n-1];
    }
  }

  /* on left end */
  switch(im)
  {
  case 1:
    return rec1d_p_WENO3_at_im(n, u, im, u_scale);
  case 0:
    {
      //return rec1d_p_1(n, u, im, u_scale);; // <--old

      double lw1 = 0.33333333333333333333333333333;
      double lw0 = 1. - lw1;
      return lw0*u[0] + lw1*u[1];
    }
  case -1:
    {
      errorexit("p_WENO3_2 not implemented at left face");
      double lw1 = -0.33333333333333333333333333333;
      double lw0 = 1. - lw1;
      return lw0*u[0] + lw1*u[1];
    }
  }
  /* should never get here: */
  errorexiti("im=%d out of range", im);
}
double rec1d_m_WENO3_2(int n, const double *u, int im, double u_scale)
{
  int lm;

  /* inside */
  if(im>0 && im<n-3)
    return rec1d_m_WENO3_uniform(n, u, im, u_scale);

  /* on left end */
  switch(im)
  {
  case 0:
    if(n>2) return rec1d_m_WENO3_at_im(n, u, im, u_scale);
    else    return rec1d_m_1(n, u, im, u_scale);
  case -1:
    {
      //errorexit("lft???");

      double lw1 = -0.33333333333333333333333333333;
      double lw0 = 1. - lw1;
      return lw0*u[0] + lw1*u[1];
    }
  }

  /* on right end: lm=0 is last midpoint lm=-1 is facepoint */
  lm = (n-2) - im;
  switch(lm)
  {
  case 1:
    return rec1d_m_WENO3_at_last_minus_l(n, u, 1, u_scale);
  case 0:
    {
      //return rec1d_m_1(n, u, im, u_scale);; // <--old

      double lw0 = 0.33333333333333333333333333333;
      double lw1 = 1. - lw0;
      return lw0*u[n-2] + lw1*u[n-1];
    }
  case -1:
    {
      errorexit("m_WENO3_2 not implemented at right face");
      double lw0 = -0.33333333333333333333333333333;
      double lw1 = 1. - lw0;
      return lw0*u[n-2] + lw1*u[n-1];
    }
  }
  /* should never get here: */
  errorexiti("im=%d out of range", im);
}


/* weight for the funny WENO3 at the outer side of the first/last midpoint */
#define WENO3_DGghost_gamma1 1.
#define WENO3_DGghost_gamma2 1.
/* Question: is this correct??? */

/* Use WENO3 inside, some 2nd order thing at the boundary, and a WENO3
   with the WENO3_DGghost_gamma1/2 weights on the outer side of the
   first/last midpoint.
   The grid points at i=0 and i=n-1 are considered to be moved in by h/4.
   The n-2 midpoints are at im=0,...,n-2
   The 2 face points are at im=-1 & im = n-1 */
double rec1d_p_WENO3_2g(int n, const double *u, int im, double u_scale)
{
  int lm;

  /* inside */
  if(im>1 && im<n-2)
    return rec1d_p_WENO3_uniform(n, u, im, u_scale);

  /* on right end: lm=0 is last midpoint lm=-1 is facepoint */
  lm = (n-2) - im;
  switch(lm)
  {
  case 0:
    if(n>2) return rec1d_p_WENO3_at_last_minus_l(n, u, 0, u_scale);
    else    return rec1d_p_1(n, u, im, u_scale);
  case -1:
    {
      //errorexit("rt???");

      double lw0 = -0.33333333333333333333333333333;
      double lw1 = 1. - lw0;
      return lw0*u[n-2] + lw1*u[n-1];
    }
  }

  /* on left end */
  switch(im)
  {
  case 1:
    return rec1d_p_WENO3_at_im(n, u, im, u_scale);
  case 0:
    {
      tWENO3weight W3[1];
      W3->lw[0][0] = -1.;
      W3->lw[0][1] = 2.;
      W3->lw[1][1] = 0.33333333333333333333333333333;
      W3->lw[1][0] = 1. - W3->lw[1][1];
      W3->optw[0] = WENO3_DGghost_gamma1;
      W3->optw[1] = WENO3_DGghost_gamma2;
      return rec1d_p_WENO3(n,u, im, u_scale, W3);

      /*old : */
      double lw1 = 0.33333333333333333333333333333;
      double lw0 = 1. - lw1;
      return lw0*u[0] + lw1*u[1];
    }
  case -1:
    {
      errorexit("p_WENO3_2g not implemented at left face");
      double lw1 = -0.33333333333333333333333333333;
      double lw0 = 1. - lw1;
      return lw0*u[0] + lw1*u[1];
    }
  }
  /* should never get here: */
  errorexiti("im=%d out of range", im);
}
double rec1d_m_WENO3_2g(int n, const double *u, int im, double u_scale)
{
  int lm;

  /* inside */
  if(im>0 && im<n-3)
    return rec1d_m_WENO3_uniform(n, u, im, u_scale);

  /* on left end */
  switch(im)
  {
  case 0:
    if(n>2) return rec1d_m_WENO3_at_im(n, u, im, u_scale);
    else    return rec1d_m_1(n, u, im, u_scale);
  case -1:
    {
      //errorexit("lft???");

      double lw1 = -0.33333333333333333333333333333;
      double lw0 = 1. - lw1;
      return lw0*u[0] + lw1*u[1];
    }
  }

  /* on right end: lm=0 is last midpoint lm=-1 is facepoint */
  lm = (n-2) - im;
  switch(lm)
  {
  case 1:
    return rec1d_m_WENO3_at_last_minus_l(n, u, 1, u_scale);
  case 0:
    {
      tWENO3weight W3[1];
      W3->lw[0][0] = 0.33333333333333333333333333333;
      W3->lw[0][1] = 1. - W3->lw[0][0];
      W3->lw[1][0] = 2.;
      W3->lw[1][1] = -1.;
      W3->optw[0] = WENO3_DGghost_gamma2;
      W3->optw[1] = WENO3_DGghost_gamma1;
      return rec1d_m_WENO3(n,u, im, u_scale, W3);

      /* old: */
      double lw0 = 0.33333333333333333333333333333;
      double lw1 = 1. - lw0;
      return lw0*u[n-2] + lw1*u[n-1];
    }
  case -1:
    {
      errorexit("m_WENO3_2g not implemented at right face");
      double lw0 = -0.33333333333333333333333333333;
      double lw1 = 1. - lw0;
      return lw0*u[n-2] + lw1*u[n-1];
    }
  }
  /* should never get here: */
  errorexiti("im=%d out of range", im);
}



/* forward=1: convert from u at face to u at 0.25h in
   forward=0: convert from u at 0.25h in to u at face */
void rec1d_uface_to_uin_1(tNode *node, tVarList *vlu, int forward)
{
  if(node->dat->info->use_fv)
  {
    int *n = node->n;
    int f;
    double c0 = 0.75;
    double c1 = 1. - c0;
    double w0,w1;

    if(forward) /* weights for linear interpolation*/
    {
      w0 = c0;
      w1 = c1;
    }
    else /* weights for linear extrapolation*/
    {
      w0 = 1./c0;
      w1 = -c1*w0;
    }

    /* loop over 6 faces */
    for(f=0; f<6; f++)
    {
      int face = forward ? f : 5-f;
      int dir = face/2;
      int top = (face%2);
      int sign = 2*top - 1;
      int pl = top*(n[dir] - 1);
      int i,j,k;

      /* no interpolation if only 1 or 2 points */
      if(n[dir]<=2) continue;

      /* loop over plane p */
      forplaneN(dir, i,j,k, n, pl)
      {
        int vi; /* field index */
        int i1 = i1_norm(i,j,k, dir); /* 1st and 2nd index in plane */
        int i2 = i2_norm(i,j,k, dir);
        int i0 = i0_norm(i,j,k, dir); /* index orthogonal to plane */
        int i0in = i0 - sign;         /* one away from from face */
        int ic,jc,kc, ccc, cccin;

        /* get index for i0 and i0in */
        ccc = Ind_n(i,j,k, n);
        ijk_inplaneN(dir, ic,jc,kc, i1,i2, i0in);
        cccin = Ind_n(ic,jc,kc, n);

        /* interpolate inward or extrapolate outward */
        forvl(vlu, vi)
        {
          int iu = Vind(vlu, vi);
          double *u = Vard_(node, iu);

          u[ccc] = w0*u[ccc] + w1*u[cccin];
        }
      } /* end plane loop */
    } /* end face-loop*/
  }
}

/* forward=1: convert from u at face to u at 0.25h in
   forward=0: convert from u at 0.25h in to u at face */
void rec1d_uface_to_uin_1_mesh(tMesh *mesh, tVarList *vlu, int forward)
{
  formylnodes(mesh)
  {
    tNode *node = MyLnode;
    rec1d_uface_to_uin_1(node, vlu, forward);
  }
}


/*********************************************************************/
/* WENOm3 is much like WENO3 but with different ideal weights */
/*********************************************************************/
#define WENOm3_3id_gamma1 0.25
#define WENOm3_3id_gamma2 0.75

/* use rec1d_p_WENO3 with weights for uniform grid */
double rec1d_p_WENOm3_uniform(int n, const double *u, int im, double u_scale)
{
  tWENO3weight W3[1];
  W3->lw[0][0] = -0.5;
  W3->lw[0][1] = 1.5;
  W3->lw[1][0] = 0.5;
  W3->lw[1][1] = 0.5;
  W3->optw[0] = WENOm3_3id_gamma1;
  W3->optw[1] = WENOm3_3id_gamma2;
  return rec1d_p_WENO3(n,u, im, u_scale, W3);
}

/* use rec1d_m_WENO3 with weights for uniform grid */
double rec1d_m_WENOm3_uniform(int n, const double *u, int im, double u_scale)
{
  tWENO3weight W3[1];
  W3->lw[0][0] = 0.5;
  W3->lw[0][1] = 0.5;
  W3->lw[1][0] = 1.5;
  W3->lw[1][1] = -0.5;
  W3->optw[0] = WENOm3_3id_gamma2;
  W3->optw[1] = WENOm3_3id_gamma1;
  return rec1d_m_WENO3(n,u, im, u_scale, W3);
}


/* Use WENOm3 inside and a reversed stencil at the boundary.
   The n-2 midpoints are at im=0,...,n-2
   The 2 face points are at im=-1 & im = n-1 */
double rec1d_p_WENOm3_2(int n, const double *u, int im, double u_scale)
{
  int lm;

  /* inside */
  if(im>0 && im<=n-2)
    return rec1d_p_WENOm3_uniform(n, u, im, u_scale);

  if(n<3) return rec1d_p_1(n, u, im, u_scale);

  /* on left end */
  switch(im)
  {
  case 0:
    /* reverse stencil and use m_WENO */
    return rec1d_m_WENOm3_uniform(n, u, im, u_scale);
  case -1:
    return u[0];
  }

  /* on right end: lm=0 is last midpoint lm=-1 is facepoint */
  lm = (n-2) - im;
  switch(lm)
  {
  case -1:
    return u[n-1];
  }

  /* should never get here: */
  errorexiti("im=%d out of range", im);
}
double rec1d_m_WENOm3_2(int n, const double *u, int im, double u_scale)
{
  int lm;

  /* inside */
  if(im>=0 && im<=n-3)
    return rec1d_m_WENO3_uniform(n, u, im, u_scale);

  if(n<3) return rec1d_m_1(n, u, im, u_scale);

  /* on right end: lm=0 is last midpoint lm=-1 is facepoint */
  lm = (n-2) - im;
  switch(lm)
  {
  case 0:
    /* reverse stencil and use p_WENO */
    return rec1d_p_WENOm3_uniform(n, u, im, u_scale);
  case -1:
    return u[n-1];
  }

  /* on left end */
  switch(im)
  {
  case -1:
    return u[0];
  }

  /* should never get here: */
  errorexiti("im=%d out of range", im);
}
