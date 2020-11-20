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

/* use rec1d_p_WENO3 and set weights uniform grid */
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

/* use rec1d_m_WENO3 and set weights uniform grid */
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
