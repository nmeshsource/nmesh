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
double rec1d_p_WENO3(int n, const double *u, int im, double u_scale)
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
double rec1d_m_WENO3(int n, const double *u, int im, double u_scale)
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


/* Use WENO3 inside and WENO1 (i.e. just copying) near the boundary.
   Even if the stencil would fit we stay one more away from the boundary
   because the grid points at i=0 and i=n-1 are considerd to be moved in
   by h/4. */
double rec1d_p_WENO3_1(int n, const double *u, int im, double u_scale)
{
  if(im>1 && im<n-2) return rec1d_p_WENO3(n, u, im, u_scale);
  else               return rec1d_p_1(n, u, im, u_scale);
}
double rec1d_m_WENO3_1(int n, const double *u, int im, double u_scale)
{
  if(im>0 && im<n-3) return rec1d_m_WENO3(n, u, im, u_scale);
  else               return rec1d_m_1(n, u, im, u_scale);
}
