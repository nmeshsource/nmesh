/* WENOinterp.c */
/* Wolfgang Tichy, 2/2024 */

#include "nmesh.h"
#include "basis.h"

#define sqr(x) ((x)*(x))


/***********************************************************************/
/* funcs for 1d WENO interpolation */
/***********************************************************************/

/* WENO6 interpolation as in
   Colin B. Macdonald, Steven J. Ruth, J Sci Comput (2008) 35: 219-240 */
double interpolate_WENO_6_ds(double x, int n, const double *x_p,
                             const double *w_interp,
                             const double *f, int ds, double f_scale)
{
  double u[] = { f[0], f[ds], f[2*ds], f[3*ds], f[4*ds], f[5*ds] };
  double fs2 = f_scale*f_scale;
  double epsl = 1e-6 * fs2;

  double p[4], C[4], IS[4], alpha[4], w[6];
  double h = x_p[1] - x_p[0]; /* assume equal spacing */
  double ooh  = 1./h;
  double ooh2 = ooh*ooh;
  double ooh3 = ooh2*ooh;

  double dx0 = x - x_p[0];   // x - x_i-2
  double dx1 = x - x_p[1];   // x - x_i-1
  double dx2 = x - x_p[2];   // x - x_i
  double dx3 = x - x_p[3];   // x - x_i+1
  double dx4 = x - x_p[4];   // x - x_i+2
  double dx5 = x - x_p[5];   // x - x_i+3
  double sumalpha, sum;
  int i;

  //i-2 == 0 ...
  p[1] = u[0] + (u[1]-u[0])*ooh*dx0 + (u[2]-2.*u[1]+u[0])*ooh2/2.*dx0*dx1 + (u[3]-3.*u[2]+3.*u[1]-u[0])*ooh3/6.*dx0*dx1*dx2;
  p[2] = u[1] + (u[2]-u[1])*ooh*dx1 + (u[3]-2.*u[2]+u[1])*ooh2/2.*dx1*dx2 + (u[4]-3.*u[3]+3.*u[2]-u[1])*ooh3/6.*dx1*dx2*dx3;
  p[3] = u[2] + (u[3]-u[2])*ooh*dx2 + (u[4]-2.*u[3]+u[2])*ooh2/2.*dx2*dx3 + (u[5]-3.*u[4]+3.*u[3]-u[2])*ooh3/6.*dx2*dx3*dx4;

  C[1] = dx4*dx5*ooh2/20.;
  C[2] =-dx5*dx0*ooh2/10.;
  C[3] = dx0*dx1*ooh2/20.;

  IS[1] = (814.*sqr(u[3]) + 4326.*sqr(u[2]) + 2976.*sqr(u[1]) + 244.*sqr(u[0]) - 3579.*u[2]*u[3] - 6927.*u[2]*u[1]
          + 1854.*u[2]*u[0] + 2634.*u[3]*u[1] - 683.*u[3]*u[0] - 1659.*u[1]*u[0])/180.;
  IS[2] = (1986.*sqr(u[3]) + 1986.*sqr(u[2]) + 244.*sqr(u[1]) + 244.*sqr(u[4]) + 1074.*u[2]*u[4] - 3777.*u[2]*u[3]
          - 1269.*u[2]*u[1] + 1074.*u[3]*u[1] - 1269.*u[4]*u[3] - 293.*u[4]*u[1])/180.;
  IS[3] = (814.*sqr(u[2]) + 4326.*sqr(u[3]) + 2976.*sqr(u[4]) + 244.*sqr(u[5]) - 683.*u[2]*u[5] + 2634.*u[2]*u[4]
            - 3579.*u[2]*u[3] - 6927.*u[3]*u[4] + 1854.*u[3]*u[5] - 1659.*u[4]*u[5])/180.;

  sumalpha = 0;
  for(i=1; i<=3; i++)
  {
    alpha[i] = C[i]/sqr(epsl + IS[i]);
    sumalpha += alpha[i];
  }

  for(i=1; i<=3; i++)
    w[i] = alpha[i]/sumalpha;

//for(i=0; i<n; i++) printf("x_p[%d]=%g ", i, x_p[i]);
//printf("\n");
/*
for(i=1; i<=3; i++) printf("w[%d]=%g ", i, w[i]);
printf("\n");
for(i=1; i<=3; i++) printf("p[%d]=%g ", i, p[i]);
printf("\n");
*/
  sum = 0;
  for(i=1; i<=3; i++)
    sum += w[i]*p[i];

  return sum;
}


/* WENO4 interpolation as in
   Colin B. Macdonald, Steven J. Ruth, J Sci Comput (2008) 35: 219-240 */
double interpolate_WENO_4_ds(double x, int n, const double *x_p,
                             const double *w_interp,
                             const double *f, int ds, double f_scale)
{
  double u[] = { f[0], f[ds], f[2*ds], f[3*ds] };
  double fs2 = f_scale*f_scale;
  double epsl = 1e-6 * fs2;

  double p[3],C[3],IS[3],alpha[3],w[4];
  double h = x_p[1] - x_p[0]; /* assume equal spacing */
  double ooh  = 1./h;
  double ooh2 = ooh*ooh;

  double dx0 = x - x_p[0];   // x - x_i-1
  double dx1 = x - x_p[1];   // x - x_i
  //double dx2 = x - x_p[2];   // x - x_i+1
  double dx3 = x - x_p[3];   // x - x_i+2
  double sumalpha, sum;
  int i;

  //i-1 == 0 ...
  p[1] = u[1] + (u[2]-u[0])*ooh/2.*dx1             + (u[2]-2.*u[1]+u[0])*ooh2/2.*dx1*dx1;
  p[2] = u[1] + (-u[3]+4.*u[2]-3.*u[1])*ooh/2.*dx1 + (u[3]-2.*u[2]+u[1])*ooh2/2.*dx1*dx1;

  C[1] =-dx3*ooh2/3.;
  C[2] = dx0*ooh2/3.;

  IS[1] = (25.*sqr(u[2]) + 64.*sqr(u[1]) + 13.*sqr(u[0]) +
          26.*u[2]*u[0] - 52.*u[1]*u[0] - 76.*u[2]*u[1])/12.;
  IS[2] = (25.*sqr(u[1]) + 64.*sqr(u[2]) + 13.*sqr(u[3]) +
          26.*u[3]*u[1] - 52.*u[3]*u[2] - 76.*u[2]*u[1])/12.;

  sumalpha = 0;
  for (i=1; i<=2; i++)
  {
    alpha[i] = C[i]/sqr(epsl + IS[i]);
    sumalpha += alpha[i];
  }

  for (i=1; i<=2; i++)
    w[i] = alpha[i]/sumalpha;

  sum = 0;
  for (i=1; i<=2; i++)
    sum += w[i]*p[i];

  return sum;
}


/* How we can get WENO3 weights:
# 3 points and their function values
pt = [ 1,2,3 ]
fv = [ 1,4,9 ]

# 2nd order Lagrange poly
P2(x) =
fv[0]*(x-pt[1])*(x-pt[2])/((pt[0]-pt[1])*(pt[0]-pt[2])) +
fv[1]*(x-pt[0])*(x-pt[2])/((pt[1]-pt[0])*(pt[1]-pt[2])) +
fv[2]*(x-pt[0])*(x-pt[1])/((pt[2]-pt[0])*(pt[2]-pt[1]))

# 2 1st order Lagrange polys
# L0(x) = fv[0]*(x-pt[1])/((pt[0]-pt[1])) + fv[1]*(x-pt[0])/((pt[1]-pt[0]))
# L1(x) = fv[1]*(x-pt[2])/((pt[1]-pt[2])) + fv[2]*(x-pt[1])/((pt[2]-pt[1]))

# coeffs for linear interpolation
a00(x) = (x-pt[1])/(pt[0]-pt[1])
a01(x) = (x-pt[0])/(pt[1]-pt[0])
a10(x) = (x-pt[2])/(pt[1]-pt[2])
a11(x) = (x-pt[1])/(pt[2]-pt[1])

# linear interpolation
L0(x) = fv[0]*a00(x) + fv[1]*a01(x)
L1(x) = fv[1]*a10(x) + fv[2]*a11(x)

# weights for L0 and L1 to reconstruct P2
g0(x) = (x-pt[2])/(pt[0]-pt[2])
g1(x) = (x-pt[0])/(pt[2]-pt[0])

# This S(x) is equal to P2(x):
S(x) = g0(x)*L0(x) + g1(x)*L1(x)
*/
/* centered WENO3 interpolation, i.e. not the one used to reconstruct fluxes
   in e.g. an upwind scheme */
double interpolate_cenWENO_3_ds(double x, int n, const double *x_p,
                                const double *w_interp,
                                const double *f, int ds, double f_scale)
{
  double u[] = { f[0], f[ds], f[2*ds] };
  double fs2 = f_scale*f_scale;
  double epsl = 1e-6 * fs2;

  /* grid spacings */
  double h0 = x_p[1]-x_p[0];
  double h1 = x_p[2]-x_p[1];
  double h0sqr = h0*h0;
  double h1sqr = h1*h1;

  /* coeffs for linear polys */
  double a00 = (x-x_p[1])/(-h0);
  double a01 = (x-x_p[0])/(+h0);
  double a10 = (x-x_p[2])/(-h1);
  double a11 = (x-x_p[1])/(+h1);
  /* 2 linear interpolations */
  double L0 = u[0]*a00 + u[1]*a01;
  double L1 = u[1]*a10 + u[2]*a11;
  /* optimal weights, to reconstruct 2nd order Lagrange poly */
  double g0 = (x-x_p[2])/(x_p[0]-x_p[2]);
  double g1 = (x-x_p[0])/(x_p[2]-x_p[0]);

  /* diffs */
  double d0 = u[1] - u[0];
  double d1 = u[2] - u[1];

  /* smoothness indicators */
  double beta0 = d0*d0;
  double beta1 = d1*d1;
  double beta0_p_eps = (beta0 + epsl)/h0sqr;
  double beta1_p_eps = (beta1 + epsl)/h1sqr;

  /* non-normalized weights */
  double omegab0 = g0/(beta0_p_eps*beta0_p_eps);
  double omegab1 = g1/(beta1_p_eps*beta1_p_eps);
  double omegab_sum = omegab0 + omegab1;

  /* normalized weights */
  double omega0 = omegab0/(omegab_sum);
  double omega1 = omegab1/(omegab_sum);

  /* final interpolation */
  return omega0*L0 + omega1*L1;
}

/***********************************************************************/
/* func to select correct WENO */
/***********************************************************************/

/* nth order WENO, Note: w_interp is ignored */
double interpolate_WENO_n_ds(double x, int n, const double *x_p,
                             const double *w_interp,
                             const double *f, int ds, double fscal)
{
  double w[2];
  switch(n)
  {
  case 1:
  case 2:
    /* set Lagrange interp. weights in case they are not in w_interp */
    Lagrange_winterp(n, x_p, w); // cheap for n=2
    return Lagrange_interp_barycentric2_ds(x, n,x_p, w, f, ds, fscal);
  case 4:
    return interpolate_WENO_4_ds(x, n,x_p, w_interp, f, ds, fscal);
  case 6:
    return interpolate_WENO_6_ds(x, n,x_p, w_interp, f, ds, fscal);
  default:
    errorexiti("order n=%d is not implemented", n);
  }
}
