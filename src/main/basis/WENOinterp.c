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
//OLD: (int N, double x, double x0, double h, double *c, double *u)
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


/***********************************************************************/
/* func to select correct WENO */
/***********************************************************************/

/* nth order WENO */
double interpolate_WENO_n_ds(double x, int n, const double *x_p,
                             const double *w_interp,
                             const double *f, int ds, double fscal)
{
  switch(n)
  {
  case 1:
  case 2:
    return Lagrange_interp_barycentric2_ds(x, n,x_p, w_interp, f, ds, fscal);
  case 4:
    return interpolate_WENO_4_ds(x, n,x_p, w_interp, f, ds, fscal);
  case 6:
    return interpolate_WENO_6_ds(x, n,x_p, w_interp, f, ds, fscal);
  default:
    errorexiti("order n=%d is not implemented", n);
  }
}
