/* Chebyshev.c */
/* Wolfgang Tichy, 01/2025
   does explicit slow Chebyshev trafos: */

#include "nmesh.h"
#include "basis.h"



/* ************************************************************************ */
/* various functions for expansions in Chebyshev polynomials                */
/* ************************************************************************ */

/* Note here we have X = [a,b] (with a=-1, b=1) and x=[1,-1]
   X = (a-b)*x/2 + (a+b)/2   so that: x=-1 => X=b,  x=1 => X=a
   x = (2X-b-a)/(a-b)
   This does not imply that the result of all Cheb. sums
   (e.g. u[j] = sum_k c[k] T_k[x_j]
         c[j] ~ sum_k u[x_k] T_j[x_k] )
   gets a minus sign compared to what is done in numerical recipes,
   even though T_k(-x) = (-1)^k T_k(x).
*/

/* get Cheb coeffs c[0...np-1] of function f(X) (X in [-1,1]) */
void Chebyshev_getcoeff(int np, double c[], double (*func)(double))
{
  int k,j;
  double fac, *f;

  f = (double*) calloc(np,sizeof(double));

  //amb = -1.; // = 0.5*(a-b);
  //bpa = 0.;  // = 0.5*(b+a);

  /* use function at zeros: */
  for(k=0; k<np; k++)
  {
    double y=cos(PI*(k+0.5)/np);
    f[k]=(*func)(-y);
  }
  fac=2.0/np;
  for(j=0; j<np; j++)
  {
    double sum=0.0;
    for(k=0; k<np; k++)
      sum += f[k]*cos(PI*j*(k+0.5)/np);
    c[j]=fac*sum;
  }
  free(f);
}


/* find value of function at X (in [a,b]) from Cheb coeffs c[0...np-1] */
double Chebyshev_eval(int np, double c[], double X)
{
  double d=0.0, dd=0.0, sv, y, y2;
  int j;

  if ((X+1.)*(X-1.) > 0.0) printf("X not in range in routine Chebyshev_eval\n");
  y = -X; // =(2.0*X-a-b)/(a-b);
  y2= 2.0*y;
  for(j=np-1; j>=1; j--)
  {
    sv=d;
    d=y2*d-dd+c[j];
    dd=sv;
  }
  return y*d-dd+0.5*c[0];
}


/* compute Cheb coeffs of deriv cder[0...n] from Cheb coeffs c[0...n] */
void Chebyshev_deriv(int np, double c[], double cder[])
{
  int n = np-1;
  int j;
  double con;

  cder[n]=0.0;
  if(n>=1) cder[n-1]=2*(n)*c[n];
  for (j=n-2;j>=0;j--)
    cder[j]=cder[j+2]+2*(j+1)*c[j+1];

  /* convert to interval [a,b] */
  con = -1.; // = 2.0/(a-b);
  for (j=0;j<n;j++)
    cder[j] *= con;
}


/* compute Cheb coeffs of integral cint[0...n] from Cheb coeffs c[0...n] */
void Chebyshev_integ(int np, double c[], double cint[])
{
  int n = np-1;
  int j;
  double sum=0.0, fac=1.0, con;

  con = -0.5; // = 0.25*(a-b);
  for (j=1;j<=n-1;j++)
  {
    cint[j]=con*(c[j-1]-c[j+1])/j;
    sum += fac*cint[j];
    fac = -fac;
  }
  cint[n]=con*c[n-1]/(n);
  sum += fac*cint[n];
  cint[0]=2.0*sum;  /* <--arbitrary const picked as in numrec */
}


/* compute Cheb coeffs c[0...n] from function u at the zeros of T_{np}(x).
   Note np=n+1 */
void Chebyshev_coeffs_fromZeros(int np, double u[], double c[])
{
  int k, j;
  double fac, sum, PIonp;

  PIonp=PI/np;
  fac=2.0/np;

  for(j=0;j<np;j++)
  {
    sum=0.0;
    for(k=0;k<np;k++)
      sum += cos(j*PIonp*(k+0.5))*u[k];
    c[j]=fac*sum;
  }
  /* we should use a FFT for eveything above this line */
}


/* compute Cheb coeffs c[0...N] from function u at the extrema of T_N(x). */
void Chebyshev_coeffs_fromExtrema(int np, double u[], double c[])
{
  int N = np-1;
  int k, j;
  double fac, sum, PIoN;

  PIoN=PI/N;
  fac=2.0/N;

  for(j=0;j<=N;j++)
  {
    sum = 0.5 * (u[0] + cos(j*PI)*u[N]);
    for(k=1;k<N;k++)
      sum += cos(j*PIoN*k)*u[k];
    c[j]=fac*sum;
  }
  c[N]*=0.5;
  /* we should use a FFT for eveything above this line */
}


/* find function u on the zeros of T_{np}(x) */
void Chebyshev_eval_onZeros(int np, double u[], double c[])
{
  int k, j;
  double sum;
  double PIonp=PI/np;
  //double fac=2.0/np;

  for(j=0;j<np;j++)
  {
    sum=0.5*c[0];
    for(k=1;k<np;k++)
      sum += c[k]*cos(k*PIonp*(j+0.5));
    u[j]=sum;
  }
  /* we should use a FFT for eveything above this line */
}


/* find function u on the extrema of T_N(X) */
void Chebyshev_eval_onExtrema(int np, double u[], double c[])
{
  int N = np-1;
  int k, j;
  double sum;
  double PIoN=PI/N;
  //double fac=2.0/N;

  for(j=0;j<=N;j++)
  {
    sum = 0.5*c[0];
    for(k=1;k<=N;k++)
      sum += c[k]*cos(k*PIoN*j);
    u[j]=sum;
  }
  /* we should use a FFT for eveything above this line */
}


/* filter: zero all c[j] with k<=j<np */
void Chebyshev_filter(double c[], int k, int np)
{
  int j;
  for(j=k; j<np; j++) c[j] = 0.;
}


/* find value of Cheb. basis function T_n at X (in [a,b]) */
double Chebyshev_basisfunc_FromSum(int n, double X, int np)
{
  double d=0.0, dd=0.0, sv, y, y2;
  int j;

  if(n==0) return 0.5;  /* in numrec T_0 / 2 is used as basisfunc # 0 */

  if ((X+1.)*(X-1.) > 0.0) printf("X not in range in routine Chebyshev_basis\n");
  y = -X; // =(2.0*X-a-b)/(a-b);
  y2=2.0*y;

  sv=d;
  d=y2*d-dd+ 1.0;
  dd=sv;

  for(j=n-1;j>=1;j--)
  {
    sv=d;
    d=y2*d-dd;
    dd=sv;
  }
  return y*d-dd;
}


/* find value of Cheb. basis function T_n at X (in [a,b]) */
double Chebyshev_basisfunc(int n, double X, int np)
{
  double y;
  if(n==0) return 0.5;  /* in numrec T_0 / 2 is used as basisfunc # 0 */

  /* y=(2.0*X-a-b)/(a-b); */
  y = -X;

  /* catch small rounding errors, where |y| is a bit larger than 1 */
  if(dequal(y, 1.0)) return 1.0;
  if(dequal(y,-1.0))
  {
    if(n%2) return -1.0;
    else    return +1.0;
  }
  //if(y>1. || y<-1.) printf("Chebyshev_basisfunc: y=%.19g out of range\n", y);

  return cos(n*acos(y));
}
