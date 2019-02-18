/* Legendre.c */
/* Wolfgang Tichy, 2/2018
   some functions to deal with Legendre polynomials, especially at the
   Legendre Gauss-Lobatto (LGL) points. */

#include "nmesh.h"
#include "basis.h"



/* associated Legendre polynomial (not really needed here) */
/* use recursion recommended in numrec to get P_l^m for positive m */
double assocLegendreP(int l, int m, double x)
{
  double fact, Pkm, Pmm, Pmp1_m, somx2;
  int i, k;

  if(m < 0 || m > l || fabs(x) > 1.0) errorexit("Bad arguments");

  /* build P_m^m */
  Pmm = 1.0;
  if(m > 0)
  {
    somx2 = sqrt((1.0-x)*(1.0+x));
    fact=1.0;
    for(i=1; i<=m; i++)
    {
      Pmm *= -fact*somx2;
      fact += 2.0;
    }
  }
  /* are we done? */
  if(l == m) return Pmm;

  /* set P_{m+1}^m = x*(2*m+1) P_m^m */
  Pmp1_m = x*(2*m+1)*Pmm;
  /* is this all we need? */
  if(l == (m+1)) return Pmp1_m;

  /* use recursion:
     P_l^m = ( x*(2*l-1)*P_{l-1}^m - (l+m-1)*P_{l-2}^m )/(l-m)   */
  for(k=m+2; k<=l; k++) /* k runs through all l we have */
  {
    Pkm    = (x*(2*k-1)*Pmp1_m - (k+m-1)*Pmm)/(k-m);
    Pmm    = Pmp1_m;
    Pmp1_m = Pkm;
  } /* now k=l */
  return Pkm;
}

/* Legendre polynomial basis function, np is number of points (not need for
   Legendre, but maybe for Fourier) */
/* use recursion recommended in numrec to get P_l */
double basis_LegendreP(int l, double x, int np)
{
  double Pk, P0, P1;
  int k;

  /* build P_m^m */
  P0 = 1.0;
  if(l == 0) return P0;

  /* set P_{m+1} = x P_m */
  P1 = x; // P1 = x*P0;

  /* is this all we need? */
  if(l == 1) return P1;

  /* use recursion:
     P_l = ( x*(2*l-1)*P_{l-1} - (l-1)*P_{l-2} )/l   */
  for(k=2; k<=l; k++) /* k runs through all l we have */
  {
    Pk = (x*(2*k-1)*P1 - (k-1)*P0)/k;
    P0 = P1;
    P1 = Pk;
  } /* now k=l */
  return Pk;
}

/* basis_LegendreP without redundant np arg */
double LegendreP(int l, double x)
{
  if(fabs(x) > 1.0) errorexit("Bad arguments");
  return basis_LegendreP(l, x, 0);
}

/* normalized Legendre poly */
double basis_normLegendreP(int i, double x, int np)
{
  int N = np-1;
  double ci, sci;

  ci = 2.0/(2*i+1);
  if( (i==N) && (N!=0) ) ci = 2.0/i;
  sci = sqrt(ci);
  return basis_LegendreP(i, x, 0)/sci;
}

/* compute the values of the Legendre polynomial P=P_l(x), dP=dP_l(x)/dx,
   Q(x) = P_{l+1}(x) - P_{l-1}(x), and of dQ=Q'(x) */
void Legendre_P_dP_Q_dQ(int l, double x,
                       double *P, double *dP, double *Q, double *dQ)
{
  int k;
  double P0 = 1., P1 = x,  dP0 = 0., dP1 = 1.;
  double Pk, dPk;

  if(l==0)
  {
    *P  = 1.;
    *dP = 0.;
    return;
  }

  if(l==1)
  {
    *P  = x;
    *dP = 1.;
    return;
  }

  /* use recursion relation for P_l and its deriv */
  for(k = 2; k <= l+1; k++)
  {
    Pk  = (P1 * x * (2*k-1) - P0 * (k-1))/k;
    dPk = P1 * (2*k-1) + dP0;
    if(k>=l+1) break;
    P0  = P1;
    P1  = Pk;
    dP0 = dP1;
    dP1 = dPk;
  }

  *P  = P1;
  *dP = dP1;
  *Q  = Pk - P0;
  *dQ = dPk - dP0;
}


/* integral from x to 1 of Legendre polynomial P_k.
   This is coming from
   \frac{d}{dx}[(1-x^{2})\frac{dP_k}{dx}] + k(k+1) P_k(x) = 0 */
double int_LegendreP_x_1(int k, double x)
{
  double I, P,dP, q,e;

  if(k == 0)
  {
    I = 1. - x;
  }
  else
  {
    Legendre_P_dP_Q_dQ(k,x, &P,&dP, &q,&e);
    I = dP * (1. - x*x) / (k*(k+1));
  }
  return I;
}

/* integral from a to b of Legendre polynomial P_k */
double int_LegendreP_a_b(int k, double a, double b)
{
  return int_LegendreP_x_1(k, a) - int_LegendreP_x_1(k, b);
}


/* ************************************************************************ */
/* various functions needed for Legendre Gauss-Lobatto points or nodes      */
/* ************************************************************************ */

/* Legendre Gauss-Lobatto nodes x_i and integration weights w_i
   N = degree, so there are N+1 points
   This is based on
    (n+1) P_{n+1} = (2n+1) x P_{n} - n P_{n-1}
   which gives x P_{n} = \frac{n+1}{2n+1} P_{n+1} + \frac{n}{2n+1} P_{n-1}
   when x P_{n} is inserted into
    \frac{x^2-1}{n}\frac{d}{dx} P_{n} = x P_{n} - P_{n-1}
   we get:
    \frac{x^2-1}{n}\frac{d}{dx} P_{n} = \frac{n+1}{2n+1} (P_{n+1} - P_{n-1})
   So for x\in(-1,1) we get:
    \frac{d}{dx} P_{n}(x) = \frac{n(n+1)}{(x^2-1)(2n+1)} Q(x)
   So \frac{d}{dx} P_{n}(x) = 0 is equivalent to: Q(x) = 0 */
void LGL_x_winteg(int npoints, double *x, double *w)
{
  int N = npoints-1;
  int i,j;
  double P, dP, Q, dQ, y, dy;
  int Newton_itmax = 100;
  double Newton_rel_dytol = 1e-15;

  if(N == 0)
  {
    x[N] = 0.;
    if(w) w[N] = 2.;
    return;
  }

  x[0] = -1.;
  x[N] =  1.;
  if(w)
  {
    w[0] = 2.0/(N*(N+1.0));
    w[N] = w[0];
  }

  if(N == 1) return;


  /* get nodes in (-1,0] */
  for(j = 1; j <= N/2; j++)
  {
    if(j == N/2.0)  /* same as: if( (N%2==0 && j==N/2) ) */
    {
      y = 0.;
    }
    else
    {
      /* guess for node */
      y = - cos(PI*(j+0.25)/N - 3/(8*N*PI*(j+0.25)));

      /* find y such that Q(y) = 0,
         use Newton method to refine guess from line above */
      for(i = 0; i < Newton_itmax; i++)
      {
        Legendre_P_dP_Q_dQ(N, y, &P, &dP, &Q, &dQ);
        dy = - Q / dQ;
        y += dy;
        if(fabs(dy) < fabs(Newton_rel_dytol * y)) break;
      }
    }

    /* set x array for [-1,0] and [0,1] */
    x[j]   =  y;
    x[N-j] = -y;

    /* really just set P as in: P =  LegendreP(N, y); */
    Legendre_P_dP_Q_dQ(N, y, &P, &dP, &Q, &dQ);

    if(w)
    {
      w[j] = 2.0/( (N*(N+1)) * P*P);
      w[N-j] = w[j];
    }
  }
}

/* Gauss or Gauss-Lobatto (GL) quadrature:
   compute I = \int_{-1}^1 dx f(x) where f(x) is known at the nodes.
   For LGL nodes, the w[i] are the integration weights from LGL_x_winteg
   This is accurate for polynomials up to degree 2n-3 for LGL. */
double Gauss_integral(int n, const double *x, const double *w, const double *f)
{
  int i;
  double I = 0.;

  for(i = 0; i < n; i++)  I += w[i] * f[i];
  return I;
}



/* Set analysis and synthesis matrices for expansions in Legendre polynmials
   at Legendre Gauss-Lobatto nodes.

   The usual Legendre polynomials are not normalized
     (P_i,P_j)_N = c_i delta_{ij}
     where:  c_i = 2/(2i+1)
     except: c_N = 2/i       which is used for numerical integration

     (P_i,P_j)_numerical = \sum_{k=0}^N w_k P_i(x_k) P_j(x_k) = c_i delta_{ij}

   Without normalization we use this:
     u_i    = S_{ij} uhat_j,  S_{ij} = P_j(x_i)
     uhat_i = A_{ij} u_j,     A_{ij} = S^{-1}_{ij} = w_j S_{ji} / c_i

   Some prefer to normalize, so with normalization:
     u_i    = Shat_{ij} uhat_j,  Shat_{ij} = P_j(x_i)/sqrt(c_j)
     uhat_i = Ahat_{ij} u_j,     Ahat_{ij} = Shat^{-1}_{ij} = w_j Shat_{ji}

   In the function below store the transposes A^T and S^T of A and S
   if AT and ST are interpreted as stored in column-major form */
void LGL_AT_ST_matrices(int n, double *x, double *w, double *AT, double *ST)
{
  int normalized = 1;
  int N = n-1;
  int i, j;
  double ci, sci;

  /* row-loop */
  for(i = 0; i < n; i++)
  {
    ci = 2.0/(2*i+1);
    if( (i==N) && (N!=0) ) ci = 2.0/i;
    sci = sqrt(ci);

    /* col-loop */
    for(j = 0; j < n; j++)
    {
      if(normalized)
      {
	/* transpose of synthesis matrix Sij */
	/* set transposed elements first because we need them below */
	ST[i + j*n] = LegendreP(i, x[j]) / sci;
	/* transpose of analysis matrix Aij */
	AT[j + i*n] = w[j] * ST[i + j*n];
      }
      else
      {
	ST[i + j*n] = LegendreP(i, x[j]);
	AT[j + i*n] = w[j] * ST[i + j*n] / ci;
      }
    }
  }
}

/* ***************************************************************** */
/* some functions for Discontinous Galerkin (DG) method              */
/* ***************************************************************** */

/* NOTE: there is more in Legendre.c of titan:Archives/DG_test*/
