/* Legendre.c */
/* Wolfgang Tichy, 2/2018
   some functions to deal with Legendre polynomials, especially at the
   Legendre Gauss-Lobatto (LGL) points. */

#include "nmesh.h"
#include "basis.h"



/* Legendre polynomial basis function, np is number of points (not need for
   basis_LegendreP, unless we normalize as in basis_normLegendreP) */
/* use recursion recommended in numrec to get P_l */
double basis_LegendreP(int l, double x, int np)
{
  double Pk=0., P0, P1;
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
  double ooci, oosci;

  ooci = 0.5*(2*i+1);                   /* ci = 2.0/(2*i+1); */
  if( (i==N) && (N!=0) ) ooci = 0.5*i;  /* ci = 2.0/i; */
  oosci = sqrt(ooci);                   /* sci = sqrt(ci); */
  return basis_LegendreP(i, x, 0)*oosci;
}

/* compute the values of the Legendre polynomial P=P_l(x), dP=dP_l(x)/dx,
   Q(x) = P_{l+1}(x) - P_{l-1}(x), and of dQ=Q'(x) */
void Legendre_P_dP_Q_dQ(int l, double x,
                        double *P, double *dP, double *Q, double *dQ)
{
  int k;
  double P0 = 1., P1 = x,  dP0 = 0., dP1 = 1.;
  double Pk=0., dPk=0.;

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
  double I, P,dP, Q,dQ;

  if(k == 0)
  {
    I = 1. - x;
  }
  else
  {
    Legendre_P_dP_Q_dQ(k,x, &P,&dP, &Q,&dQ);
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

/* Note: All points and weights in this file are computed for the standard
   interval [-1,1]. If we have an interval [a,b] instead, we have to scale
   the standard points and weights like this:
     x^[a,b] = m*x^[-1,1] + (a+b)/2,    where m:=(b-a)/2
     w^[a,b]_k = m w^[-1,1]_k                              */

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
void LGL_x_wquad(int npoints, double *x, double *w)
{
  int N = npoints-1;
  int i,j;
  double P, dP, Q, dQ, y, dy;
  int Newton_itmax = 100;
  double Newton_rel_dytol = 1e-15;

  /* special case for just 1 point */
  if(N == 0)
  {
    x[N] = 0.;
    if(w) w[N] = 2.;
    return;
  }

  /* end points are local extrema in [-1,1] */
  x[0] = -1.;
  x[N] =  1.;
  if(w)
  {
    w[0] = 2.0/(N*(N+1.0));
    w[N] = w[0];
  }

  /* if we have only 2 points we are done */
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

      /* find y such that Q(y) = 0, (recall Q=0 <==> dP/dx=0),
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
    x[N-j] = -y;
    x[j]   =  y;

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
   For LGL nodes, the w[i] are the integration weights from LGL_x_wquad
   This is accurate for polynomials up to degree 2n-3 for LGL. */
double Gauss_integral(int n, const double *w, const double *f)
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
void LGL_AT_ST_matrices(int n, const double *x, const double *w,
                        double *AT, double *ST)
{
  int normalized = 1;
  int N = n-1;
  int i, j;
  double ooci, oosci;

  /* row-loop */
  for(i = 0; i < n; i++)
  {
    ooci = 0.5*(2*i+1);                   /* ci = 2.0/(2*i+1); */
    if( (i==N) && (N!=0) ) ooci = 0.5*i;  /* ci = 2.0/i; */
    oosci = sqrt(ooci);                   /* sci = sqrt(ci); */

    /* col-loop */
    for(j = 0; j < n; j++)
    {
      if(normalized)
      {
	/* transpose of synthesis matrix Shat_{ij} */
	/* ST is used in AT below */
	ST[i + j*n] = LegendreP(i, x[j]) * oosci;
	/* transpose of analysis matrix Ahat_{ij} */
	AT[j + i*n] = w[j] * ST[i + j*n];
      }
      else
      {
	ST[i + j*n] = LegendreP(i, x[j]);
	AT[j + i*n] = w[j] * ST[i + j*n] * ooci;
      }
    }
  }
}


/* ************************************************************************ */
/* various functions needed for Legendre Gauss-points or nodes      */
/* ************************************************************************ */

/* Legendre Gauss nodes x_i (zeros of P_np(x)) and integration weights w_i
   (see Gauss' formula on p. 887 of Milton Abramowitz and Irene A. Stegun)
   N = degree, so there are np = N+1 points */
void LG_x_wquad(int np, double *x, double *w)
{
  int N = np-1;
  double np2 = np*np;
  int i,j;
  double P, dP, Q, dQ, y, dy;
  int Newton_itmax = 100;
  double Newton_rel_dytol = 1e-15;

  /* special case for just 1 point */
  if(N == 0)
  {
    x[N] = 0.;
    if(w) w[N] = 2.;
    return;
  }

  /* get nodes in (-1,0] */
  for(j = 0; j <= N/2; j++)
  {
    if(j == N/2.0)  /* same as: if( (N%2==0 && j==N/2) ) */
    {
      y = 0.;
    }
    else
    {
      /* guess for node */
      y = -( 1. - 1./(8.*np2) + 1./(8.*np2*np) ) * cos((4*j+3)*PI/(4*np+2.));
      //printf("guess y=%g", y);

      /* find y such that P(y) = 0,
         use Newton method to refine guess from line above */
      for(i = 0; i < Newton_itmax; i++)
      {
        Legendre_P_dP_Q_dQ(np, y, &P, &dP, &Q, &dQ);
        dy = - P / dP;
        y += dy;
        if(fabs(dy) < fabs(Newton_rel_dytol * y)) break;
      }
      //printf(" => y=%g\n", y);
    }

    /* set x array for [0,1] and [-1,0] */
    x[N-j] = -y;
    x[j]   =  y;

    /* set dP */
    Legendre_P_dP_Q_dQ(np, y, &P, &dP, &Q, &dQ);
    //printf("      y=%g P=%g dP=%g\n", y, P, dP);

    if(w)
    {
      w[j] = 2.0/((1. - y*y) * dP*dP);
      w[N-j] = w[j];
    }
  }
}

/* Integral I over a sphere:
   I = \int_0^{2\pi} d\phi [ \int_0^{\pi} d\theta \sin(\theta) f(\theta,\phi) ]
     ~ (\pi/n) \sum_{j=0}^{2n-1} \sum_{i=0}^{n-1} w_i f(\theta_i,\phi_j)
   where we use 2n points in the phi-dir:
     \phi_j = \pi j/n,      j \in [0,2n-1]
   and n points in the theta-dir:
     \theta_i = acos(x_i),  i \in [0,n-1]
     (the x_i are the Legendre Gauss nodes with weights w_i)
   Note: f[i + j*n] = f(\theta_i,\phi_j,),
   We need to get w[i] via LG_x_wquad(n,x,w); */
double LG_integrate_f_theta_phi(int n, const double *w,
                                const double f[2*n * n])
{
  int i, j;
  double I = 0.;

  for(j=0; j<2*n; j++)
    for(i=0; i<n; i++)
      I += f[j*n + i] * w[i];

  I *= PI/n;
  return I;
}


/* ************************************************************************ */
/* various functions needed for equally spaced points or nodes              */
/* ************************************************************************ */

/* compute the Gaussian quadrature weights using Gaussian integration itself:
   w_j = \int_{-1}^{1} dx l_j(x) = \sum_i wg_i l_j(xg_i)
   where l_j(x) is the interpolating polynomial made from the
   grid points x_i: l_j(x) = Lagrange_of_x(j,x,...).
   We assume that the x_i are symmetric about x=0.
   NOTE: w_j = \int_{-1}^{1} dx l_j(x) may not be positive!!! */
void Gauss_wquad_from_symm_x(int npoints, const double *x, double *w)
{
  int N = npoints-1;
  /* Here wg_i and xg_i are LGL weights and points that we pick such that
     the numerical integral is exact. */
  int ng = 1 + (npoints+3)/2; // ng = (npoints+3)/2; is probably enough
  double *xg = dmalloc(ng);
  double *wg = dmalloc(ng);
  double *fg = dmalloc(ng);
  double *wI = dmalloc(npoints);
  int i, j;

  /* get points and weights for ng points */
  LGL_x_wquad(ng, xg, wg);

  /* get interpolation weights wI for l_j(x) for npoints */
  Lagrange_winterp(npoints, x, wI);

  /* now get w_j by integration */
  for(j = 0; j <= N/2; j++)
  {
    /* set fg_i = l_j(xg_i) */
    for(i=0; i<ng; i++)
      fg[i] = Lagrange_of_x(j, xg[i], npoints, x, wI);

    /* w_j = \int_{-1}^{1} dx l_j(x) */
    w[j] = Gauss_integral(ng, wg, fg);
    w[N-j] = w[j];
  }

  free(wI);
  free(fg);
  free(wg);
  free(xg);
}

/* grid points for equally spaced i.e. uniform points from -1 to +1. */
void uniform_x(int npoints, double *x)
{
  int N = npoints-1;
  int j;
  double y, h;

  /* special case for just 1 point */
  if(N == 0)
  {
    x[N] = 0.;
    return;
  }

  /* end points */
  x[0] = -1.;
  x[N] =  1.;

  /* if we have only 2 points we are done */
  if(N == 1)
    return;

  /* grid spacing */
  h = 2./N;

  /* get nodes in (-1,0] */
  for(j = 1; j <= N/2; j++)
  {
    if(j == N/2.0)  /* same as: if( (N%2==0 && j==N/2) ) */
      y = 0.;
    else
      y = -1. + j*h;

    /* set x array for [-1,0] and [0,1] */
    x[N-j] = -y;
    x[j]   = +y;
  }
}

/* Nodes x_i and Gauss quadrature weights w_i for equally spaced i.e.
   uniform points.
   N = degree, so there are N+1 points
   This is using Gaussian quad to get the weights. */
void uniform_x_wGaussquad(int npoints, double *x, double *w)
{
  int N = npoints-1;

  /* set uniform points */
  uniform_x(npoints, x);

  /* special case for just 1 point */
  if(N == 0)
  {
    if(w) w[N] = 2.;
    return;
  }

  /* if we have only 2 points we are done */
  if(N == 1)
  {
    if(w)
    {
      w[N] = w[0] = 2.0/(N*(N+1.0));
    }
    return;
  }

  /* now compute the weights using Gaussian integration
     w_j = \int_{-1}^{1} dx l_j(x) = \sum_i wg_i l_j(xg_i)
     where l_j(x) is the interpolating polynomial made from the
     grid points x_i: l_j(x) = Lagrange_of_x(j,x,...). */
  if(w)
    Gauss_wquad_from_symm_x(npoints, x, w);
}

/* uniform points and weights for trapezoidal integration rule */
void uniform_x_wTrapez(int npoints, double *x, double *w)
{
  int N = npoints-1;
  int j;

  /* set uniform points */
  uniform_x(npoints, x);

  /* special case for just 1 point */
  if(npoints == 1)
  {
    if(w) w[N] = 2.;
    return;
  }

  if(w)
  {
    w[N] = w[0] = 1./N;
    for(j=1; j<N; j++) w[j] = 2./N;
  }
}

/* first set ST using LGL_AT_ST_matrices, which gets the correct ST for
   for any grid point, then set AT simply as the inverse of ST */
void Legendre_AT_ST_matrices(int n, const double *x, const double *w,
                             double *AT, double *ST)
{
  LGL_AT_ST_matrices(n, x, w, AT, ST);
  M_to_Minv_gaussjordan(n, ST, AT);
}


/* ***************************************************************** */
/* some functions for Discontinous Galerkin (DG) method              */
/* ***************************************************************** */

/* NOTE: there is more in Legendre.c of titan:Archives/DG_test */
