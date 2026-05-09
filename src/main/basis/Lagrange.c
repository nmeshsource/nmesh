/* Lagrange.c */
/* Wolfgang Tichy, 2/2018
   some functions to deal with Lagrange interpolation polynomials */

#include "nmesh.h"
#include "basis.h"


/* ************************************************************************ */
/* various functions needed for Lagrange interpolation                      */
/* ************************************************************************ */

/* get Lagrange interpolation weights w_interp from the n points in x,
   this is coming from the denominator in Lagrange interpolation only */
void Lagrange_winterp(int n, const double *x, double *w_interp)
{
  int m, i;
  double denom;

  for(i = 0; i < n; i++)
  {
    denom = 1.;
    for(m = 0; m < n; m++)
      if(m != i) denom *= (x[i] - x[m]);

    w_interp[i] = 1./denom;
  }
}

/* Standard barycentric 1d Lagrange interpolation of data f[k] to point x.
   Also called "second (true) form of the barycentric formula." */
#define BARY_DIFF_MIN 1e-200
double Lagrange_interp_barycentric2(double x, int n, const double *x_p,
                                    const double *w_interp, const double *f)
{
  double denom, numer, W, diff;
  int k;
  denom = numer = 0.;
  for(k=0; k<n; k++)
  {
    diff = x - x_p[k];
    if(fabs(diff) <= BARY_DIFF_MIN) return f[k];

    W = w_interp[k] / diff;
    numer += W * f[k];
    denom += W;
  }
  return numer / denom;
}
/* Standard barycentric 1d Lagrange interpolation of data f to point x,
   where data has a stride ds in memory, fscal is a dummy arg */
double Lagrange_interp_barycentric2_ds(double x, int n, const double *x_p,
                                       const double *w_interp,
                                       const double *f, int ds, double fscal)
{
  double denom, numer, W, diff;
  int k;

  if(ds==1) return Lagrange_interp_barycentric2(x, n, x_p, w_interp, f);
  /*
  PRF;printf(": n=%d", n);
  printf("  x_p=");
  for(k=0; k<n; k++) printf(" %g", x_p[k]);
  printf("\n");
  printf("  w_i=");
  for(k=0; k<n; k++) printf(" %g", w_interp[k]);
  printf("\n");
  printf("  f  =");
  for(k=0; k<n; k++) printf(" %g", f[ds*k]);
  printf("\n");
  */
  denom = numer = 0.;
  for(k=0; k<n; k++)
  {
    diff = x - x_p[k];
    if(fabs(diff) <= BARY_DIFF_MIN) return f[ds*k];

    W = w_interp[k] / diff;
    numer += W * f[ds*k];
    denom += W;
  }
  //PRF;printf(": x=%g  numer / denom = %g\n", x, numer / denom);
  return numer / denom;
}


/* WT's Notes:
   Def.  k_j := \Pi_{m\neq j}(x_j - x_m),  c_j := 1/k_j
   Then  l_j(x) := \Pi_{m\neq j}(x - x_m) / k_j = c_j \Pi_{m\neq j}(x - x_m)
   Note: l_j(x_i) = \delta_{ji}
   The deriv matrix is D_{ij} := d_x l_j(x_i).

   d_x l_j(x) = c_j \sum_{l\neq j} \Pi_{m\neq j,l}(x - x_m)  (*)

   If x \neq x_j (*) becomes
   d_x l_j(x) = c_j \sum_{l\neq j} [\Pi_{m\neq l}(x - x_m)] / (x - x_j)
              = [c_j / (x - x_j)]  \sum_{l\neq j} l_l(x) k_l
   Now set x = x_i \neq x_j
   d_x l_j(x_i) = [c_j / (x_i - x_j)] \sum_{l\neq j} \delta_{li} k_l
    =  D_{ij}   = c_j /[c_i (x_i - x_j)] = (c_j/c_i) / (x_i - x_j)

   If x = x_j we change (*) like this:
   d_x l_j(x) = c_j [\Pi_{m\neq j}(x - x_m)] \sum_{l\neq j} 1 / (x - x_l)
              = l_j(x) \sum_{l\neq j} 1 / (x - x_l)
   So if x = x_j
   D_{jj} = d_x l_j(x_j) = \sum_{l\neq j} 1 / (x_j - x_l)

   Note also
   \sum_{j\neq i} D_{ij} = \sum_{j\neq i} (c_j/c_i)/(x_i - x_j)
   Now define k_{ij} := \Pi_{m\neq i,j}(x_i - x_m)
   (c_j/c_i)/(x - x_j) = (k_i/k_j)/(x_i - x_j) = (k_i/k_j) k_{ij}/k_i
                       = k_{ij}/k_j = ? <-- Not 1/ (x_j - x_i) !!!

   On the other hand we know that we must have \sum_{j} D_{ij} = 0
   Thus D_{ii} = - \sum_{j\neq i} D_{ij}

   Why is \sum_{j} D_{ij} = 0? Well let's define f(x) := \sum_{j} l_j(x).
   Then f(x_i) = 1 for all i, since l_j(x_i) = \delta_{ji}.
   All l_j(x) are polynomials of degree N = n-1. So f(x) is a polynomial
   of degree n-1 that is 1 on n points. Thus f(x) = 1 everywhere.
   Hence d_x f(x) = d_x \sum_{j} l_j(x) = 0, and thus \sum_{j} D_{ij} = 0. */

/* find matrix D for taking derivatives: D_{ij} = \partial_x l_j(x_i),
   this sets the transpose D^T if DT is interpreted as stored in
   column-major form */
void Lagrange_DT(int n, const double *x, const double *w_interp, double *DT)
{
  int i, j;
  double Dii, Dij;

  for(i = 0; i < n; i++)
  {
    Dii = 0;
    for(j = 0; j < n; j++)
    {
      if(i != j)
      {
        Dij  = (w_interp[j] / w_interp[i]) / (x[i] - x[j]);
        Dii -= Dij;
        DT[j + i*n] = Dij;
        /* NOTE: this DT is D_{ij} in row-major form or its transpose
                 D_{ji} in column-major form */
      }
    }
    DT[i*n + i] = Dii;
  }
}

/* get Lagrange basis function l_k(x),
   here x_p are the grid points, w_interp the interp. weights */
double Lagrange_of_x(int k, double x, int np,
                     const double *x_p, const double *w_interp)
{
  int m;
  double prod = 1.;

  for(m=0; m<np; m++) if(m!=k) prod *= (x - x_p[m]);
  return prod * w_interp[k];
}

/* helpers: */
/* product where we omit (x - x_p[l]), i.e. product in in Lagrange_of_x */
double Lagrange_prod1(int l, double x, int np, const double *x_p)
{
  int q;
  double prod = 1.;

  for(q=0; q<np; q++) if(q!=l) prod *= (x - x_p[q]);
  return prod;
}

/* product where we omit (x - x_p[l]) and (x - x_p[m])*/
double Lagrange_prod2(int l, int m, double x, int np, const double *x_p)
{
  int q;
  double prod = 1.;

  for(q=0; q<np; q++) if(q!=l && q!=m) prod *= (x - x_p[q]);
  return prod;
}


/* Set matrix DT for finite differences using a stencil of size 2sr+1 on a
   uniform grid with n gridpoints in [-1,1]. sr is stencil radius */
void fd_deriv_DT_uniform__old(int n, const double *x, int sr, double *DT)
{
  int ssz = 2*sr + 1; /* stencil size */
  double *w_interp = malloc(ssz * sizeof(w_interp[0]));

  if(ssz<=0) errorexit("stencil size ssz needs to be possitive");

  /* for small n use Lagrange_DT for n points */
  if(n <= ssz)
  {
    Lagrange_winterp(n, x, w_interp);
    Lagrange_DT(n, x, w_interp, DT);
  }
  else /* use Lagrange_DT for ssz points */
  {
    double *xs = malloc(ssz * sizeof(xs[0]));
    double *Dt = malloc(ssz*ssz * sizeof(w_interp[0]));
    double h  = 2./(n-1);
    double hs = 2./(ssz-1);
    double fac = hs/h;
    int i,j, is,js;

    /* first zero DT */
    for(i=0; i<n; i++)
      for(j=0; j<n; j++)
         DT[i*n + j] = 0.;

    /* put ssz equally spaced points into xs */
    for(i=0; i<ssz; i++) xs[i] = -1. + hs*i;

    /* get diff mat Dt of size ssz*ssz */
    Lagrange_winterp(ssz, xs, w_interp);
    Lagrange_DT(ssz, xs, w_interp, Dt);

    /* loop over DT entries */
    for(i=0; i<n; i++)
    {
      int jmin;

      if(i<sr) /* left end */
      {
        is = i;
        jmin = 0;
      }
      else if(n-i <= sr) /* right end */
      {
        is = ssz - (n-i);
        jmin = n - ssz;
      }
      else /* middle */
      {
        is = sr;
        jmin = i - sr;
      }

      for(j=jmin; j<jmin+ssz; j++)
      {
        js = j - jmin;
        DT[i*n + j] = Dt[is*ssz + js] * fac;
      }
    }

    free(Dt);
    free(xs);
  }

  free(w_interp);
}

/* Set matrix DT for finite differences using a stencil of size ssz on a
   uniform grid with n gridpoints in [-1,1]. Shift stencil by lop to
   the right (forward differencing) or left if lop<0.
   Also set non-zero range of DT, in case range is non-trivial. If range
   is the trival [0,n) range is not set here. */
void fd_lopderiv_DT_uniform(int n, const double *x, int ssz, int lop,
                            double *DT, int *range[2])
{
  int odd = (ssz%2);  /* odd=1 if stencil size ssz is odd */
  int sszo2 = ssz/2;  /* stencil radius without forward or backward shift */
  int sdl, sdr;       /* stencil radius on left and right after shift */
  double *w_interp;

  if(ssz<=0) errorexit("stencil size ssz needs to be possitive");

  if(abs(lop)>sszo2)
    errorexiti("stencil size ssz is too small for lop=%d", lop);

  if(odd)
  {
    sdl = sszo2 - lop; /* stencil radius on the left after shift by lop */
    sdr = sszo2 + lop; /* stencil radius on the right after shift by lop */

  }
  else
  {
    int pos = (lop>0);
    if(lop==0) errorexit("for an even stencil size we need lop!=0");
    sdl = sszo2 - lop - !pos;
    sdr = sszo2 + lop -  pos;
  }

  /* mem for interpolation weights */
  w_interp = malloc(ssz * sizeof(w_interp[0]));

  /* for small n use Lagrange_DT for n points, here we cannot shift at all */
  if(n <= ssz)
  {
    Lagrange_winterp(n, x, w_interp);
    Lagrange_DT(n, x, w_interp, DT);
    /* range of DT is trivial [0,n). So we do not set range[0], range[1] */
  }
  else /* use Lagrange_DT for ssz points */
  {
    double *xs = malloc(ssz * sizeof(xs[0]));
    double *Dt = malloc(ssz*ssz * sizeof(w_interp[0]));
    double h  = 2./(n-1);
    double hs = 2./(ssz-1);
    double fac = hs/h;
    int i,j, is,js;

    /* first zero DT */
    for(i=0; i<n; i++)
      for(j=0; j<n; j++)
         DT[i*n + j] = 0.;

    /* put ssz equally spaced points into xs */
    for(i=0; i<ssz; i++) xs[i] = -1. + hs*i;

    /* get diff mat Dt of size ssz*ssz */
    Lagrange_winterp(ssz, xs, w_interp);
    Lagrange_DT(ssz, xs, w_interp, Dt);

    /* loop over DT entries */
    for(i=0; i<n; i++)
    {
      int jmin;

      if(i<sdl) /* left end */
      {
        is = i;
        jmin = 0;
      }
      else if(n-i <= sdr) /* right end */
      {
        is = ssz - (n-i);
        jmin = n - ssz;
      }
      else /* middle */
      {
        is = sdl;
        jmin = i - sdl;
      }

      /* set DT */
      for(j=jmin; j<jmin+ssz; j++)
      {
        js = j - jmin;
        DT[i*n + j] = Dt[is*ssz + js] * fac;
      }

      /* save j-range of DT */
      if(range[0]) range[0][i] = jmin;
      if(range[1]) range[1][i] = jmin+ssz;
    }

    free(Dt);
    free(xs);
  }

  free(w_interp);
}
