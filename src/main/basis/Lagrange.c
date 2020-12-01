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
