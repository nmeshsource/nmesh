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
  int m, j;
  double denom;

  for(j = 0; j < n; j++)
  {
    denom = 1.;
    for(m = 0; m < n; m++)
      if(m != j) denom *= x[j] - x[m];

    w_interp[j] = 1./denom;
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
