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
   her x_p are the grid points, w_interp the interp. weights */
double Lagrange_of_x(int k, double x, int np,
                     const double *x_p, const double *w_interp)
{
  int m;
  double prod = 1.;

  for(m=0; m<np; m++) if(m!=k) prod *= (x - x_p[m]);
  return prod * w_interp[k];
}

/***********************************************************************/
/* interpolate */
/***********************************************************************/

/* 3d interpolation:
   interpolate to the point (Xb[0],Xb[1],Xb[2]) for variable in array var
   Note: for Lagrange interpolation the coeffs are simply the function
         values at grid points */
double Lagrange_array_interpolate(tNode *node, tArray *var, double Xb[3])
{
  int *n = node->n;
  double *xp0 = node->Xb[0]->d; /* points */
  double *xp1 = node->Xb[1]->d;
  double *xp2 = node->Xb[2]->d;
  double *w0 = node->WL[0]->d;  /* weights */
  double *w1 = node->WL[1]->d;
  double *w2 = node->WL[2]->d;
  double *restrict B0 = dmalloc(n[0]);   /* basis */
  double *restrict B1 = dmalloc(n[1]);
  double *restrict B2 = dmalloc(n[2]);
  int k;
  double sum;

  /* save basis func values at (Xb[0],Xb[1],Xb[2]) in B0,... */
  for(k=0; k<n[0]; k++) B0[k] = Lagrange_of_x(k, Xb[0], n[0], xp0, w0);
  for(k=0; k<n[1]; k++) B1[k] = Lagrange_of_x(k, Xb[1], n[1], xp1, w1);
  for(k=0; k<n[2]; k++) B2[k] = Lagrange_of_x(k, Xb[2], n[2], xp2, w2);

  /* interpolate to (Xb[0],Xb[1],Xb[2]) */
  sum = 0.;
  //SGRID_LEVEL3_Pragma(omp parallel for reduction(+:sum))
  for(k = n[2]-1; k >=0; k--)
  {
    int j,i;
    for(j = n[1]-1; j >=0; j--)
    for(i = n[0]-1; i >=0; i--)
      sum += var->d[Ind_n(i,j,k, n)] * B0[i] * B1[j] * B2[k];
  }

  free(B2);
  free(B1);
  free(B0);
  return sum;
}

/* 2d interpolation:
   interpolate to the point (Cb1, Cb2) for variable in array var
   in plane p orthogonal to direction dir
   NOTE: We can set node=neighbor when we call this, even if the var is not
         on neighbor. We can use this to interpolate on a surface that was
         copied from a neighbor node! */
double Lagrange_array_interpolate2d(tNode *node, tArray *var,
                                    int dir, int p, double Cb1, double Cb2)
{
  int *n = node->n;
  double *xp0 = node->Xb[0]->d; /* points */
  double *xp1 = node->Xb[1]->d;
  double *xp2 = node->Xb[2]->d;
  double *w0 = node->WL[0]->d;  /* weights */
  double *w1 = node->WL[1]->d;
  double *w2 = node->WL[2]->d;
  double *restrict B1 = dmalloc(n[1]);
  double *restrict B2 = dmalloc(n[2]);
  int i,j,k;
  double sum;

  switch(dir)
  {
  case 1:
    /* save basis func values */
    for(k=0; k<n[1]; k++) B1[k] = Lagrange_of_x(k, Cb1, n[1], xp1, w1);
    for(k=0; k<n[2]; k++) B2[k] = Lagrange_of_x(k, Cb2, n[2], xp2, w2);

    /* interpolate */
    sum = 0.;
    for(k = n[2]-1; k >=0; k--)
    for(j = n[1]-1; j >=0; j--)
      sum += var->d[Ind_n(p,j,k, n)] * B1[j] * B2[k];
    break;
  case 2:
    /* save basis func values */
    for(k=0; k<n[0]; k++) B1[k] = Lagrange_of_x(k, Cb1, n[0], xp0, w0);
    for(k=0; k<n[2]; k++) B2[k] = Lagrange_of_x(k, Cb2, n[2], xp2, w2);

    /* interpolate */
    sum = 0.;
    for(k = n[2]-1; k >=0; k--)
    for(i = n[0]-1; i >=0; i--)
      sum += var->d[Ind_n(i,p,k, n)] * B1[i] * B2[k];
    break;
  case 3:
    /* save basis func values */
    for(k=0; k<n[0]; k++) B1[k] = Lagrange_of_x(k, Cb1, n[0], xp0, w0);
    for(k=0; k<n[1]; k++) B2[k] = Lagrange_of_x(k, Cb2, n[1], xp1, w1);

    /* interpolate */
    sum = 0.;
    for(j = n[1]-1; j >=0; j--)
    for(i = n[0]-1; i >=0; i--)
      sum += var->d[Ind_n(i,j,p, n)] * B1[i] * B2[j];
    break;
  default:
    errorexit("dir must be 0,1,2");
  }
  free(B2);
  free(B1);
  return sum;
}
