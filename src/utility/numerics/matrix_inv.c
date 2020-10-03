/* matrix_inv.c,  Wolfgang Tichy 10/2020 */

#include "nmesh.h"
#include "numerics.h"

/* for gaussjordan_inv */
#define SWAP(a, b, temp) { (temp)=(a); (a)=(b); (b)=(temp); }


/* replace n*n matrix A_{jk} = a[j + k*n] by its inverse.
   similar to gaussj from numrec but with the RHS b removed, all numrec
   allocs are replaced, and array indices shifted by -1.
   numrec's a[j'][k'] = a[j + k*n] , where j'=j+1, k'=k+1 */
int gaussjordan_inv(int n, double a[])
{
  int *indxc = imalloc(n);
  int *indxr = imalloc(n);
  int *ipiv  = imalloc(n);
  int i,icol,irow,j,k,l,ll;
  double big,dum,pivinv,temp;

  for(j=0;j<n;j++) ipiv[j]=0;
  for(i=0;i<n;i++)
  {
    big=0.0;
    for(j=0;j<n;j++)
      if(ipiv[j] != 1)
        for(k=0;k<n;k++)
        {
          if(ipiv[k] == 0)
          {
            if(fabs(a[j+k*n]) >= big)
            {
              big=fabs(a[j+k*n]);
              irow=j;
              icol=k;
            }
          }
        }
    (ipiv[icol])++;

    if(irow != icol)
    {
      for(l=0;l<n;l++) SWAP(a[irow+l*n], a[icol+l*n], temp)
    }
    indxr[i]=irow;
    indxc[i]=icol;
    if(a[icol+icol*n] == 0.0) errorexit("Matrix inversion: Singular Matrix");
    pivinv=1.0/a[icol+icol*n];
    a[icol+icol*n]=1.0;
    for(l=0;l<n;l++) a[icol+l*n] *= pivinv;
    for(ll=0;ll<n;ll++)
    if(ll != icol)
    {
      dum=a[ll+icol*n];
      a[ll+icol*n]=0.0;
      for(l=0;l<n;l++) a[ll+l*n] -= a[icol+l*n]*dum;
    }
  }
  for(l=n-1;l>=0;l--)
  {
    if(indxr[l] != indxc[l])
    for(k=0;k<n;k++)
    SWAP(a[k+indxr[l]*n], a[k+indxc[l]*n], temp);
  }
  free(ipiv);
  free(indxr);
  free(indxc);
  return 0;
}

/* write inverse of n*n matrix M into Minv */
int Minv_from_M_gaussjordan(int n, const double M[], double Minv[])
{
  int i;
  for(i=0; i<n*n; i++) Minv[i] = M[i];
  return gaussjordan_inv(n, Minv);
}
