/* basis.c */
/* Wolfgang Tichy, 2/2019 */

#include "nmesh.h"
#include "basis.h"



/***********************************************************************/
/* get derivatives of a variable */
/***********************************************************************/

/* get deriv in direction dir of array var, result goes into dvar */
void basis_array_deriv1(tNode *node, int dir, tArray *var, tArray *dvar)
{
  tArray *Dt = node->Dt[dir];
  mm_array_indir(Dt, var, dir, dvar);
}
/* get deriv in direction dir of var with index vi, result goes to var dvi */
int basis_var_deriv1(tNode *node, int dir, int vi, int dvi)
{
  tArray *var, *dvar;
  tDat *dat = node->dat;
  if(dat)
  {
    var  = dat->v[vi];
    dvar = dat->v[dvi];
    basis_array_deriv1(node, dir, var, dvar);
    return 1;
  }
  else
    return 0;
}

/* compute 1st derivs in all 3 dirs, result goes into vars dvi, dvi+1, dvi+2 */
int basis_var_derivs(tNode *node, int vi, int dvi)
{
  tArray *var, *dvar;
  tDat *dat = node->dat;
  if(dat)
  {
    var  = dat->v[vi];
    int dir;
    for(dir=0; dir<3; dir++)
    {
      dvar = dat->v[dvi + dir];
      basis_array_deriv1(node, dir, var, dvar);
    }
    return 1;
  }
  else
    return 0;
}

/***********************************************************************/
/* functions for analysis and synthesis */
/***********************************************************************/

/* get coeffs in direction dir of array var, coeffs goes into c */
void basis_array_analysis1(tNode *node, int dir, tArray *var, tArray *c)
{
  tArray *At = node->At[dir];
  mm_array_indir(At, var, dir, c);
}
/* get coeffs in direction dir of variable vi, coeffs goes into ci */
int basis_var_analysis1(tNode *node, int dir, int vi, int ci)
{
  tArray *var, *c;
  tDat *dat = node->dat;
  if(dat)
  {
    var = dat->v[vi];
    c   = dat->v[ci];
    basis_array_analysis1(node, dir, var, c);
    return 1;
  }
  else
    return 0;
}

/* get array from coeffs c in direction dir, array goes into var */
void basis_array_synthesis1(tNode *node, int dir, tArray *var, tArray *c)
{
  tArray *St = node->St[dir];
  mm_array_indir(St, c, dir, var);
}
/* get var vi from coeffs in var ci in direction dir */
int basis_var_synthesis1(tNode *node, int dir, int vi, int ci)
{
  tArray *var, *c;
  tDat *dat = node->dat;
  if(dat)
  {
    var = dat->v[vi];
    c   = dat->v[ci];
    basis_array_synthesis1(node, dir, var, c);
    return 1;
  }
  else
    return 0;
}

/* get coeffs c=c_ijk of array u */
void basis_array_analysis3(tNode *node, tArray *u, tArray *c)
{
  basis_array_analysis1(node, 0, u, c);
  basis_array_analysis1(node, 1, c, c); //should work because it uses mm_array1
  basis_array_analysis1(node, 2, c, c);
}

/* get array u from coeffs c=c_ijk */
void basis_array_synthesis3(tNode *node, tArray *u, tArray *c)
{
  basis_array_synthesis1(node, 0, u, c);
  basis_array_synthesis1(node, 1, u, u); //should work because it uses mm_array1
  basis_array_synthesis1(node, 2, u, u);
}

/***********************************************************************/
/* interpolate */
/***********************************************************************/

/* 3d interpolation:
   interpolate to the point (Xb[0],Xb[1],Xb[2]) using coeffs in array coef */
double basis_array_interpolate(tNode *node, tArray *coef, double Xb[3])
{
  tPat *pat = node->pat;
  int *n = coef->n;
  double *B0 = dmalloc(n[0]);
  double *B1 = dmalloc(n[1]);
  double *B2 = dmalloc(n[2]);
  int k;
  double sum;

  /* save basis func values at (Xb[0],Xb[1],Xb[2]) in B0,... */
  for(k=0; k<n[0]; k++) B0[k] = pat->basis[0](k, Xb[0], n[0]);
  for(k=0; k<n[1]; k++) B1[k] = pat->basis[1](k, Xb[1], n[1]);
  for(k=0; k<n[2]; k++) B2[k] = pat->basis[2](k, Xb[2], n[2]);

  /* interpolate to (Xb[0],Xb[1],Xb[2]) */
  sum = 0.;
  //SGRID_LEVEL3_Pragma(omp parallel for reduction(+:sum))
  for(k = n[2]-1; k >=0; k--)
  {
    int j,i;
    for(j = n[1]-1; j >=0; j--)
    for(i = n[0]-1; i >=0; i--)
      sum += coef->d[Ind_n(i,j,k, n)] * B0[i] * B1[j] * B2[k];
  }

  free(B2);
  free(B1);
  free(B0);
  return sum;
}
