/* basis.c */
/* Wolfgang Tichy, 2/2019 */

#include "nmesh.h"
#include "basis.h"



/***********************************************************************/
/* get derivatives of a variable */
/***********************************************************************/

/* get deriv in direction dir of array var, result goes into dvar */
void basis_array_deriv(tNode *node, int dir, tArray *var, tArray *dvar)
{
  tArray *Dt = node->Dt[dir];
  mm_array_indir(Dt, var, dir, dvar);
}
/* get deriv in direction dir of var with index vi, result goes to var dvi */
int basis_var_deriv(tNode *node, int dir, int vi, int dvi)
{
  tArray *var, *dvar;
  tDat *dat = node->dat;
  if(dat)
  {
    var  = dat->v[vi];
    dvar = dat->v[dvi];
    basis_array_deriv(node, dir, var, dvar);
    return 1;
  }
  else
    return 0;
}

/***********************************************************************/
/* functions for analysis and synthesis */
/***********************************************************************/

/* get coeffs in direction dir of array var, coeffs goes into c */
void basis_array_analysis(tNode *node, int dir, tArray *var, tArray *c)
{
  tArray *At = node->At[dir];
  mm_array_indir(At, var, dir, c);
}
/* get coeffs in direction dir of variable vi, coeffs goes into ci */
int basis_var_analysis(tNode *node, int dir, int vi, int ci)
{
  tArray *var, *c;
  tDat *dat = node->dat;
  if(dat)
  {
    var = dat->v[vi];
    c   = dat->v[ci];
    basis_array_analysis(node, dir, var, c);
    return 1;
  }
  else
    return 0;
}

/* get array from coeffs c in direction dir, array goes into var */
void basis_array_synthesis(tNode *node, int dir, tArray *var, tArray *c)
{
  tArray *St = node->St[dir];
  mm_array_indir(St, c, dir, var);
}
/* get var vi from coeffs in var ci in direction dir */
int basis_var_synthesis(tNode *node, int dir, int vi, int ci)
{
  tArray *var, *c;
  tDat *dat = node->dat;
  if(dat)
  {
    var = dat->v[vi];
    c   = dat->v[ci];
    basis_array_synthesis(node, dir, var, c);
    return 1;
  }
  else
    return 0;
}
