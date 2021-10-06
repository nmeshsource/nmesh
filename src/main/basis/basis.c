/* basis.c */
/* Wolfgang Tichy, 2/2019 */

#include "nmesh.h"
#include "basis.h"


#define PR 0

/* frequently used global vars */
tbasis basis[1];

/* initialize coordinates in each patch */
int basis_init_globals(tMesh *mesh)
{
  PRF;printf(":\n");

  /* set some global vars */
  basis->expfilter_JacobianPower = Par("basis_expfilter_JacobianPower");
  basis->filter_fv = Par("basis_filter_fv");

  return 0;
}

/***********************************************************************/
/* get derivatives of a variable */
/***********************************************************************/

/* get deriv in direction dir of array var, result goes into dvar */
void basis_array_deriv1(tNode *node, int dir, tArray *var, tArray *dvar,
                        tDerivOpt *opt)
{
  tArray *Dt;
  if(!opt)
  {
    Dt = node_Dt(node,dir);
  }
  else /* check options */
  {
    switch(opt->lop)
    {
    case -1:
      errorexit("backward");
      Dt = node_Dt(node,dir);
      break;
    case 1:
      errorexit("forward");
      Dt = node_Dt(node,dir);
      break;
    default:
      Dt = node_Dt(node,dir);
    }
  }
  /* use Dt diff matrix */
  mm_array_indir(Dt, var, dir, dvar);
}

/* compute 1st derivs of array var in all 3 dirs,
   result goes into arrays dvar[0], dvar[1], dvar[2] */
void basis_array_derivs(tNode *node, tArray *var, tArray *dvar[3],
                        tDerivOpt *opt)
{
  int dir;
  for(dir=0; dir<3; dir++)
    basis_array_deriv1(node, dir, var, dvar[dir], opt);
}

/* get deriv in direction dir of var with index vi, result goes to var dvi */
int basis_var_deriv1(tNode *node, int dir, int vi, int dvi, tDerivOpt *opt)
{
  tArray *var, *dvar;
  tDat *dat = node->dat;
  if(dat)
  {
    var  = dat->v[vi];
    dvar = dat->v[dvi];
    basis_array_deriv1(node, dir, var, dvar, opt);
    return 1;
  }
  else
    return 0;
}

/* compute 1st derivs in all 3 dirs,
   result goes into vars dvi[0], dvi[1], dvi[2] */
int basis_var_derivs(tNode *node, int vi, int dvi[3], tDerivOpt *opt)
{
  tArray *var, *dvar;
  tDat *dat = node->dat;
  if(dat)
  {
    var  = dat->v[vi];
    int dir;
    for(dir=0; dir<3; dir++)
    {
      dvar = dat->v[dvi[dir]];
      basis_array_deriv1(node, dir, var, dvar, opt);
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
  tArray *At = node_At(node,dir);
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
  tArray *St = node_St(node,dir);
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

/* get coeffs c=c_ijk of var ui into var ci */
int basis_var_analysis3(tNode *node, int ui, int ci)
{
  tArray *u, *c;
  tDat *dat = node->dat;
  if(dat)
  {
    u = dat->v[ui];
    c = dat->v[ci];
    basis_array_analysis3(node, u, c);
    return 1;
  }
  else
    return 0;
}

/* get var ui from coeffs in var ci */
int basis_var_synthesis3(tNode *node, int ui, int ci)
{
  tArray *u, *c;
  tDat *dat = node->dat;
  if(dat)
  {
    u = dat->v[ui];
    c = dat->v[ci];
    basis_array_synthesis3(node, u, c);
    return 1;
  }
  else
    return 0;
}

/***********************************************************************/
/* interpolate */
/***********************************************************************/

/* 3d interpolation:
   interpolate to the point (Xb[0],Xb[1],Xb[2]) using coeffs in array coef */
double basis_array_interpolate(tNode *node, tArray *coef, double Xb[3])
{
  int *n = coef->n;
  double *B0 = dmalloc(n[0]);
  double *B1 = dmalloc(n[1]);
  double *B2 = dmalloc(n[2]);
  int k;
  double sum;

  /* save basis func values at (Xb[0],Xb[1],Xb[2]) in B0,... */
  for(k=0; k<n[0]; k++) B0[k] = node_basis(node,0, k, Xb[0], n[0]);
  for(k=0; k<n[1]; k++) B1[k] = node_basis(node,1, k, Xb[1], n[1]);
  for(k=0; k<n[2]; k++) B2[k] = node_basis(node,2, k, Xb[2], n[2]);

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

/* 3d interpolation:
   interpolate var vi to the point (Xb[0],Xb[1],Xb[2]) locally */
double basis_var_interpolate_local(tNode *node, int vi, double Xb[3])
{
  tArray *v, *c;
  double val;

  /* set coeffs of var vi in c */
  v = VarA(node, vi);
  if(!v) return 0.; /* return 0 as interp value if var vi has no storage */
  c = alloc_array(node->n);
  basis_array_analysis3(node, v, c);

  /* interp var vi to Xb in node */
  val = basis_array_interpolate(node, c, Xb);
  free_array(c);
  return val;
}

/* 3d interpolation:
   call basis_var_interpolate_local and then send interp. val around */
double basis_var_interpolate(tNode *node, int vi, double Xb[3])
{
  double Val, val=0.;
  int Haveval, haveval=0;

  if(node->dat)
  {
    val = basis_var_interpolate_local(node, vi, Xb);
    haveval = 1;
  }
  Val = val;
  Haveval = haveval;

  /* find out how many have a value, and add all of them */
  nMPI_Allreduce(&haveval, &Haveval, 1, nMPI_INT, nMPI_SUM);
  nMPI_Allreduce(&val, &Val, 1, nMPI_DOUBLE, nMPI_SUM);
  if(!Haveval) errorexit("one MPI proc should have this node");
  Val = Val/Haveval;

  //PRF;printf(": Val=%g Haveval=%d\n", Val, Haveval);
  return Val;
}

/* 3d interpolation:
   interpolate var vi to the point (x[0],x[1],x[2]).
   out: val
   returns: node if success, or NULL if failure to find x */
tNode *basis_var_interpolate_mesh(tMesh *mesh, int vi, const double x[3],
                                  double *val)
{
  double X[3], Xb[3];
  tNode *node = node_XYZ_of_xyz_mesh(mesh, X, x);

  /* return NULL if node and X are not found */
  if(!node) return NULL;

  /* set Xb in node */
  XbYbZb_of_XYZ(node, Xb, X);

  /* interp var vi to Xb in node */
  *val = basis_var_interpolate(node, vi, Xb);

  return node;
}


/***********************************************************************/
/* integrate using Gauss-Lobatto points */
/***********************************************************************/

/* get integral in direction dir of array var, result goes into Ivar */
tArray *array_GLquadrature1(tNode *node, int dir, tArray *var, tArray *Ivar)
{
  tArray *Wq = node_Wq(node,dir);
  int Ivar_realloc;

  if(PR)
  {
    PRF;printf(": dir=%d\n", dir);
    printf("Wq");printarray(Wq);
    printf("var");printarray(var);
  }

  Ivar_realloc = redimension_array(Ivar, Arrn(var));
  if(Ivar_realloc) errorexit("Ivar was too small");

  /* multiply weights Wq and var */
  mm_array_indir(Wq, var, dir, Ivar);

  /* re-dim Ivar array to 1 in the direction we just integrated */
  Ivar_realloc = redim_array(Ivar, dir==0, dir==1, dir==2);
  if(Ivar_realloc) errorexit("Ivar was too small");

  return Ivar;
}

/* put 2d integral in directions perpendicular to norm into Ivar */
tArray *array_GLquadrature2(tNode *node, int norm, tArray *var, tArray *Ivar)
{
  switch(norm)
  {
  case 0:
    array_GLquadrature1(node, 1, var, Ivar);
    array_GLquadrature1(node, 2, Ivar, Ivar); //should work because it uses mm_array2
    break;
  case 1:
    array_GLquadrature1(node, 0, var, Ivar);
    array_GLquadrature1(node, 2, Ivar, Ivar); //should work because it uses mm_array2
    break;
  case 2:
    array_GLquadrature1(node, 0, var, Ivar);
    array_GLquadrature1(node, 1, Ivar, Ivar); //should work because it uses mm_array1
    break;
  default:
    errorexit("dir must be 0,1,2");
  }
  return Ivar;
}

/* compute 3d integral (\int d^3Xb var) of array var */
double array_GLquadrature3(tNode *node, tArray *var)
{
  double I;
  tArray *Ivar = alloc_array(node->n);

  array_GLquadrature1(node, 0, var, Ivar);
  array_GLquadrature1(node, 1, Ivar, Ivar); //should work because it uses mm_array1
  array_GLquadrature1(node, 2, Ivar, Ivar);
  I = Ivar->d[0];
  free_array(Ivar);

  return I;
}

/* compute average of array var */
double array_nodeaverage(tNode *node, tArray *var)
{
  /* in Xb coords node volume is 2*2*2=8, so divide by 8 */
  return array_GLquadrature3(node, var) * 0.125;
}

/* compute 3d integral (\int d^3X var) of array var */
double array_GLquadratureXYZ3(tNode *node, tArray *var)
{
  double dXdXb[3];
  dXYZ_dXbYbZb(node, dXdXb);
  return array_GLquadrature3(node, var) * fabs(dXdXb[0]*dXdXb[1]*dXdXb[2]);
}

/* compute 3d integral (\int d^3Xb u) of var ui */
double var_GLquadrature3(tNode *node, int ui)
{
  tArray *u;
  tDat *dat = node->dat;
  if(dat)
  {
    u = dat->v[ui];
    return array_GLquadrature3(node, u);
  }
  else
    errorexit("no dat on this node");
}

/* compute average of var ui */
double var_nodeaverage(tNode *node, int ui)
{
  /* in Xb coords node volume is 2*2*2=8, so divide by 8 */
  return var_GLquadrature3(node, ui) * 0.125;
}

/* compute 3d integral (\int d^3X u) of var ui */
double var_GLquadratureXYZ3(tNode *node, int ui)
{
  double dXdXb[3];
  dXYZ_dXbYbZb(node, dXdXb);
  return var_GLquadrature3(node, ui) * fabs(dXdXb[0]*dXdXb[1]*dXdXb[2]);
}
