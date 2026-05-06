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
  intList *il = intList_alloc();
  int k;

  PRF;printf(":\n");

  /* set some global vars */
  basis->expfilter_JacobianPower = Par("basis_expfilter_JacobianPower");
  basis->filter_fv = Par("basis_filter_fv");
  str_to_intList(Gets(Par("basis_boundary_interp_order")), " ", il);
  forList(il, k)
  {
    if(k>=N_WENO_BOUNDARY_ORDERS)
      errorexit("N_WENO_BOUNDARY_ORDERS is too small");
    basis->boundary_interp_order[k] = ListEntry(il, k);
  }

  /* print global vars */
  printf(" basis->expfilter_JacobianPower = par_%04d :  "
         "Gets(basis->expf...Power) = %s\n",
         basis->expfilter_JacobianPower, Gets(basis->expfilter_JacobianPower));
  printf(" basis->filter_fv = par_%04d :  Gets(basis->filter_fv) = %s\n",
         basis->filter_fv, Gets(basis->filter_fv));
  printf(" basis->boundary_interp_order =");
  for(k=0; k<N_WENO_BOUNDARY_ORDERS; k++)
  {
    int ord = basis->boundary_interp_order[k];
    if(ord>0) printf(" %d", ord);
  }
  printf("\n");

  intList_free(il);

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
      Dt = node_Dmt(node,dir);
      break;
    case 1:
      Dt = node_Dpt(node,dir);
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

/* get coeffs c=c_ijk of array u, using ana. matrices in At[3] */
void basis_array_analysis3_At(tArray *At[3], tArray *u, tArray *c)
{
  mm_array_indir(At[0], u, 0, c);
  mm_array_indir(At[1], c, 1, c); //should work because it uses mm_array1
  mm_array_indir(At[2], c, 2, c);
}

/* get array u from coeffs c=c_ijk, using syn. matrices in St[3] */
void basis_array_synthesis3_St(tArray *St[3], tArray *u, tArray *c)
{
  mm_array_indir(St[0], c, 0, u);
  mm_array_indir(St[1], u, 1, u); //should work because it uses mm_array1
  mm_array_indir(St[2], u, 2, u);
}

/***********************************************************************/
/* Interpolate using the basis funcs node_basis for each node.
   NOTE: These funcs may be inaccurate on uniform grids, since
   init_gridpoints seems to set basis_normLegendreP for all cases...  */
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

  PRF;printf(": WARNING: this func is deprecated!\n");

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

  errorexit("use interp_var_local instead");

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
   call basis_var_interpolate_local and then send interp. val around
   Returns: 1 if success
            0 if if all MPI procs have node=NULL */
int basis_var_interpolate_ok(tNode *node, int vi, double Xb[3],
                             double *vinterp)
{
  double Val, val=0.;
  int Haveval, haveval=0;

  errorexit("use interpolate_var_ok instead");

  if(node) if(node->dat)
  {
    val = basis_var_interpolate_local(node, vi, Xb);
    haveval = 1;
  }
  Val = val;
  Haveval = haveval;

  /* find out how many have a value, and add all of them */
  MCK( nMPI_Allreduce(&haveval, &Haveval, 1, nMPI_INT, nMPI_SUM) );
  MCK( nMPI_Allreduce(&val, &Val, 1, nMPI_DOUBLE, nMPI_SUM) );
  if(!Haveval) return 0; /* could not get interp on any node */
  Val = Val/Haveval;

  //PRF;printf(": Val=%g Haveval=%d\n", Val, Haveval);
  *vinterp = Val;
  return 1; /* got interp value */
}

/* 3d interpolation: call basis_var_interpolate_ok */
double basis_var_interpolate(tNode *node, int vi, double Xb[3])
{
  double Val;
  int Haveval = basis_var_interpolate_ok(node, vi, Xb, &Val);

  errorexit("use interpolate_var instead");

  if(!Haveval) errorexit("one MPI proc should have this node");

  //PRF;printf(": Val=%g Haveval=%d\n", Val, Haveval);
  return Val;
}

/* 3d interpolation:
   interpolate var vi to the point (x[0],x[1],x[2]).
   out: val
   returns: 1 if success, or 0 if failure to find x */
int basis_var_interpolate_mesh(tMesh *mesh, int vi, const double x[3],
                               double *val)
{
  int Haveval;
  double X[3], Xb[3];
  tNode *node = node_XYZ_of_xyz_mesh(mesh, X, x);

  errorexit("use interpolate_var_mesh instead");

  /* set Xb in node */
  if(node) XbYbZb_of_XYZ(node, Xb, X);

  /* interp var vi to Xb in node */
  Haveval = basis_var_interpolate_ok(node, vi, Xb, val);
  return Haveval;
}


/***********************************************************************/
/* integrate using Gauss-Lobatto points */
/***********************************************************************/

/* get integral in direction dir of array var, result goes into Ivar */
tArray *array_GLquadrature1(tNode *node, int dir, tArray *var, tArray *Ivar)
{
  tArray *Wq = node_Wq(node,dir);
  int Ivar_realloc;
  int nr[] = {-1,-1,-1};

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
  nr[dir] = 1;
  Ivar_realloc = redimension_array(Ivar, nr);
  if(Ivar_realloc) errorexit("Ivar was too small");

  return Ivar;
}

/* put 2d integral (\int d^2Xb var) in directions perpendicular to norm
   into Ivar */
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

/* put 2d integral (\int d^2X var) in directions perpendicular to norm
   into Ivar  */
tArray *array_GLquadrature2X(tNode *node, int norm, tArray *var, tArray *Ivar)
{
  int k;
  double dXdXb[3];
  double jac; /* Jacobian between Xb and X is const */
  dXYZ_dXbYbZb(node, dXdXb);
  switch(norm)
  {
  case 0:
    jac = fabs(dXdXb[1]*dXdXb[2]);
    break;
  case 1:
    jac = fabs(dXdXb[0]*dXdXb[2]);
    break;
  case 2:
    jac = fabs(dXdXb[0]*dXdXb[1]);
    break;
  default:
    errorexit("dir must be 0,1,2");
  }
  /* get \int d^2Xb var */
  array_GLquadrature2(node, norm, var, Ivar);
  /* and then multiply by Jacobian to get \int d^2X var */
  forarray(Ivar,k) Arrd_(Ivar)[k] *= jac;

  return Ivar;
}

/* compute 3d integral (\int d^3Xb var) of array var */
double array_GLquadrature3(tNode *node, tArray *var)
{
  double Integ;
  tArray *Ivar = alloc_array(node->n);

  array_GLquadrature1(node, 0, var, Ivar);
  array_GLquadrature1(node, 1, Ivar, Ivar); //should work because it uses mm_array1
  array_GLquadrature1(node, 2, Ivar, Ivar);
  Integ = Ivar->d[0];
  free_array(Ivar);

  return Integ;
}

/* compute average of array var */
double array_nodeaverage(tNode *node, tArray *var)
{
  /* in Xb coords node volume is 2*2*2=8, so divide by 8 */
  return array_GLquadrature3(node, var) * 0.125;
}

/* compute 3d integral (\int d^3X var) of array var */
double array_GLquadrature3X(tNode *node, tArray *var)
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
double var_GLquadrature3X(tNode *node, int ui)
{
  double dXdXb[3];
  dXYZ_dXbYbZb(node, dXdXb);
  return var_GLquadrature3(node, ui) * fabs(dXdXb[0]*dXdXb[1]*dXdXb[2]);
}


/***********************************************************************/
/* get sub arrays out of arrays, e.g. for lower order interpolation */
/***********************************************************************/

/* Find an index range of size n around Xb0 in direction dir.
   If CenterOnXb0=1, we center exactly on Xb0, and shrink the range if it
   wouldn't fit into the node. Otherwise, we just push the range inside the
   node.
   Out: i0 = start of range , ni = size of range */
void IndexRange_Xb0_get__old(tNode *node, int dir, double Xb0, int n,
                        int CenterOnXb0, int *i0, int *ni)
{
  double *Xb = node_Xb(node,dir)->d; /* get point coords in node */
  int nn = node->n[dir];
  int ind0, b0,b1, nb;

  /* the default for nb is n */
  nb = n;

  /* find node-point ind0 closest to Xb0 */
  ind0 = nearest_i0_of_Xb_indir(node, dir, Xb0);
  if(ind0==-1)      ind0 = 0;
  else if(ind0<-1)  ind0 = nn-1;

  /* move ind0 to the left of Xb0 if needed */
  if( (n%2 == 0) && (Xb[ind0] > Xb0) ) ind0--;

  /* find index range start b0 and end b1 */
  b0 = ind0 - (nb-1)/2;
  b1 = b0 + (nb-1);

  /* check if range with ends b0 and b1 fits */
  if(CenterOnXb0) /* shorten range if it does not fit */
  {
    if( (n%2) && (ind0<=0) ) /* on left end */
    {
      if(b0 < 0)   { b1 -= -b0-1;    b0 = 0; }
      if(b1 >= nn) { b0 += b1-nn+1;  b1 = nn-1; }
    }
    else if(ind0>=nn-1) /* on right end */
    {
      if(b0 < 0)   { b1 -= -b0;      b0 = 0; }
      if(b1 >= nn) { b0 += b1-nn;    b1 = nn-1; }
    }
    else /* more in the middle */
    {
      if(b0 < 0)   { b1 -= -b0;      b0 = 0; }
      if(b1 >= nn) { b0 += b1-nn+1;  b1 = nn-1; }
    }
    /* cut off pieces outside node */
    if(b0 >= nn) b0 = nn-1;
    if(b1 < 0)   b1 = 0;
  }
  else /* push range inside and shorten it, if it does not fit */
  {
    if(b0 < 0)   { b1 += -b0;      b0 = 0; }
    if(b1 >= nn) { b0 -= b1-nn+1;  b1 = nn-1; }
    /* cut off left piece if needed */
    if(b0 < 0) b0 = 0;
  }

  /* set i0 and ni */
  *i0 = b0;
  *ni = b1 - b0 + 1;
}


/* Consider n (e.g. 7) grid points marked by o:
   0       1       2       3       4       5       6    <-- gridpoint number
   o   |   o   |   o   |   o   |   o   |   o   |   o    <-- gridpoints
   | 0 | 1 | 2 | 3 | 4 | 5 | 6 | 7 | 8 | 9 | 10| 11|    <-- half zone number

   HalfZoneOf_left returns the half zone number x is in for the points xp[] */
int HalfZoneOf_left(int n, const double xp[], double x)
{
  int i, h;

  /* special case for n=1 */
  if(n<2)
  {
    if(x < xp[0]) return 0;
    else          return 1;
  }

  /* find zone */
  h = 0;
  for(i=0; i<n-1; i++)
  {
    double beg = xp[i];
    double end = xp[i+1];
    double mid = (beg + end)*0.5;
    if(x<mid) break;
    h++;
    if(x<end) break;
    h++;
  }
  /* check if we are at or beyond last zone */
  if(h >= (n-1)*2) h = (n-1)*2 - 1;
  return h;
}

/* same as HalfZoneOf_left but start from right end */
int HalfZoneOf_back(int n, const double xp[], double x)
{
  int i, h;

  /* special case for n=1 */
  if(n<2)
  {
    if(x < xp[0]) return 0;
    else          return 1;
  }

  /* find zone */
  h = (n-1)*2 - 1;
  for(i=n-1; i>0; i--)
  {
    double beg = xp[i-1];
    double end = xp[i];
    double mid = (beg + end)*0.5;
    if(x>mid) break;
    h--;
    if(x>beg) break;
    h--;
  }
  /* check if we are at or beyond first zone */
  if(h < 0) h = 0;
  return h;
}

/* return 0 in last zone, 1 in 2nd to last, etc. */
int HalfZoneOf_right(int n, const double xp[], double x)
{
  int h = HalfZoneOf_back(n, xp, x);
  return (n-1)*2 - 1 - h;
}

/* Find an index range of size n around Xb0 in direction dir.
   If CenterOnXb0=1, we center exactly on Xb0, and shrink the range if it
   wouldn't fit into the node. Otherwise, we just push the range inside the
   node.
   Out: i0 = start of range , ni = size of range */
int IndexRange_Xb0_get(tNode *node, int dir, double Xb0, int n,
                       int CenterOnXb0, int *i0, int *ni)
{
  double *Xb = node_Xb(node,dir)->d; /* get point coords in node */
  int nn = node->n[dir];
  int ind0, b0,b1, nb;
  int hr_or_hl_is_zero = 0;
  double diff;

  /* the default for nb is n */
  nb = n;

  /* find node-point ind0 closest to Xb0 */
  ind0 = nearest_i0_of_Xb_indir_diff(node, dir, Xb0, &diff);

  /* return range of size one if Xb0 is on a grid point */
  if(fabs(diff) < dequaleps) { *i0 = ind0;  *ni = 1;  return 0; }

  /* ensure ind0 is inside elm */
  if(ind0==-1)      ind0 = 0;
  else if(ind0<-1)  ind0 = nn-1;

  /* move ind0 to the left of Xb0 if needed */
  if( (n%2 == 0) && (Xb[ind0] > Xb0) ) ind0--;

  /* find index range start b0 and end b1 */
  b0 = ind0 - (nb-1)/2;
  b1 = b0 + (nb-1);

  /* push range inside and shorten it, if it does not fit */
  if(b0 < 0)   { b1 += -b0;      b0 = 0; }
  if(b1 >= nn) { b0 -= b1-nn+1;  b1 = nn-1; }
  /* cut off left piece if needed */
  if(b0 < 0) b0 = 0;

  if(CenterOnXb0) /* shorten range if it does not fit */
  {
    if(b0==0)
    {
      int hl = HalfZoneOf_left(nn, Xb, Xb0);
      int nl = basis->boundary_interp_order[hl];
      if(nl)
      {
        nb = b1 - b0 + 1;
        if(nb > nl) b1 -= nb-nl; //pull in right end
      }
      if(hl==0) hr_or_hl_is_zero |= 1;
    }

    if(b1==nn-1)
    {
      int hr = HalfZoneOf_right(nn, Xb, Xb0);
      int nr = basis->boundary_interp_order[hr];
      if(nr)
      {
        nb = b1 - b0 + 1;
        if(nb > nr) b0 += nb-nr; //pull in left end
      }
      if(hr==0) hr_or_hl_is_zero |= 2;
    }
  }

  /* set i0 and ni */
  *i0 = b0;
  *ni = b1 - b0 + 1;
  return hr_or_hl_is_zero;
}


/* Find an index box of size n[3] around Xb0.
   If CenterOnXb0=1 we center exactly on Xb0, and shrink the box if it
   wouldn't fit into the node. Otherwise, we just push the box inside the
   node.
   Out: b0 = lower corner of box, nb = size of box nb */
void IndexBox_Xb0_get(tNode *node, const double Xb0[3], const int n[3],
                      int CenterOnXb0[3], int b0[3], int nb[3])
{
  int d;
  for(d=0; d<3; d++)
    IndexRange_Xb0_get(node, d, Xb0[d], n[d], CenterOnXb0[d],
                       &(b0[d]), &(nb[d]));
}

/* Copy node point coordinates into Xb[3] */
void IndexBox_set_Xb(tNode *node, int b0[3], int nb[3], tArray *Xb[3])
{
  int d, k;

  /* copy node points nXb into Xb */
  for(d=0; d<3; d++)
  {
    redim_array(Xb[d], nb[d],1,1);
    for(k=0; k<nb[d]; k++) Xb[d]->d[k] = node_Xb(node,d)->d[k + b0[d]];
  }
}

/* Copy data from ar inside a box of size nb[0]*nb[1]*nb[2] with origin
   corner at ijk=(b0[0],b0[1],b0[2]). */
void IndexBox_fill_subarray(tArray *ar, int b0[3], int nb[3], tArray *subar)
{
  int *na = Arrn(ar);
  int i,j,k;
  errorexit("this funs is not tested!");

  redimension_array(subar, nb);

  /* copy data into subar */
  for(k = 0; k < nb[2]; k++)
  for(j = 0; j < nb[1]; j++)
  for(i = 0; i < nb[0]; i++)
  {
    int ia = i + b0[0];
    int ja = j + b0[1];
    int ka = k + b0[2];
    Arrd_(subar)[Ind_n(i,j,k, nb)] = Arrd_(ar)[Ind_n(ia,ja,ka, na)];
  }
}
