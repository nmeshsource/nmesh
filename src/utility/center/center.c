/* center.c */
/* Wolfgang Tichy 1/2024 */

#include "nmesh.h"
#include "center.h"

/* global pars for center */
tcenter center[1];


/* init global center vars */
int center_init_globals(tMesh *mesh)
{
  PRFs(":\n");
  center->center0_x = Par("center0_x");
  center->center0_y = Par("center0_y");
  center->center0_z = Par("center0_z");
  center->center1_x = Par("center1_x");
  center->center1_y = Par("center1_y");
  center->center1_z = Par("center1_z");
  center->center2_x = Par("center2_x");
  center->center2_y = Par("center2_y");
  center->center2_z = Par("center2_z");

  /* make sure some pars are saved in checkpoint */
  if(Geti(Par("center0_track_method")))
    checkpoint_save_pars_append(mesh, "center0_x center0_y center0_z");
  if(Geti(Par("center1_track_method")))
    checkpoint_save_pars_append(mesh, "center1_x center1_y center1_z");
  if(Geti(Par("center2_track_method")))
    checkpoint_save_pars_append(mesh, "center2_x center2_y center2_z");

  return 0;
}

/* update position of center N */
int centerN_update(tMesh *mesh, int N)
{
  double minmove = Getd(Par("center_track_minmove"));
  double m1 = Getd(Par("center1_mass"));
  double m2 = Getd(Par("center2_mass"));
  int meth, var, findMax, dir;
  double xold[3], xnew[3], x1[3], x2[3];
  double h;
  char pname[1000];

  /* read current center N location pars into xold and xnew */
  for(dir=0; dir<3; dir++)
    xold[dir] = xnew[dir] = Getd(center->center0_x + 3*N + dir);

  sprintf(pname, "center%d_track_method", N);
  meth = Geti(Par(pname));
  switch(meth)
  {
  case 0:
    break;
  case 1: /* track max */
  case 2: /* track min */
    sprintf(pname, "center%d_track_var", N);
    var = Ind( Gets(Par(pname)) );
    findMax = (meth==1) ?  1 : 0;
    h = average_grid_spacing(mesh, xold);
    center_track_extremum(mesh, h, var, findMax, xold, minmove, xnew);
    break;
  case 3: /* track CM computed from centers 1 and 2 */
    for(dir=0; dir<3; dir++)
    {
      x1[dir] = Getd(center->center0_x + 3*1 + dir);
      x2[dir] = Getd(center->center0_x + 3*2 + dir);
      xnew[dir] = (m1*x1[dir] + m2*x2[dir])/(m1 + m2);
    }
    break;
  default:
    errorexiti("unknown center track method", meth);
  }

  /* write xnew into current center N location pars */
  for(dir=0; dir<3; dir++)
    Setd(center->center0_x + 3*N + dir, xnew[dir]);

  return 0;
}

/* update positions of centers */
int center_update(tMesh *mesh)
{
  //PRFs(":\n");

  /* update center 1 and 2 */
  centerN_update(mesh, 1);
  centerN_update(mesh, 2);

  /* update center 0, which can depend on center 1 and 2 */
  centerN_update(mesh, 0);

  return 0;
}


/* get grid spacing near point x */
double average_grid_spacing(tMesh *mesh, double x[3])
{
  tElm0 elm0[1];
  ulong eid, elmindex;
  int elmrank;
  double X[3], h[3];
  int d;

  set_elm0_XYZ_of_xyz_mesh(mesh, elm0, &eid,&elmindex,&elmrank, X, x);
  for(d=0; d<3; d++)
    h[d] = (elm0->bbox[2*d+1] - elm0->bbox[2*d])/elm0->n[d];
  return max3(h[0],h[1],h[2]);
}

/* find maximum by fitting 1D polynomial:
  we use a 2nd order function to track extremum
  f = ax^2+bx+c
  we use f = [vm, v0, vp]  at  x = [-1, 0, 1]
  a = 0.5*(vp+vm) - v0;
  b = 0.5*(vp-vm);
  c = v0;
  f' == 0      =>  x = -b/(2.*a);
  we set ds = -b/(2.*a)
  the actual step towards the extremum is then delx = ds*h
  where h is the grid point spacing.
  Thus we return ds*h. */
double center_extremum_step(double vm, double v0, double vp, int findMax,
                            double h)
{
  double a = 0.5*(vp+vm) - v0;
  double ds;

  if(!finit(vm) || !finit(v0) || !finit(vp)) return 0.;

  if(fabs(a) <= 1e-10)
  {
    if      ((vp>v0) && (v0>vm)) ds = 1.;
    else if ((vm>v0) && (v0>vp)) ds =-1.;
    else                         ds = 0.;
    ds *= (findMax)?(1.):(-1.);
  }
  else
  {
    /* we step toward the place with f'=0, i.e. we ignore findMax */
    ds = 0.25*(vm-vp)/a;

    /* disallow large steps */
    ds = (ds> 1.)? 1.:ds;
    ds = (ds<-1.)?-1.:ds;
  }

  return ds*h;
}

/* track center by finding the approx location of the extremum position */
int center_track_extremum(tMesh *mesh, double h, int var, int findMax,
                          const double xold[3], double minmove,
                          double xnew[3])
{
  double v0,vm,vp, v;
  double dx,dy,dz;
  double x0,y0,z0, x1,y1,z1;
  int pr = 0;

  /* previous coordinates */
  x0 = xold[0];
  y0 = xold[1];
  z0 = xold[2];

  // value at the old puncture
  v0 = basis_var_interp_x_y_z(mesh, var, x0,y0,z0, Lagrange_of_x);

  /* find var to left and right in x-dir */
  vm = basis_var_interp_x_y_z(mesh, var, x0-h,y0,z0, Lagrange_of_x);
  vp = basis_var_interp_x_y_z(mesh, var, x0+h,y0,z0, Lagrange_of_x);
  if (pr) printf("   %2.2e   %2.2e   %2.2e\n",vm,v0,vp);
  dx = center_extremum_step(vm,v0,vp, findMax, h);

  // y-direction
  vm = basis_var_interp_x_y_z(mesh, var, x0,y0-h,z0, Lagrange_of_x);
  vp = basis_var_interp_x_y_z(mesh, var, x0,y0+h,z0, Lagrange_of_x);
  if (pr) printf("   %2.2e   %2.2e   %2.2e\n",vm,v0,vp);
  dy = center_extremum_step(vm,v0,vp, findMax, h);

  // z- direction
  vm = basis_var_interp_x_y_z(mesh, var, x0,y0,z0-h, Lagrange_of_x);
  vp = basis_var_interp_x_y_z(mesh, var, x0,y0,z0+h, Lagrange_of_x);
  if (pr) printf("   %2.2e   %2.2e   %2.2e\n",vm,v0,vp);
  dz = center_extremum_step(vm,v0,vp, findMax, h);

  if (pr) printf("center_track_extremum:  %e %e %e\n",dx,dy,dz);

  v = basis_var_interp_x_y_z(mesh, var, x0+dx,y0+dy,z0+dz, Lagrange_of_x);
  if( !finit(v) || !finit(v0) ||
      ((findMax) && (v<=v0))  ||
      ((!findMax) && (v>=v0)) ||
      (sqrt(dx*dx + dy*dy + dz*dz) < minmove*h) )
    dx = dy = dz = 0.;

  x1 = x0 + dx;
  y1 = y0 + dy;
  z1 = z0 + dz;
  //if ((pintoxyplane == 1) || (pintoxyplane == 3)) z1 = z0 ;

  /* save the result */
  xnew[0] = x1;
  xnew[1] = y1;
  xnew[2] = z1;

  if(0)
  {
    PRFs("\n");
    printf("  %.9f ->  %.9f\n", xold[0], xnew[0]);
    printf("  %.9f ->  %.9f\n", xold[1], xnew[1]);
    printf("  %.9f ->  %.9f\n", xold[2], xnew[2]);
  }

  return 1;
}
