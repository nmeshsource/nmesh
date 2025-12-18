/* center.c */
/* Wolfgang Tichy 1/2024 */

#include "nmesh.h"
#include "center.h"

/* global pars for center */
tcenter center[1];

/* globals from other modules */
extern tOutput output[1];


/* init global center vars */
int center_init_globals(tMesh *mesh)
{
  PRFs(":\n");
  center->cx[0][0] = Par("center0_x");
  center->cx[0][1] = Par("center0_y");
  center->cx[0][2] = Par("center0_z");
  center->cx[1][0] = Par("center1_x");
  center->cx[1][1] = Par("center1_y");
  center->cx[1][2] = Par("center1_z");
  center->cx[2][0] = Par("center2_x");
  center->cx[2][1] = Par("center2_y");
  center->cx[2][2] = Par("center2_z");

  /* make sure some pars are saved in checkpoint */
  if(!Getv(Par("center0_track"),"no"))
    checkpoint_save_pars_append(mesh, "center0_x center0_y center0_z "
                                      "center1_mass center2_mass");
  if(!Getv(Par("center1_track"),"no"))
    checkpoint_save_pars_append(mesh, "center1_x center1_y center1_z");
  if(!Getv(Par("center2_track"),"no"))
    checkpoint_save_pars_append(mesh, "center2_x center2_y center2_z");

  return 0;
}

/* update position of center N */
int centerN_update(tMesh *mesh, int N)
{
  double minmove = Getd(Par("center_track_minmove"));
  double m1 = Getd(Par("center1_mass"));
  double m2 = Getd(Par("center2_mass"));
  int track, var, findMax, dir, setCenter;
  double xold[3], xnew[3], x1[3], x2[3];
  double h;
  char pname[1000];

  /* read current center N location pars into xold and xnew */
  for(dir=0; dir<3; dir++)
    xold[dir] = xnew[dir] = Getd(center->cx[N][dir]);

  setCenter = 0;
  sprintf(pname, "center%d_track", N);
  track = Par(pname);
  if(Getv(track, "no"))
  {
    setCenter = 0;
  }
  else if(Getv(track, "max") || Getv(track, "min"))
  {
    findMax = Getv(track, "max") ?  1 : 0;
    sprintf(pname, "center%d_track_var", N);
    var = Ind( Gets(Par(pname)) );
    h = average_grid_spacing(mesh, xold);
    center_track_extremum(mesh, h, var, findMax, xold, minmove, xnew);
    setCenter = 1;
  }
  else if(Getv(track, "circle"))
  {
    double dx, d, omega, mo, M=m1+m2, t=mesh->time;
    static double d2 = -1.;
    if(d2<0.) //compute d2 only once to ensure both centers get same d2!
      for(d2=0., dir=0; dir<3; dir++)
      {
        x1[dir] = Getd(center->cx[1][dir]);
        x2[dir] = Getd(center->cx[2][dir]);
        dx = x1[dir] - x2[dir];
        d2 += dx*dx;
      }
    d = sqrt(d2);
    /* Kepler: M = omega^2 d^3 */
    omega = sqrt(M/(d*d2));
    mo = (N==1 ? m2 : -m1); /* mass of other object */
    xnew[0] = d*cos(omega*t) * mo/M; /* circ. orbit centered on 0 */
    xnew[1] = d*sin(omega*t) * mo/M;
    xnew[2] = 0.;
    setCenter = 1;
  }
  else if(Getv(track, "CM"))
  {
    /* track CM computed from centers 1 and 2 */
    for(dir=0; dir<3; dir++)
    {
      x1[dir] = Getd(center->cx[1][dir]);
      x2[dir] = Getd(center->cx[2][dir]);
      xnew[dir] = (m1*x1[dir] + m2*x2[dir])/(m1 + m2);
    }
    setCenter = 1;
  }
  else
  {
    PRF;printf(": %s = %s\n", pname, Gets(track));
    errorexit("parameter value unknown");
  }

  /* write xnew into current center N location pars */
  if(setCenter)
    for(dir=0; dir<3; dir++)
      Setd(center->cx[N][dir], xnew[dir]);

  return 0;
}

/* update positions of centers */
int center_update(tMesh *mesh)
{
  int N, dir;
  //PRFs(":\n");

  /* update center 1 and 2 */
  centerN_update(mesh, 1);
  centerN_update(mesh, 2);

  /* update center 0, which can depend on center 1 and 2 */
  centerN_update(mesh, 0);

  /* print center locations */
  if(Getb(Par("center_verbose")))
  {
    for(N=0; N<3; N++)
      for(dir=0; dir<3; dir++)
        printf("  center%d_%c = %g\n", N, 'x'+dir,
               Getd(center->cx[N][dir]));
  }

  /* set pt output vars if center position vars are not empty */
  for(N=0; N<min2(3, output->Noutpt); N++)
    for(dir=0; dir<3; dir++)
    {
      int centerN_xi = center->cx[N][dir];
      if(strlen( Gets(centerN_xi) ))
        output->xpt[N][dir] = Getd(centerN_xi);
    }

  return 0;
}


/* get grid spacing near point x */
double average_grid_spacing(tMesh *mesh, double x[3])
{
  tPat *pat;
  tElm0 elm0[1];
  ulong eid, elmindex;
  int elmrank;
  double X[3], H[3], h[3];
  int d, m,n;

  set_elm0_XYZ_of_xyz_mesh(mesh, elm0, &eid,&elmindex,&elmrank, X, x);
  for(d=0; d<3; d++)
  {
    if(elm0->n[d] > 1)
      H[d] = (elm0->bbox[2*d+1] - elm0->bbox[2*d])/(elm0->n[d] - 1);
    else
      H[d] = (elm0->bbox[2*d+1] - elm0->bbox[2*d])/elm0->n[d];
  }

  /* transform spacing H to Cartesian coords h if needed */
  pat = mesh->pat[Elm_p(elm0)];
  if(pat->dXYZ_dxyz)
  {
    double xx[3], dXdx[3][3], dxdX[3][3];

    pat->dXYZ_dxyz(pat,NULL,-1, X, xx, dXdx);
    inv3Dmat_from_3Dmat(dXdx, dxdX);

    for(m=0; m<3; m++)
      for(n=0; n<3; n++)
        h[m] = dxdX[m][n]*H[n];
  }
  else
  {
    for(m=0; m<3; m++) h[m] = H[m];
  }
  //printf("h[0],h[1],h[2]=%g %g %g\n", h[0],h[1],h[2]);

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
  /* center_extremum_step uses 3 points, so iord>3 may be unnecessary */
  int iord = Geti(Par("center_track_iorder"));
  int pr = 0;

  if(pr) { PRF;printf(": h=%g findMax=%d\n", h, findMax); }

  /* previous coordinates */
  x0 = xold[0];
  y0 = xold[1];
  z0 = xold[2];

  // value at the old center
  v0 = interp_var_x_y_z(mesh, var, x0,y0,z0, iord, INTERP_LAGRANGE, 1.);

  /* find var to left and right in x-dir */
  vm = interp_var_x_y_z(mesh, var, x0-h,y0,z0, iord, INTERP_LAGRANGE, 1.);
  vp = interp_var_x_y_z(mesh, var, x0+h,y0,z0, iord, INTERP_LAGRANGE, 1.);
  if(pr) printf("   %2.2e   %2.2e   %2.2e\n",vm,v0,vp);
  dx = center_extremum_step(vm,v0,vp, findMax, h);

  // y-direction
  vm = interp_var_x_y_z(mesh, var, x0,y0-h,z0, iord, INTERP_LAGRANGE, 1.);
  vp = interp_var_x_y_z(mesh, var, x0,y0+h,z0, iord, INTERP_LAGRANGE, 1.);
  if(pr) printf("   %2.2e   %2.2e   %2.2e\n",vm,v0,vp);
  dy = center_extremum_step(vm,v0,vp, findMax, h);

  // z- direction
  vm = interp_var_x_y_z(mesh, var, x0,y0,z0-h, iord, INTERP_LAGRANGE, 1.);
  vp = interp_var_x_y_z(mesh, var, x0,y0,z0+h, iord, INTERP_LAGRANGE, 1.);
  if(pr) printf("   %2.2e   %2.2e   %2.2e\n",vm,v0,vp);
  dz = center_extremum_step(vm,v0,vp, findMax, h);

  if(pr) printf("  dx=%g  dy=%g  dz=%g\n", dx,dy,dz);

  v = interp_var_x_y_z(mesh, var, x0+dx,y0+dy,z0+dz, iord, INTERP_LAGRANGE, 1.);
  if(pr) printf("  v(x0,y0,z0)=%g  v(x0+dx,y0+dy,z0+dz)=%g\n", v0, v);

  if( !finit(v) || !finit(v0) ||
      ((findMax) && (v<=v0))  ||
      ((!findMax) && (v>=v0)) ||
      (sqrt(dx*dx + dy*dy + dz*dz) < minmove*h) )
    dx = dy = dz = 0.;

  if(pr) printf("  ==> using: dx=%g  dy=%g  dz=%g\n", dx,dy,dz);

  x1 = x0 + dx;
  y1 = y0 + dy;
  z1 = z0 + dz;
  //if ((pintoxyplane == 1) || (pintoxyplane == 3)) z1 = z0 ;

  /* save the result */
  xnew[0] = x1;
  xnew[1] = y1;
  xnew[2] = z1;

  if(pr)
  {
    printf("  %.9f ->  %.9f\n", xold[0], xnew[0]);
    printf("  %.9f ->  %.9f\n", xold[1], xnew[1]);
    printf("  %.9f ->  %.9f\n", xold[2], xnew[2]);
  }

  return 1;
}
