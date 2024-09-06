/* center.c */
/* Wolfgang Tichy 9/2024 */

#include "nmesh.h"
#include "center.h"



/* global pars for center */
extern tcenter center[1];


/* compute coord distance between x1 to x2 */
double coord_distance(const double x1[3], const double x2[3])
{
  double length;
  int d;

  length = 0.;
  for(d=0; d<3; d++)
  {
    double dd = x2[d] - x1[d];
    length += dd*dd;
  }
  length = sqrt(length);
  return length;
}

/* compute proper length of straight coord line from x1 to x2:
   length = \int_{-1}^1 d\lambda \sqrt{g_{ij} dx^i/d\lambda dx^j/d\lambda}
   here   x^i(\lambda) = m^i \lambda +  c^i,  with \lambda \in [-1,1]
   where: m^i = (x2^i - x1^i)/2,  c^i = (x2^i + x1^i)/2
   We use np points for \lambda \in [-1,1], and use interp. of order iord
   to find the metric at these points. */
double proper_length_of_coordline(const double x1[3], const double x2[3],
                                  tVarList *vl_3metric, int iord, int np)
{
  tMesh *mesh = vl_3metric->mesh;
  double m[3], c[3];
  double *lambda = dmalloc(np);
  double *w      = dmalloc(np);
  tArray *xp[] = {alloc_array1d(np), alloc_array1d(np), alloc_array1d(np)};
  double *x[]  = {Arrd(xp[0]),       Arrd(xp[1]),       Arrd(xp[2])};
  tArray *val = alloc_array1d(np*VLn(vl_3metric));
  double *f = dmalloc(np);
  double length;
  int d,e, i;

  /* set m and c */
  for(d=0; d<3; d++)
  {
    m[d] = 0.5 * (x2[d] - x1[d]);
    c[d] = 0.5 * (x2[d] + x1[d]);
  }

  /* set lambdas and weights w for np points */
  LG_x_wquad(np, lambda, w);

  /* set xp, i.e. set x */
  for(d=0; d<3; d++)
    for(i=0; i<np; i++)
      x[d][i] = m[d] * lambda[i] + c[d];

  /* get metric at points xp, using Lagrange interp. of order iord */
  interp_VL_xp(mesh, vl_3metric, xp, iord, INTERP_LAGRANGE, 1., val);

  /* set integrand f = \sqrt{g_{de} m^d m^e} */
  for(i=0; i<np; i++)
  {
    f[i] = 0.;

    for(d=0; d<3; d++)
      for(e=0; e<3; e++)
      {
        int vl_index = index_symmmat(3, d,e);
        double g_de = Arrd(val)[i + np*vl_index];
        f[i] += g_de * m[d]*m[e];
      }
    if(f[i]<0.) {PRF;printf("f[i]=%g\n", f[i]);}
    f[i] = sqrt(f[i]);
  }

  /* PRFs(": ");
  printarray(val);
  pr3v("m", m);
  pr3v("c", c);
  pr3v("lambda", lambda);
  pr3v("w", w);
  pr3v("x[0]", x[0]);
  pr3v("x[1]", x[1]);
  pr3v("x[2]", x[2]);
  pr3v("val[0]", Arrd(val));
  pr3v("f", f);
  printf("\n"); */

  /* get path length */
  length = Gauss_integral(np, w, f);

  free(f);
  free_array(val);
  free_3_arrays(xp);
  free(w);
  free(lambda);

  return length;
}


/* output distance between center1 and center2 */
int center_1_2_distance_output(tMesh *mesh)
{
  int di = -1;
  double dt = Getd(Par("center_distance_output_time"));

  /* //test with particular const metric:
  tVarList *u = vlalloc(mesh);
  vlpushone(u, Ind("ADM_gxx"));
  vlsetconstant(u, 2.);
  vldropn(u, 1);
  vlpushone(u, Ind("ADM_gxy"));
  vlsetconstant(u, -0.1);
  vldropn(u, 1);
  vlpushone(u, Ind("ADM_gxz"));
  vlsetconstant(u, -0.2);
  vldropn(u, 1);
  vlpushone(u, Ind("ADM_gyy"));
  vlsetconstant(u, 3.);
  vldropn(u, 1);
  vlpushone(u, Ind("ADM_gyz"));
  vlsetconstant(u, 0.3);
  vldropn(u, 1);
  vlpushone(u, Ind("ADM_gzz"));
  vlsetconstant(u, 4.);
  vldropn(u, 1);
  vlfree(u);
  */

  /* write only if it's time */
  if(TimeIsAt_di_dt(mesh, di,dt))
  {
    double time = mesh->time;
    double x1[3], x2[3];
    tVarList *vl_3metric;
    double r1, r2, coord_dist, prop_dist;
    int iord, np, dir;

    /* get positions of center 1 & 2 */
    for(dir=0; dir<3; dir++)
    {
      x1[dir] = Getd(center->cx[1][dir]);
      x2[dir] = Getd(center->cx[2][dir]);
    }

    /* Move x1 and x2 away from centers along line connecting them.
       This is important for punctures where metric diverges. */
    r1 = Getd(Par("center_distance_radius1"));
    r2 = Getd(Par("center_distance_radius2"));
    coord_dist = coord_distance(x1, x2);
    for(dir=0; dir<3; dir++)
    {
      x1[dir] = x1[dir] + (x2[dir] - x1[dir]) * r1/coord_dist;
      x2[dir] = x2[dir] - (x2[dir] - x1[dir]) * r2/coord_dist;
    }
    //PRFs(": ");pr3v("x1", x1);pr3v("x2", x2);printf("\n");

    /* calc distance between x1 and x2 in several ways */
    coord_dist = coord_distance(x1, x2);
    vl_3metric = vlalloc(mesh);
    vlpush(vl_3metric, Ind(Gets(Par("center_distance_metric"))));
    iord = Geti(Par("center_distance_iorder"));
    np   = Geti(Par("center_distance_npoints"));
    prop_dist = proper_length_of_coordline(x1, x2, vl_3metric, iord, np);
    vlfree(vl_3metric);

    /* write distances into file */
    if(Rank0)
    {
      FILE *fp;
      char fname[1000];
      sprintf(fname, "%s/%s", Gets(Par("outdir")), "center_1_2_distance.t");
      //PRF;printf(": %s %.15g %.15g %.15g\n", fname, time, coord_dist, prop_dist);
      fp = fopen(fname, "a");
      if(time==0.) fprintf(fp, "# time coord_dist prop_dist\n");
      fprintf(fp, "%.15g %.15g %.15g\n", time, coord_dist, prop_dist);
      fclose(fp);
    }
  }

  return 0;
}
