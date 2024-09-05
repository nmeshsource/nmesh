/* center.c */
/* Wolfgang Tichy 9/2024 */

#include "nmesh.h"
#include "center.h"



/* compute coord distance between x1 to x2 */
double coord_distance(const double x1[3], const double x2[3])
{
  double m[3];
  double length;
  int d;

  length = 0.;
  for(d=0; d<3; d++)
  {
    double dd = x2[d] - x1[d];
    length += dd*dd
  }
  length = sqrt(length);
  return length
}

/* compute proper length of straight coord line from x1 to x2:
   length = \int_{-1}^1 d\lambda \sqrt{g_{ij} dx^i/d\lambda dx^j/d\lambda}
   here   x^i(\lambda) = m^i \lambda +  c^i,  with \lambda \in [-1,1]
   where: m^i = (x2^i - x1^i)/2,  c^i = (x2^i + x1^i)/2
   We use np points for \lambda \in [-1,1], and use interp. of order iord
   to find the metric at these points. */
double proper_length_of_coordline(const double x1[3], const double x2[3],
                                  tVarList *vl_3metric, int np, int iord)
{
  double m[3], c[3];
  double *lambda = dmalloc(np);
  double *w      = dmalloc(np);
  tArray *xp[] = {alloc_array1d(np), alloc_array1d(np), alloc_array1d(np)};
  double *x[]  = {Arrd(xp[0]),       Arrd(xp[1]),       Arrd(xp[2])}
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
        double g_de = Arrd(val)[i + np*vl_index]
        f[i] += g_de * m[d]*m[e];
      }
    if(f[i]<0.) {PRF;printf("f[i]=%g\n", f[i]);}
    f[i] = sqrt(f[i]);
  }

  /* get path length */
  length = Gauss_integral(np, w, f);

  free(f);
  free_array(val);
  free_3_arrays(xp);
  free(w);
  free(lambda);

  return length;
}
