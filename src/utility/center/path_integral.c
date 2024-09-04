
/* compute length of straight coor line from x1 to x2 */
double proper_length_of_coordline(const double x1[3], const double x2[3],
                                  tVarList *vl_metric)
{
  double m[3], c[3];
  int d, i;
  int np = 20;
  double *lambda = dmalloc(np);
  double *w      = dmalloc(np);
  tArray *xp[] = {alloc_array1d(np), alloc_array1d(np), alloc_array1d(np)};
  double *x[]  = {Arrd(xp[0]),       Arrd(xp[1]),       Arrd(xp[2])}


  /* fill xp: x^i(\lambda) = m^i \lambda +  c^i
     where: m^i = (x2^i - x1^i)/2
            c^i = (x2^i + x1^i)/2  */

  /* set m and c */
  for(d=0; d<3; d++)
  {
    m[d] = 0.5 * (x2[d] - x1[d]);
    c[d] = 0.5 * (x2[d] + x1[d]);
  }

  /* set lambda and w */
  LG_x_wquad(np, lambda, w);

  /* set xp, i.e. set x */
  for(d=0; d<3; d++)
    for(i=0; i<np; i++)
      x[d][i] = m[d] * lambda[i] + c[d];

  /* get metric at points xp */
  npts = 4;
  interp_VL_xp(mesh, vl_metric, xp, npts, INTERP_LAGRANGE, 1., val);

  /* set integrand f = sqrt(g_{ij} m^i m^j) */
  ...

  /* get path length */
  length = Gauss_integral(np, w, f);


  free_3_arrays(xp);
  free(w);
  free(lambda);

  return length;
}