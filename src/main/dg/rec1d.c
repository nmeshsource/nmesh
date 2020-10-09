/* rec1d.c */
/* Wolfgang Tichy, October 2020 */

/* several reconstruction methods in 1d */

#include "nmesh.h"
#include "dg.h"




/* Interpolate a field u to midpoint with index im.
   Here we interpolate in the positive direction (p) from the left of the
   midpoint to the midpoint to obtain umid_p. */
double rec1d_p_0(int n, const double *u, int im)
{
  return u[im]; // one sided 0-th order interpolation
}

/* Interpolate a field u to midpoint with index im.
   Here we interpolate in the negative direction (m) from the right of the
   midpoint to the midpoint to obtain umid_m. */
double rec1d_m_0(int n, const double *u, int im)
{
  return u[im+1]; // one sided 0-th order interpolation
}
