/* main/coordinates/integrals.c */
/* Wolfgang Tichy, 4/2019 */

#include "nmesh.h"
#include "coordinates.h"




/* compute volume integral \int dx dy dz v(x,y,z) of var v with
   index vind in a node. Here volume Jacobian is included. */
double NodeVolumeIntegral(tNode *node, int vind)
{
  tPat *pat = node->pat;
  double *var = Vard(node,vind);
  double VolInt;

  if(pat->dXYZ_dxyz) /* not Cartesian coords */
  {
    tMesh *mesh = pat->mesh;
    int i;
    double *det_dXbdx = Vard(node, Ind("det_dXbdx"));
    tArray *IntegA = alloc_array(node->n);
    double *Integ = Arrd(IntegA);

    forpoints(node, i)
    {
      double jac;
      double det = det_dXbdx[i];

      if(det!=0.0) jac = 1.0/fabs(det);
      /* if det=0 jac should really be infinite, but we hope that the integrand
         goes to zero quickly enough that jac=0 makes no difference! */
      else jac = 0.0;

      /* include Jacobian in integrand */
      Integ[i] = var[i] * jac;
    }
    /* integrate with Jac. */
    VolInt = array_GLquadrature3(node, IntegA);

    free_array(IntegA);
  }
  else /* integrate without jac */
    VolInt = var_GLquadratureXYZ3(node, vind);

  return VolInt;
}
