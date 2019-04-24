/* main/coordinates/integrals.c */
/* Wolfgang Tichy, 4/2019 */

#include "nmesh.h"
#include "coordinates.h"


enum { ASIS=0, ABS=1 };

/* compute volume integral \int dx dy dz v(x,y,z) of var v with
   index vind in a node. Here volume Jacobian is included. */
double NodeVolumeIntegral(tNode *node, int vind, double power, int mode)
{
  tPat *pat = node->pat;
  double *var = Vard(node,vind);
  double VolInt;

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

      /* integrand without Jacobian */
      Integ[i] = var[i];

      /* transform integrand */
      switch(mode)
      {
      case ABS:
        Integ[i] = fabs(Integ[i]);
      }
      Integ[i] = pow(Integ[i], power);

      /* include Jacobian in integrand */
      Integ[i] = Integ[i] * jac;
    }

    /* integrate with Jac. */
    VolInt = array_GLquadrature3(node, IntegA);

    free_array(IntegA);
  }

  return VolInt;
}
