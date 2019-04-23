/* norms.c */
/* Wolfgang Tichy, 4/2019 */

#include "nmesh.h"



/* compute volume integral \int dx dy dz v(x,y,z) of var v with
   index vind over a patch */
double PatchVolumeIntegral(tPat *pat, int vind)
{
  tMesh *mesh = pat->mesh;
  int myid;
  double sum, VolInt = 0.;

  formylnodes(mesh, myid)
  {
    tNode *node;

    if(node->pat != pat) continue;

    node = MyNode(mesh, myid);
    VolInt += NodeVolumeIntegral(node, vind);
  }

  nMPI_Allreduce(&VolInt, &sum, 1, nMPI_DOUBLE, nMPI_SUM);

  return sum;
}

/* compute volume integral \int dx dy dz v(x,y,z) of var v with
   index vind over entire mesh */
double MeshVolumeIntegral(tMesh *mesh, int vind)
{
  int myid;
  double sum, VolInt = 0.;

  formylnodes(mesh, myid)
  {
    tNode *node = MyNode(mesh, myid);
    VolInt += NodeVolumeIntegral(node, vind);
  }

  nMPI_Allreduce(&VolInt, &sum, 1, nMPI_DOUBLE, nMPI_SUM);

  return sum;
}
