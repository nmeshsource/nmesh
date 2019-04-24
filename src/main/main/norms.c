/* norms.c */
/* Wolfgang Tichy, 4/2019 */

#include "nmesh.h"



/* compute volume integral  \int dx dy dz v(x,y,z)^power  of var v with
   index vind over a patch, or the entire mesh if pat=0 */
double MeshVolumeIntegral(tMesh *mesh, tPat *pat, int vind,
                          double power, int mode)
{
  double sum, VolInt = 0.;
  int myid;

  formylnodes(mesh, myid)
  {
    tNode *node = MyNode(mesh, myid);

    if(pat && node->pat != pat) continue;

    VolInt += NodeVolumeIntegral(node, vind, power, mode);

  }
  sum = VolInt;
  nMPI_Allreduce(&VolInt, &sum, 1, nMPI_DOUBLE, nMPI_SUM);

  return sum;
}


/* compute max of var with index vind over a patch or mesh */
double MeshMax(tMesh *mesh, tPat *pat, int vind)
{
  double Max, max = -1e300;
  int myid;

  formylnodes(mesh, myid)
  {
    tNode *node = MyNode(mesh, myid);
    double nmax;
    int ijk;

    if(pat && node->pat != pat) continue;

    nmax = max_array(VarA(node,vind), &ijk);
    if(nmax > max) max = nmax;
  }
  Max = max;
  nMPI_Allreduce(&max, &Max, 1, nMPI_DOUBLE, nMPI_MAX);

  return Max;
}


/* compute min of var with index vind over a patch or mesh */
double MeshMin(tMesh *mesh, tPat *pat, int vind)
{
  double Min, min = 1e300;
  int myid;

  formylnodes(mesh, myid)
  {
    tNode *node = MyNode(mesh, myid);
    double nmin;
    int ijk;

    if(pat && node->pat != pat) continue;

    nmin = min_array(VarA(node,vind), &ijk);
    if(nmin < min) min = nmin;
  }
  Min = min;
  nMPI_Allreduce(&min, &Min, 1, nMPI_DOUBLE, nMPI_MIN);

  return Min;
}
