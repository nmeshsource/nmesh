/* norms.c */
/* Wolfgang Tichy, 4/2019 */

#include "nmesh.h"



/* compute volume integral  \int dx dy dz v(x,y,z)^power  of var v with
   index vind over a patch, or the entire mesh if pat=0 */
double MeshVolumeIntegral(tMesh *mesh, tPat *pat, int vind,
                          double power, int mode)
{
  double sum, VolInt = 0.;

  formylnodes_noomp(mesh)
  {
    tNode *node = MyLnode;

    if(pat && node->pat != pat) continue;

    VolInt += NodeVolumeIntegral(node, vind, power, mode);

  }
  sum = VolInt;
  nMPI_Allreduce(&VolInt, &sum, 1, nMPI_DOUBLE, nMPI_SUM);

  return sum;
}


/* compute MPI-proc local max of var with index vind over a patch or mesh,
   set Mnode, Mijk to the node and point-index with the Max */
double MeshMaxLoc_local(tMesh *mesh, tPat *pat, int vind,
                        tNode *Mnode, int *Mijk)
{
  double max = -DBL_MAX;

  Mnode = NULL;
  formylnodes_noomp(mesh)
  {
    tNode *node = MyLnode;
    double nmax;
    int ijk;

    if(pat && node->pat != pat) continue;

    nmax = max_array(VarA(node,vind), &ijk);
    if(nmax > max)
    {
      max = nmax;
      Mnode = node;
      *Mijk = ijk;
    }
  }
  return max;
}


/* compute MPI-proc local min of var with index vind over a patch or mesh,
   set Mnode, Mijk to the node and point-index with the Min */
double MeshMinLoc_local(tMesh *mesh, tPat *pat, int vind,
                        tNode *Mnode, int *Mijk)
{
  double min = DBL_MAX;

  Mnode = NULL;
  formylnodes_noomp(mesh)
  {
    tNode *node = MyLnode;
    double nmin;
    int ijk;

    if(pat && node->pat != pat) continue;

    nmin = min_array(VarA(node,vind), &ijk);
    if(nmin < min)
    {
      min = nmin;
      Mnode = node;
      *Mijk = ijk;
    }
  }
  return min;
}


/* compute max of var with index vind over a patch or mesh */
double MeshMax(tMesh *mesh, tPat *pat, int vind)
{
  double Max, max;
  tNode *Mnode=NULL;
  int Mijk;

  max = MeshMaxLoc_local(mesh, pat, vind, Mnode, &Mijk);
  Max = max;
  nMPI_Allreduce(&max, &Max, 1, nMPI_DOUBLE, nMPI_MAX);

  return Max;
}


/* compute min of var with index vind over a patch or mesh */
double MeshMin(tMesh *mesh, tPat *pat, int vind)
{
  double Min, min;
  tNode *Mnode=NULL;
  int Mijk;

  min = MeshMinLoc_local(mesh, pat, vind, Mnode, &Mijk);
  Min = min;
  nMPI_Allreduce(&min, &Min, 1, nMPI_DOUBLE, nMPI_MIN);

  return Min;
}


/* compute max of var with index vind over a patch or mesh,
   set Mnode, Mijk to the node and point-index with the Max */
double MeshMaxLoc(tMesh *mesh, tPat *pat, int vind, tNode *Mnode, int *Mijk)
{
  struct { /* max and rank where max is */
    double max;
    int rank;
  } mr[1], Mr[1];

  struct Loc { /* location info */
    tNode *node;
    int ijk;
  };

  union { /* union to convert Loc to char array */
    struct Loc loc[1];
    char bytes[sizeof(struct Loc)];
  } uloc[1];

  /* get local max and rank into mr and Mr */
  mr->max  = MeshMaxLoc_local(mesh, pat, vind, Mnode, Mijk);
  mr->rank = nMPI_rank();
  Mr->max  = mr->max;
  Mr->rank = mr->rank;

  /* get global max and rank into Mr */
  nMPI_Allreduce(mr, Mr, 1, nMPI_DOUBLE_INT, nMPI_MAXLOC);

  /* now we know rank and value, so broadcast node and point index to all */
  nMPI_Bcast(&(uloc->bytes), sizeof(struct Loc), nMPI_CHAR, Mr->rank);

  /* set location */
  Mnode = uloc->loc->node;
  *Mijk = uloc->loc->ijk;

  return Mr->max;
}


/* compute min of var with index vind over a patch or mesh,
   set Mnode, Mijk to the node and point-index with the Min */
double MeshMinLoc(tMesh *mesh, tPat *pat, int vind, tNode *Mnode, int *Mijk)
{
  struct { /* min and rank where min is */
    double min;
    int rank;
  } mr[1], Mr[1]; /* for use as MPI_DOUBLE_INT */

  struct Loc { /* location info */
    tNode *node;
    int ijk;
  };

  union { /* union to convert struct Loc to char array */
    struct Loc loc[1];
    char bytes[sizeof(struct Loc)];
  } uloc[1];

  /* get local min and rank into mr and Mr */
  mr->min  = MeshMinLoc_local(mesh, pat, vind, Mnode, Mijk);
  mr->rank = nMPI_rank();
  Mr->min  = mr->min;
  Mr->rank = mr->rank;

  /* get global min and rank into Mr */
  nMPI_Allreduce(mr, Mr, 1, nMPI_DOUBLE_INT, nMPI_MINLOC);

  /* now we know rank and value, so broadcast node and point index to all */
  nMPI_Bcast(&(uloc->bytes), sizeof(struct Loc), nMPI_CHAR, Mr->rank);

  /* set location */
  Mnode = uloc->loc->node;
  *Mijk = uloc->loc->ijk;

  return Mr->min;
}
