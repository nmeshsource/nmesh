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

    if(VarA(node, vind)) /* only count nodes where var vind has storage */
      VolInt += NodeVolumeIntegral(node, vind, power, mode);
  }
  sum = VolInt;
  nMPI_Allreduce(&VolInt, &sum, 1, nMPI_DOUBLE, nMPI_SUM);

  return sum;
}


/* Compute MPI-proc local extremum of var with index vind inside a sphere of
   radius r. The sphere center is in xc, given in x,y,z coords.
   If pat!=NULL we only look inside patch pat.
   If xc==NULL we do not check if extremum location is within sphere.
   Set Mnode, Mijk to the node and point-index with the extremum. */
double SphereExtremumLoc_local(tMesh *mesh, tPat *pat,
                               const double *xc, double r,
                               int vind, int findMax,
                               tNode **Mnode, int *Mijk)
{
  double extr;
  if(findMax) extr = -DBL_MAX;
  else        extr = +DBL_MAX;

  *Mnode = NULL;
  *Mijk  = 0;
  formylnodes_noomp(mesh)
  {
    tNode *node = MyLnode;
    double nextr;
    int found, ijk, d;

    if(pat && node->pat != pat) continue;

    if(VarA(node,vind) == NULL) continue;

    /* check if there is a new extremum */
    found = 0;
    if(findMax)
    {
      nextr = max_array(VarA(node,vind), &ijk);
      if(nextr > extr) found = 1;
    }
    else
    {
      nextr = min_array(VarA(node,vind), &ijk);
      if(nextr < extr) found = 1;
    }

    if(found)
    {
      if(xc)
      {
        double X[3], x[3], xs[3];

        /* check if new extremum loc is within sphere */
        /* get X and x for ijk */
        XYZ_of_ind(node, ijk, X);
        set_xyz(NULL, node, ijk, X, x);

        /* continue to next node if |xs| = |x - xc| > r */
        for(d=0; d<3; d++) xs[d] = x[d] - xc[d];
        if(magnitude_xyz(xs) > r) continue;
      }

      extr = nextr;
      *Mnode = node;
      *Mijk = ijk;
    }
  }
  return extr;
}

/* compute MPI-proc local max of var with index vind over a patch or mesh,
   set Mnode, Mijk to the node and point-index with the Max */
double MeshMaxLoc_local(tMesh *mesh, tPat *pat, int vind,
                        tNode **Mnode, int *Mijk)
{
  return SphereExtremumLoc_local(mesh, pat, NULL, DBL_MAX,
                                 vind, 1, Mnode, Mijk);
}


/* compute MPI-proc local min of var with index vind over a patch or mesh,
   set Mnode, Mijk to the node and point-index with the Min */
double MeshMinLoc_local(tMesh *mesh, tPat *pat, int vind,
                        tNode **Mnode, int *Mijk)
{
  return SphereExtremumLoc_local(mesh, pat, NULL, DBL_MAX,
                                 vind, 0, Mnode, Mijk);
}


/* compute max of var with index vind over a patch or mesh */
double MeshMax(tMesh *mesh, tPat *pat, int vind)
{
  double Max, max;
  tNode *Mnode;
  int Mijk;

  max = MeshMaxLoc_local(mesh, pat, vind, &Mnode, &Mijk);
  Max = max;
  nMPI_Allreduce(&max, &Max, 1, nMPI_DOUBLE, nMPI_MAX);

  return Max;
}


/* compute min of var with index vind over a patch or mesh */
double MeshMin(tMesh *mesh, tPat *pat, int vind)
{
  double Min, min;
  tNode *Mnode;
  int Mijk;

  min = MeshMinLoc_local(mesh, pat, vind, &Mnode, &Mijk);
  Min = min;
  nMPI_Allreduce(&min, &Min, 1, nMPI_DOUBLE, nMPI_MIN);

  return Min;
}


/* compute max/min of var with index vind over a patch or mesh
   input: mesh, pat, vind, findMax
   output: Mp, Mnodeloc, Mijk, MX[3], Mx[3]
           output has max/min location */
double MeshExtremumLoc__old(tMesh *mesh, tPat *pat, int vind, int findMax,
                       int *Mp, char Mnodeloc[104], int *Mijk,
                       double *MX, double *Mx)
{
  tNode *Mnode=NULL;
  double Xb[3];

  struct { /* extremum and rank where extr. is */
    double extr;
    int rank;
  } mr[1], Mr[1];

  struct Loc { /* location info */
    int p;
    char nodeloc[104]; /* node location string */
    int ijk;
    double X[3];
    double x[3];
  };

  union { /* union to convert Loc to char array */
    struct Loc loc[1];
    char bytes[sizeof(struct Loc)];
  } uloc[1];

  /* write local extr and rank into mr and Mr */
  if(findMax) mr->extr  = MeshMaxLoc_local(mesh, pat, vind, &Mnode, Mijk);
  else        mr->extr  = MeshMinLoc_local(mesh, pat, vind, &Mnode, Mijk);
  mr->rank = nMPI_rank();
  Mr->extr = mr->extr;
  Mr->rank = mr->rank;
  //printf("mr->extr=%g\n", mr->extr);
  //printf("Mnode=%p *Mijk=%d\n", Mnode, *Mijk);

  if(Mnode)
  {
    /* write local patch coords into MX and uloc, if we found a node */
    XbYbZb_of_ind(Mnode, *Mijk, Xb);
    XYZ_of_XbYbZb(Mnode, Xb, MX);
    set_xyz(NULL, Mnode, *Mijk, MX, Mx);
    uloc->loc->p = Mnode->pat->p;
    node_location_str(Mnode, uloc->loc->nodeloc, 103);

    /* write local results into uloc */
    uloc->loc->ijk  = *Mijk;
    uloc->loc->X[0] = MX[0];
    uloc->loc->X[1] = MX[1];
    uloc->loc->X[2] = MX[2];
    uloc->loc->x[0] = Mx[0];
    uloc->loc->x[1] = Mx[1];
    uloc->loc->x[2] = Mx[2];
  }
  else
  {
    /* If we can't find a node just set uloc to zero, since in that case
       another MPI proc must have found something... */
    memset(&(uloc[0]), 0, sizeof(uloc[0]));
  }

  /* get global extr and rank into Mr */
  if(findMax) nMPI_Allreduce(mr, Mr, 1, nMPI_DOUBLE_INT, nMPI_MAXLOC);
  else        nMPI_Allreduce(mr, Mr, 1, nMPI_DOUBLE_INT, nMPI_MINLOC);

  /* now we have rank and value in Mr,
     so broadcast local results from Mr->rank to all */
  nMPI_Bcast(&(uloc->bytes[0]), sizeof(struct Loc), nMPI_CHAR, Mr->rank);

  /* set location */
  *Mp = uloc->loc->p;
  strncpy(Mnodeloc, uloc->loc->nodeloc, 104);
  *Mijk = uloc->loc->ijk;
  MX[0] = uloc->loc->X[0];
  MX[1] = uloc->loc->X[1];
  MX[2] = uloc->loc->X[2];
  Mx[0] = uloc->loc->x[0];
  Mx[1] = uloc->loc->x[1];
  Mx[2] = uloc->loc->x[2];

  return Mr->extr;
}


/* compute max/min of var with index vind over a patch or mesh
   within sphere |x-xc|<=r
   input: mesh, pat, xc, r, vind, findMax
   output: Mp, Mnodeloc, Mijk, MX[3], Mx[3]
           output has max/min location */
double SphereExtremumLoc(tMesh *mesh, tPat *pat, const double *xc, double r,
                         int vind, int findMax, int *Mp, char Mnodeloc[104],
                         int *Mijk, double *MX, double *Mx)
{
  tNode *Mnode=NULL;

  struct { /* extremum and rank where extr. is */
    double extr;
    int rank;
  } mr[1], Mr[1];

  struct Loc { /* location info */
    int p;
    char nodeloc[104]; /* node location string */
    int ijk;
    double X[3];
    double x[3];
  };

  union { /* union to convert Loc to char array */
    struct Loc loc[1];
    char bytes[sizeof(struct Loc)];
  } uloc[1];

  /* write local extr and rank into mr and Mr */
  mr->extr = SphereExtremumLoc_local(mesh,pat, xc,r,
                                     vind, findMax, &Mnode, Mijk);
  if(Mnode)  mr->rank = nMPI_rank();
  else       mr->rank = -1; /* set rank to -1 in case nothing was found */
  Mr->extr = mr->extr;
  Mr->rank = mr->rank;
  //printf("mr->extr=%g\n", mr->extr);
  //printf("Mnode=%p *Mijk=%d\n", Mnode, *Mijk);

  if(Mnode)
  {
    /* write local patch coords into MX and uloc, if we found a node */
    XYZ_of_ind(Mnode, *Mijk, MX);
    set_xyz(NULL, Mnode, *Mijk, MX, Mx);
    uloc->loc->p = Mnode->pat->p;
    node_location_str(Mnode, uloc->loc->nodeloc, 103);

    /* write local results into uloc */
    uloc->loc->ijk  = *Mijk;
    uloc->loc->X[0] = MX[0];
    uloc->loc->X[1] = MX[1];
    uloc->loc->X[2] = MX[2];
    uloc->loc->x[0] = Mx[0];
    uloc->loc->x[1] = Mx[1];
    uloc->loc->x[2] = Mx[2];
  }
  else
  {
    /* If we can't find a node just set most of uloc to zero. In that case
       another MPI proc could have found something... */
    memset(&(uloc[0]), 0, sizeof(uloc[0]));
    uloc->loc->p = -1; /* invalid patch signals that no node was found */
  }

  /* get global extr and rank into Mr */
  if(findMax) nMPI_Allreduce(mr, Mr, 1, nMPI_DOUBLE_INT, nMPI_MAXLOC);
  else        nMPI_Allreduce(mr, Mr, 1, nMPI_DOUBLE_INT, nMPI_MINLOC);

  /* if we have rank and value in Mr,
     broadcast local results from Mr->rank to all */
  if(Mr->rank >= 0)
    nMPI_Bcast(&(uloc->bytes[0]), sizeof(struct Loc), nMPI_CHAR, Mr->rank);

  /* set location */
  *Mp = uloc->loc->p;
  strncpy(Mnodeloc, uloc->loc->nodeloc, 104);
  *Mijk = uloc->loc->ijk;
  MX[0] = uloc->loc->X[0];
  MX[1] = uloc->loc->X[1];
  MX[2] = uloc->loc->X[2];
  Mx[0] = uloc->loc->x[0];
  Mx[1] = uloc->loc->x[1];
  Mx[2] = uloc->loc->x[2];

  return Mr->extr;
}

/* compute max/min of var with index vind over a patch or mesh
   input: mesh, pat, vind, findMax
   output: Mp, Mnodeloc, Mijk, MX[3], Mx[3]
           output has max/min location */
double MeshExtremumLoc(tMesh *mesh, tPat *pat, int vind, int findMax,
                       int *Mp, char Mnodeloc[104], int *Mijk,
                       double *MX, double *Mx)
{
  /* we pass NULL for xc, so that any sphere is not checked */
  return SphereExtremumLoc(mesh,pat, NULL,DBL_MAX, vind, findMax,
                           Mp, Mnodeloc, Mijk, MX, Mx);
}


/* Compute CRC of all pars */
void crc_MeshPars(tMesh *mesh, uint64_t *crc, size_t *cnt)
{
  int i;
  for(i=0; i<mesh->npdb; i++)
  {
    char *name = MeshParGetName(mesh, i);
    char *val  = Gets(i);

    crc64_continue_counters(name, strlen(name), crc, cnt);
    crc64_continue_counters(val, strlen(val), crc, cnt);
  }
}

/* Compute CRC of all pats */
void crc_pats(tMesh *mesh, uint64_t *crc, size_t *cnt)
{
  int i;
  for(i=0; i<mesh->npats; i++)
  {
    tPat *pat = mesh->pat[i];
    tCoordInfo *CI = pat->CI;

    /* include everything before mesh in crc */
    crc64_continue_counters(pat, offsetof(tPat, mesh), crc, cnt);
    /* include everything before FSurf in crc */
    crc64_continue_counters(CI, offsetof(tCoordInfo, FSurf), crc, cnt);
  }
}

/* Compute CRC of elms on this rank */
void crc_elms_local(void *Mesh, uint64_t *crc, size_t *cnt)
{
  tMesh *mesh = Mesh;
  formyelms_noomp(mesh)
  {
    tElm *elm = MyElm;
    tElm0 *elm0;
    union { const tElm *elm; tElm0 *elm0; } e2e0;

    /* inlcude only elm header in crc */
    e2e0.elm = elm;
    elm0 = e2e0.elm0;
    crc64_continue_counters(elm0, sizeof(elm0), crc, cnt);
  }
}

/* Compute CRC of varlist on this rank */
void crc_VarList_local(void *Vlist, uint64_t *crc, size_t *cnt)
{
  tVarList *vl = Vlist;
  tMesh *mesh = vl->mesh;
  int vli;

  formyelms_noomp(mesh)
  {
    tElm *elm = MyElm;
    forvl(vl, vli)
    {
      tArray *arr = VarA(elm, Vind(vl, vli));
      int N = ArrN(arr);
      double *d = Arrd(arr);
      crc64_continue_counters(d, N, crc, cnt);
    }
  }
}

/* get global CRC for elms */
void crc_elms(tMesh *mesh, uint64_t *crc, size_t *cnt)
{
  crc64_0start_global(mesh, crc_elms_local, mesh, crc, cnt);
}

/* get global CRC for Varlist */
void crc_VarList(tVarList *vl, uint64_t *crc, size_t *cnt)
{
  tMesh *mesh = vl->mesh;
  crc64_0start_global(mesh, crc_VarList_local, vl, crc, cnt);
}
