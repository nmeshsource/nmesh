/* grid.c */
/* Wolfgang Tichy, 10/2019 */

#include "nmesh.h"
#include "basis.h"

#define PR 0


/* global grid point structure */
tGridPoints gridpoints[1];

/* init gridpoints structure */
int init_gridpoints(tMesh *mesh)
{
  int nmax = Geti(Par("amr_nmax"));
  int ni, typ;

  /* save nmax */
  gridpoints->nmax = nmax;

  /* get mem. for grid points, diff. matrices, ... */
  gridpoints->Xb = calloc(nmax+1, sizeof(gridpoints->Xb[0]));
  if(!(gridpoints->Xb) )
    errorexit("out of memory for points");
  gridpoints->Wq = calloc(nmax+1, sizeof(gridpoints->Wq[0]));
  if(!(gridpoints->Wq) )
    errorexit("out of memory for integr. weights");

  gridpoints->WL = calloc(nmax+1, sizeof(gridpoints->WL[0]));
  if(!(gridpoints->WL) )
    errorexit("out of memory for Lagrange interp. weights");

  gridpoints->Dt = calloc(nmax+1, sizeof(gridpoints->Dt[0]));
  if(!(gridpoints->Dt) )
    errorexit("out of memory for diff. matrices");

  gridpoints->At = calloc(nmax+1, sizeof(gridpoints->At[0]));
  if(!(gridpoints->At) )
    errorexit("out of memory for ana. matrices");
  gridpoints->St = calloc(nmax+1, sizeof(gridpoints->St[0]));
  if(!(gridpoints->St) )
    errorexit("out of memory for syn. matrices");

  /* allocate arrays */
  for(ni=1; ni<=nmax; ni++)
  {
    int n[3];

    n[0] = n[1] = ni;
    n[2] = 1;
    for(typ=0; typ<P_NTYPES; typ++)
    {
      gridpoints->Dt[ni][typ] = alloc_array(n);
      gridpoints->At[ni][typ] = alloc_array(n);
      gridpoints->St[ni][typ] = alloc_array(n);
    }
    n[0] = ni;
    n[1] = n[2] = 1;
    for(typ=0; typ<3; typ++)
    {
      gridpoints->Xb[ni][typ] = alloc_array(n);
      gridpoints->Wq[ni][typ] = alloc_array(n);
      gridpoints->WL[ni][typ] = alloc_array(n);
    }
  }

  /* set points, weights, diff, and other matrices */
  for(typ=0; typ<P_NTYPES; typ++)
  {
    for(ni=1; ni<=nmax; ni++)
    {
      double *Xb = gridpoints->Xb[ni][typ]->d;
      double *Wq = gridpoints->Wq[ni][typ]->d;
      double *WL = gridpoints->WL[ni][typ]->d;
      double *DT = gridpoints->Dt[ni][typ]->d;
      double *AT = gridpoints->At[ni][typ]->d;
      double *ST = gridpoints->St[ni][typ]->d;

      /* set up desired points */
      switch(typ)
      {
      case P_UNIFORM: /* set equally spaced points and their weights */
        //uniform_x_wGaussquad(ni, Xb, Wq);
        uniform_x_wTrapez(ni, Xb, Wq);

        /* get analysis & synthesis matrix for Legendre basis */
        Legendre_AT_ST_matrices(ni, Xb, Wq, AT, ST);
        break;

      default: /* set Legendre Gauss-Lobatto points, weights, ... */
        LGL_x_wquad(ni, Xb, Wq);
        //Gauss_wquad_from_symm_x(npoints, x, w); // test Gauss_wquad_from_symm_x

        /* get analysis & synthesis matrix for Legendre basis,
           could be useful for filtering, but not needed for interpolation */
        LGL_AT_ST_matrices(ni, Xb, Wq, AT, ST);
      }

      /* diff matrix DT for Lagrange interp. poly basis */
      Lagrange_winterp(ni, Xb, WL);
      Lagrange_DT(ni, Xb, WL, DT);
    }
    /* set Legendre polys as basis since AT and ST are for Legendre basis */
    gridpoints->basis[typ] = basis_normLegendreP;
  }
  return 0;
}

/* free all arrays with grid points and such */
int free_gridpoints(tMesh *mesh)
{
  int ni, typ;

  /* free points, diff matrices, and such */
  for(ni=1; ni<=gridpoints->nmax; ni++)
  {
    for(typ=0; typ<P_NTYPES; typ++)
    {
      free_array(gridpoints->Dt[ni][typ]);
      free_array(gridpoints->At[ni][typ]);
      free_array(gridpoints->St[ni][typ]);
      free_array(gridpoints->Xb[ni][typ]);
      free_array(gridpoints->Wq[ni][typ]);
      free_array(gridpoints->WL[ni][typ]);
    }
  }
  free(gridpoints->Dt);
  free(gridpoints->At);
  free(gridpoints->St);
  free(gridpoints->Xb);
  free(gridpoints->Wq);
  free(gridpoints->WL);

  return 0;
}

/*
tArray *get_gridpoints_Dt(int ni, int typ)
{
  return gridpoints->Dt[ni][typ];
}
tArray *get_gridpoints_At(int ni, int typ)
{
  return gridpoints->At[ni][typ];
}
tArray *get_gridpoints_St(int ni, int typ)
{
  return gridpoints->St[ni][typ];
}
tArray *get_gridpoints_Xb(int ni, int typ)
{
  return gridpoints->Xb[ni][typ];
}
tArray *get_gridpoints_Wq(int ni, int typ)
{
  return gridpoints->Wq[ni][typ];
}
tArray *get_gridpoints_WL(int ni, int typ)
{
  return gridpoints->WL[ni][typ];
}
*/
