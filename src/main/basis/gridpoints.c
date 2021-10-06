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

  for(typ=0; typ<P_NTYPES; typ++)
  {
    gridpoints->Xb[typ] = calloc(nmax+1, sizeof(gridpoints->Xb[typ][0]));
    if(!(gridpoints->Xb[typ]))
      errorexit("out of memory for points");
    gridpoints->Wq[typ] = calloc(nmax+1, sizeof(gridpoints->Wq[typ][0]));
    if(!(gridpoints->Wq[typ]))
      errorexit("out of memory for integr. weights");

    gridpoints->WL[typ] = calloc(nmax+1, sizeof(gridpoints->WL[typ][0]));
    if(!(gridpoints->WL[typ]))
      errorexit("out of memory for Lagrange interp. weights");

    gridpoints->Dt[typ] = calloc(nmax+1, sizeof(gridpoints->Dt[typ][0]));
    if(!(gridpoints->Dt[typ]))
      errorexit("out of memory for diff. matrices");

    gridpoints->Dpt[typ] = calloc(nmax+1, sizeof(gridpoints->Dpt[typ][0]));
    if(!(gridpoints->Dpt[typ]))
      errorexit("out of memory for forward diff. matrices");
    gridpoints->Dmt[typ] = calloc(nmax+1, sizeof(gridpoints->Dmt[typ][0]));
    if(!(gridpoints->Dmt[typ]))
      errorexit("out of memory for backward diff. matrices");

    gridpoints->At[typ] = calloc(nmax+1, sizeof(gridpoints->At[typ][0]));
    if(!(gridpoints->At[typ]))
      errorexit("out of memory for ana. matrices");
    gridpoints->St[typ] = calloc(nmax+1, sizeof(gridpoints->St[typ][0]));
    if(!(gridpoints->St[typ]))
      errorexit("out of memory for syn. matrices");
  }

  /* allocate arrays */
  for(ni=1; ni<=nmax; ni++)
  {
    int n[3];

    n[0] = n[1] = ni;
    n[2] = 1;
    for(typ=0; typ<P_NTYPES; typ++)
    {
      gridpoints->Dt[typ][ni] = alloc_array(n);
      gridpoints->Dpt[typ][ni] = alloc_array(n);
      gridpoints->Dmt[typ][ni] = alloc_array(n);
      gridpoints->At[typ][ni] = alloc_array(n);
      gridpoints->St[typ][ni] = alloc_array(n);
    }
    n[0] = ni;
    n[1] = n[2] = 1;
    for(typ=0; typ<P_NTYPES; typ++)
    {
      gridpoints->Xb[typ][ni] = alloc_array(n);
      gridpoints->Wq[typ][ni] = alloc_array(n);
      gridpoints->WL[typ][ni] = alloc_array(n);
    }
  }

  /* Set points, weights, diff, and other matrices.
     Points in patch coords are then
     X=0.5*((a+b)+(b-a)*Xb),   where a=bbox[0], b=bbox[1] */
  for(typ=0; typ<P_NTYPES; typ++)
  {
    for(ni=1; ni<=nmax; ni++)
    {
      double *Xb = gridpoints->Xb[typ][ni]->d;
      double *Wq = gridpoints->Wq[typ][ni]->d;
      double *WL = gridpoints->WL[typ][ni]->d;
      double *DT = gridpoints->Dt[typ][ni]->d;
      double *DpT = gridpoints->Dpt[typ][ni]->d;
      double *DmT = gridpoints->Dmt[typ][ni]->d;
      double *AT = gridpoints->At[typ][ni]->d;
      double *ST = gridpoints->St[typ][ni]->d;

      /* set up desired points */
      switch(typ)
      {
      case P_UNIFORM: /* set equally spaced points and their weights */
        //uniform_x_wGaussquad(ni, Xb, Wq); //can give negative rms
        uniform_x_wTrapez(ni, Xb, Wq);

        /* get analysis & synthesis matrix for Legendre basis */
        Legendre_AT_ST_matrices(ni, Xb, Wq, AT, ST);

        /* Lagrange interp. weights WL and diff matrix DT for Lagrange
           interp. poly basis */
        Lagrange_winterp(ni, Xb, WL);
        //Lagrange_DT(ni, Xb, WL, DT); // very inaccurate for large ni
        fd_lopderiv_DT_uniform(ni, Xb, 3,0, DT);
        //printarray(gridpoints->Dt[typ][ni]);
        fd_lopderiv_DT_uniform(ni, Xb, 3,+1, DpT);
        fd_lopderiv_DT_uniform(ni, Xb, 3,-1, DmT);

        break;

      default: /* set Legendre Gauss-Lobatto points, weights, ... */
        LGL_x_wquad(ni, Xb, Wq);
        //Gauss_wquad_from_symm_x(npoints, x, w); // test Gauss_wquad_from_symm_x

        /* get analysis & synthesis matrix for Legendre basis,
           could be useful for filtering, but not needed for interpolation */
        LGL_AT_ST_matrices(ni, Xb, Wq, AT, ST);

        /* Lagrange interp. weights WL and diff matrix DT for Lagrange
           interp. poly basis */
        Lagrange_winterp(ni, Xb, WL);
        Lagrange_DT(ni, Xb, WL, DT);
        Lagrange_DT(ni, Xb, WL, DpT);
        Lagrange_DT(ni, Xb, WL, DmT);
      }
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
  for(typ=0; typ<P_NTYPES; typ++)
  {
    for(ni=1; ni<=gridpoints->nmax; ni++)
    {
      free_array(gridpoints->Dt[typ][ni]);
      free_array(gridpoints->Dpt[typ][ni]);
      free_array(gridpoints->Dmt[typ][ni]);
      free_array(gridpoints->At[typ][ni]);
      free_array(gridpoints->St[typ][ni]);
      free_array(gridpoints->Xb[typ][ni]);
      free_array(gridpoints->Wq[typ][ni]);
      free_array(gridpoints->WL[typ][ni]);
    }
    free(gridpoints->Dt[typ]);
    free(gridpoints->Dpt[typ]);
    free(gridpoints->Dmt[typ]);
    free(gridpoints->At[typ]);
    free(gridpoints->St[typ]);
    free(gridpoints->Xb[typ]);
    free(gridpoints->Wq[typ]);
    free(gridpoints->WL[typ]);
  }

  /* now set all in gridpoints back to 0 */
  memset(gridpoints, 0, sizeof(gridpoints[0]));

  return 0;
}

/*
tArray *get_gridpoints_Dt(int typ, int ni)
{
  return gridpoints->Dt[typ][ni];
}
tArray *get_gridpoints_At(int typ, int ni)
{
  return gridpoints->At[typ][ni];
}
tArray *get_gridpoints_St(int typ, int ni)
{
  return gridpoints->St[typ][ni];
}
tArray *get_gridpoints_Xb(int typ, int ni)
{
  return gridpoints->Xb[typ][ni];
}
tArray *get_gridpoints_Wq(int typ, int ni)
{
  return gridpoints->Wq[typ][ni];
}
tArray *get_gridpoints_WL(int typ, int ni)
{
  return gridpoints->WL[typ][ni];
}
*/
