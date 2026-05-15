/* grid.c */
/* Wolfgang Tichy, 10/2019 */

#include "nmesh.h"
#include "basis.h"

#define PR 0


/* global grid point structure */
tGridPoints gridpoints[1];

/* init gridpoints structure */
int gridpoints_init(tMesh *mesh)
{
  int stencilsize = Geti(Par("fd_stencilsize"));
  int lopsidesize = Geti(Par("fd_lopsidesize"));
  int nmax, ni, typ;

  /* get mem. for grid points, diff. matrices, ... */
  gridpoints_alloc(mesh);
  nmax = gridpoints->nmax;

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
      tArray *gDt  = gridpoints->Dt[typ][ni];
      tArray *gDpt = gridpoints->Dpt[typ][ni];
      tArray *gDmt = gridpoints->Dmt[typ][ni];
      double *DT = gDt->d;
      double *DpT = gDpt->d;
      double *DmT = gDmt->d;
      double *AT = gridpoints->At[typ][ni]->d;
      double *ST = gridpoints->St[typ][ni]->d;

      /* set up desired points */
      switch(typ)
      {
      case P_UNIFORM: /* set equally spaced points and their weights */
        //uniform_x_wGaussquad(ni, Xb, Wq); //can give negative rms
        //printarray(gridpoints->Wq[typ][ni]);
        //if(ni==20) errorexit("quad weights above");
        uniform_x_wTrapez(ni, Xb, Wq);

        /* get analysis & synthesis matrix for Legendre basis */
        Legendre_AT_ST_matrices(ni, Xb, Wq, AT, ST);
        //printarray(gridpoints->At[typ][ni]);
        //printarray(gridpoints->St[typ][ni]);
        //tArray *u = alloc_array1d(ni);
        //for(int i=0; i<ni; i++) u->d[i] = Xb[i]*Xb[i]*Xb[i]*Xb[i]*Xb[i];
        //mm_array0(gridpoints->At[typ][ni], u, gridpoints->Wq[typ][ni]);
        //free_array(u);
        //printarray(gridpoints->Wq[typ][ni]);
        //if(ni==9) errorexit("ana and syn matr. are above");

        /* Lagrange interp. weights WL */
        Lagrange_winterp(ni, Xb, WL);

        /* diff matrix DT for Lagrange interp. poly basis */
        //Lagrange_DT(ni, Xb, WL, DT);
        /* Lagrange_DT is very inaccurate for large ni on uniform grids and
           can also lead to instabilities in the evolution */

        /* if ni > stencilsize the fd matrices are sparse so we set
           the non-zero ranges */
        if(ni > stencilsize)
        {
          alloc_2darray_irange_of_j(gDt);
          alloc_2darray_irange_of_j(gDpt);
          alloc_2darray_irange_of_j(gDmt);
        }

        /* finite difference diff matrix DT */
        fd_lopderiv_DT_uniform(ni, Xb, stencilsize,0, DT, gDt->range);
        fd_lopderiv_DT_uniform(ni, Xb, stencilsize,+lopsidesize, DpT, gDpt->range);
        fd_lopderiv_DT_uniform(ni, Xb, stencilsize,-lopsidesize, DmT, gDmt->range);
        //printarray(gridpoints->Dt[typ][ni]);
        //printarray(gridpoints->Dpt[typ][ni]);
        //printarray(gridpoints->Dmt[typ][ni]);
        //if(ni==9) errorexit("fd diff. matrices are above");

        break;

      case P_CHEBEXTR: /* set Cheb. extrema points, weights, ... */
        ChebyshevExtrema_x(ni, Xb);
        Gauss_wquad_from_symm_x(ni, Xb, Wq);
        //printarray(gridpoints->Wq[typ][ni]);
        //tArray *u = alloc_array1d(ni);
        //for(int i=0; i<ni; i++) u->d[i] = 5.*(Xb[i]-1.)*(Xb[i]-1.)*(Xb[i]-1.)*(Xb[i]-1.);
        //double Gauss_Integral(tArray *Wq, tArray *func);
        //printf("int(u)=%.16g   u", Gauss_Integral(gridpoints->Wq[typ][ni], u));
        //printarray(u);
        //free_array(u);
        //if(ni==6) errorexit("Wq matrices are above");

        /* get analysis & synthesis matrix */
        ChebyshevExtrema_AT_ST(ni, AT, ST);
        //printarray(gridpoints->At[typ][ni]);
        //printarray(gridpoints->St[typ][ni]);
        //if(ni==3) errorexit("AT ST matrices are above");

        /* Lagrange interp. weights WL */
        Lagrange_winterp(ni, Xb, WL);

        /* diff matrices */
        //Lagrange_DT(ni, Xb, WL, DT); // Lagrange_DT seems less accurate
        //Lagrange_DT(ni, Xb, WL, DpT);
        //Lagrange_DT(ni, Xb, WL, DmT);
        //printarray(gridpoints->Dt[typ][ni]);
        ChebyshevExtrema_DT(ni, DT);
        ChebyshevExtrema_DT(ni, DpT);
        ChebyshevExtrema_DT(ni, DmT);
        //printarray(gridpoints->Dt[typ][ni]);
        //if(ni==9) errorexit("diff. matrices are above");

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

    /* set basis function for each grid type */
    switch(typ)
    {
    case P_CHEBEXTR:
      gridpoints->basis[typ] = Chebyshev_basisfunc;
      break;
    case P_LGL:
    case P_UNIFORM:
    default: /* set Legendre poly. as basis if AT,ST are for Legendre basis */
      gridpoints->basis[typ] = basis_normLegendreP;
    }
  }

  /* set interp. matrices */
  for(ni=1; ni<=nmax; ni++)
  {
    tArray *Xb   = gridpoints->Xb[P_UNIFORM][ni];
    //tArray *rq   = gridpoints->Wq[P_UNIFORM][ni];
    tArray *Rt   = gridpoints->UNI_to_nLGLt[ni];
    tArray *Rto2 = gridpoints->UNI_to_no2LGLt[ni];

    tArray *X    = gridpoints->Xb[P_LGL][ni];
    tArray *WL   = gridpoints->WL[P_LGL][ni];
    //tArray *wq   = gridpoints->Wq[P_LGL][ni];
    tArray *Pt   = alloc_array2d(ni, ni);

    tArray *Xo2  = gridpoints->Xb[P_LGL][ni/2];
    tArray *WLo2 = gridpoints->WL[P_LGL][ni/2];
    //tArray *wqo2 = gridpoints->Wq[P_LGL][ni/2];
    tArray *Pto2 = alloc_array2d(ni/2, ni);

    /* set matrix from LGL with ni to UNIFORM */
    Lagrange_InterpMatT(X, WL, Xb, Pt);
    /* set matrix from LGL with ni/2 to UNIFORM */
    Lagrange_InterpMatT(Xo2, WLo2, Xb, Pto2);
    //PRFs(": Pto2");printarray_matrix0(Pto2);

    /* Now calc UNI_to_nLGLt=Rt and UNI_to_no2LGLt=Rto2 */
    //Inverse_InterpMatT_rq(Pt, wq, rq, Rt);
    //Inverse_InterpMatT_rq(Pto2, wqo2, rq, Rto2);
    Inverse_InterpMatT_best_rq(Pt, Rt);
    Inverse_InterpMatT_best_rq(Pto2, Rto2);

    //tArray *RP = alloc_array2d(ni/2, ni/2);
    //array_transpose01_inplace(Pto2);
    //mm_array_indir(Rto2, Pto2, 0, RP);
    //PRFs(": P");printarray_matrix0(Pto2);
    //PRFs(": Rt");printarray_matrix0(Rto2);
    //PRFs(": RP");printarray_matrix0(RP);
    //free_array(RP);
    //if(ni==6) errorexit("dasdsadsa");

    free_array(Pto2);
    free_array(Pt);
  }

  return 0;
}

/* allocate room for stuff in tGridPoints gridpoints[1] */
int gridpoints_alloc(tMesh *mesh)
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
  /* mem. for interp. matrices */
  gridpoints->UNI_to_nLGLt = calloc(nmax+1, sizeof(gridpoints->UNI_to_nLGLt[0]));
  if(!(gridpoints->UNI_to_nLGLt))
    errorexit("out of memory for interp. matrices UNI_to_nLGLt");
  gridpoints->UNI_to_no2LGLt = calloc(nmax+1, sizeof(gridpoints->UNI_to_no2LGLt[0]));
  if(!(gridpoints->UNI_to_no2LGLt))
    errorexit("out of memory for interp. matrices UNI_to_no2LGLt");

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

    gridpoints->UNI_to_nLGLt[ni]   = alloc_array2d(ni, ni);
    gridpoints->UNI_to_no2LGLt[ni] = alloc_array2d(ni, ni/2); //this is transpose
  }

  return 0;
}

/* free all arrays with grid points and such */
int gridpoints_free(tMesh *mesh)
{
  int ni, typ;

  /* free things that have typ: points, diff matrices, and such */
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

  /* free UNI_to_nLGLt, UNI_to_no2LGLt */
  for(ni=1; ni<=gridpoints->nmax; ni++)
  {
    free_array(gridpoints->UNI_to_nLGLt[ni]);
    free_array(gridpoints->UNI_to_no2LGLt[ni]);
  }
  free(gridpoints->UNI_to_nLGLt);
  free(gridpoints->UNI_to_no2LGLt);

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
