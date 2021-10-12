/* coordinates.c */
/* Wolfgang Tichy, 2/2019 */

#include "nmesh.h"
#include "coordinates.h"

#define PR 0


/* frequently used global vars */
tcoordinates coordinates[1];


/* try to enable all coord vars, return 1 if we have dat */
int coordinates_coordvars_enabled(tNode *node)
{
  tPat *pat = node->pat;
  tMesh *mesh = pat->mesh;
  tDat *dat = node->dat;
  tCoordInfo *CI = pat->CI;
  int iX, ix, idXdx, idet_dXbdx;
  int isqrtdet2g_o_det3gamma0, isqrtgdiagx;
  int surface_metric, sqrtdet2g_o_det3gamma, sqrtgdiag;
  int itmp;
  int f, d;

  /* do nothing if this is not my node */
  if(!dat) return 0;

  /* if coords are set already do nothing */
  if(dat->coords_set) return 1;

  /* set some indices */
  iX = Ind("X");
  ix = Ind("x");
  idXdx = Ind("dXdx");
  idet_dXbdx = Ind("det_dXbdx");
  isqrtdet2g_o_det3gamma0 = Ind("sqrtdet2g_o_det3gamma0");
  isqrtgdiagx = Ind("sqrtgdiagx");
  itmp = Ind("coordinates_tmp1");

  /* which surface info do we set */
  surface_metric = Par("coordinates_surface_metric");
  sqrtdet2g_o_det3gamma  = Getv(surface_metric, "sqrtdet2g_o_det3gamma");
  sqrtgdiag              = Getv(surface_metric, "sqrtgdiag");

  /* give all these memory: */
  enablevar_innode(node, iX);
  enablevar_innode(node, iX+1);
  enablevar_innode(node, iX+2);
  enablevar_innode(node, ix);
  enablevar_innode(node, ix+1);
  enablevar_innode(node, ix+2);
  enablevar_innode(node, idXdx);
  enablevar_innode(node, idXdx+3);
  enablevar_innode(node, idXdx+6);
  enablevar_innode(node, idet_dXbdx);
  enablevar_innode(node, itmp);
  enablevar_innode(node, itmp+1);
  enablevar_innode(node, itmp+2);
  if(1 || sqrtdet2g_o_det3gamma)
    enablevar_innode(node, isqrtdet2g_o_det3gamma0);
  if(1 || sqrtgdiag)
    enablevar_innode(node, isqrtgdiagx);

  /* give oC surface coords memory if node has corresponding surface */
  for(f=0; f<6; f++)
  {
    tBface *bface0 = node->pat->bfaces[f];
    int nnb = node->nfnb[f];
    int ioC0;

    if(bface0) ioC0 = bface0->ioC0_0;
    else       ioC0 = -1;

    if(ioC0>0 && node->patface[f] && nnb)
    {
      enablevarcomp_innode(node, ioC0+f);
      enablevarcomp_innode(node, ioC0+6+f);
    }
  }

  /* enable extra vars for cubed spheres */
  for(f=0; f<6; f++)
    if(CI->iSurf[f] > 0)
    {
      enablevar_innode(node, CI->iSurf[f]);
      for(d=0; d<3; d++)
        if(CI->idSurfdX[f][d] > 0)
          enablevar_innode(node, CI->idSurfdX[f][d]);
    }

  /* enable extra vars to store stuff in the middle between grid points */
  if(Getb(Par("coordinates_midpoint_data")))
  {
    int iXm_dXdx = Ind("Xm_dXdx");
    int iYm_dXdx = Ind("Ym_dXdx");
    int iZm_dXdx = Ind("Zm_dXdx");
    int iXm_sqrtgdiagx = Ind("Xm_sqrtgdiagx");
    int iYm_sqrtgdiagx = Ind("Ym_sqrtgdiagx");
    int iZm_sqrtgdiagx = Ind("Zm_sqrtgdiagx");

    enablevar_innode(node, iXm_dXdx);
    enablevar_innode(node, iXm_dXdx+3);
    enablevar_innode(node, iXm_dXdx+6);
    enablevar_innode(node, iYm_dXdx);
    enablevar_innode(node, iYm_dXdx+3);
    enablevar_innode(node, iYm_dXdx+6);
    enablevar_innode(node, iZm_dXdx);
    enablevar_innode(node, iZm_dXdx+3);
    enablevar_innode(node, iZm_dXdx+6);

    enablevar_innode(node, iXm_sqrtgdiagx);
    enablevar_innode(node, iYm_sqrtgdiagx);
    enablevar_innode(node, iZm_sqrtgdiagx);
  }

  return 1;
}

/* (re)initialize coordinates in a node */
int coordinates_init_node(tNode *node)
{
  tPat *pat = node->pat;
  tMesh *mesh = pat->mesh;
  tDat *dat = node->dat;
  tCoordInfo *CI = pat->CI;
  int *n = node->n;
  int i,j,k, d,e, f;
  int vars_on = coordinates_coordvars_enabled(node);
  int iX = Ind("X");
  int ix = Ind("x");
  int idXdx = Ind("dXdx");
  int surface_metric, sqrtdet2g_o_det3gamma, sqrtgdiag;
  int i3metric = MeshVarIndLax(mesh, Gets(Par("coordinates_3metric")));
  double *pX[] = { Vard(node,iX), Vard(node,iX+1), Vard(node,iX+2) };
  double *px[] = { Vard(node,ix), Vard(node,ix+1), Vard(node,ix+2) };
  double *pdXdx[3][3]
            = { {Vard(node,idXdx),   Vard(node,idXdx+1), Vard(node,idXdx+2)},
                {Vard(node,idXdx+3), Vard(node,idXdx+4), Vard(node,idXdx+5)},
                {Vard(node,idXdx+6), Vard(node,idXdx+7), Vard(node,idXdx+8)} };
  double *det_dXbdx = Vard(node, Ind("det_dXbdx"));
  double dXbdX[3];
  double det_dXbYbZb_dXYZ;

  /* do nothing if coords are set already or if vars are off */
  if(!dat) return 0;
  if(dat->coords_set) return 0;
  if(!vars_on) return 0;

  if(PR) { PRF;printf(": nid%ld\n", node->nid); }

  /* which surface info do we set */
  surface_metric = Par("coordinates_surface_metric");
  sqrtdet2g_o_det3gamma  = Getv(surface_metric, "sqrtdet2g_o_det3gamma");
  sqrtgdiag              = Getv(surface_metric, "sqrtgdiag");

  /* get det of dXb/dX */
  dXbYbZb_dXYZ(node, dXbdX);
  det_dXbYbZb_dXYZ = dXbdX[0] * dXbdX[1] * dXbdX[2];

  /* set coords */
  forijk(i,j,k, n)
  {
    double Xb[] = { node_Xb(node,0)->d[i], node_Xb(node,1)->d[j],
                    node_Xb(node,2)->d[k] };
    double X[3], x[3], dXdx[3][3];
    int ijk = Ind_n(i,j,k, n);

    /* get X from Xb */
    XYZ_of_XbYbZb(node, Xb, X);
    for(d=0; d<3; d++) pX[d][ijk] = X[d];

    /* now set x, dXdx, det(dXb/dx) */
    if(pat->dXYZ_dxyz)
    {
      pat->dXYZ_dxyz(pat, node, -1, X, x, dXdx);
      if(0) if( !finit(x[0]) || !finit(x[1]) || !finit(x[2]) )
      {
        printpatch(pat);
        printCI(pat);
        printnode(node);
        pr3v("X", X);
        pr3v("x", x);
        errorexit("x is NAN");
      }

      for(d=0; d<3; d++)
      {
        px[d][ijk] = x[d];
        for(e=0; e<3; e++) pdXdx[d][e][ijk] = dXdx[d][e];
      }
      det_dXbdx[ijk] = det_dXbYbZb_dXYZ * det_3Dmatrix(dXdx);
    }
    else /* assume X,Y,Z are Cartesian*/
    {
      for(d=0; d<3; d++)
      {
        px[d][ijk] = pX[d][ijk];
        pdXdx[d][d][ijk] = 1.;
      }
      det_dXbdx[ijk] = det_dXbYbZb_dXYZ;
    }
  }

  /* set sqrtdet2g_o_det3gamma on node faces */
  if(sqrtdet2g_o_det3gamma)
    coordinates_set_sqrtdet2g_o_det3gamma_var(node, i3metric,
                                              Ind("sqrtdet2g_o_det3gamma0"));
  /* set sqrtgdiag */
  if(sqrtgdiag)
    coordinates_set_sqrtgdiag_var(node, idXdx, i3metric, Ind("sqrtgdiagx"));

  /* set oC surface coords, for now this is off */
  for(f=0; f<6; f++)
  {
    tBface *bface0 = node->pat->bfaces[f];
    tPat *opat;
    tBface *obface;
    int nnb = node->nfnb[f];
    int ioC0;

    if(bface0) ioC0 = bface0->ioC0_0;
    else       ioC0 = -1;

    if(ioC0>0 && node->patface[f] && nnb>0)
    {
      int dir = f/2;
      int pl = (n[dir]-1)*(f%2);
      int d1 = Dir1_norm(dir);
      int d2 = Dir2_norm(dir);
      double *oC[] = { Vard(node,ioC0+f), Vard(node,ioC0+6+f) };

      forplaneN(dir, i,j,k, n, pl)
      {
        int ijk = Ind_n(i,j,k, n);
        int ind = Ind_n_norm(i,j,k, n, dir);
        double x[] = { px[0][ijk], px[1][ijk], px[2][ijk] };
        double X[] = { pX[0][ijk], pX[1][ijk], pX[2][ijk] };
        double C[] = { X[d1], X[d2] };
        double oX[3];
        int odir, od1, od2, pi;

        obface = first_obface_of_bface_containing_point(pat, f, C);
        errorexit("Wrong strategy: We should find neighbors from node->fnb "
                  "and not use bfaces. We should then save oC as well as the "
                  "neighbor index in gridvars on this node. But even then we "
                  "still haven't solved the case where the point is in two or "
                  "more nighbors... So let's give up on saving oC!!!");
        opat = obface->pat;
        odir = obface->f/2;
        od1 = Dir1_norm(odir);
        od2 = Dir2_norm(odir);

        pi = p_XYZ_of_xyz(opat, oX, x);
        if(pi<0) errorexit("x should be be in opat!!!");
        oC[0][ind] = oX[od1];
        oC[1][ind] = oX[od2];
      }
    }
  }

  /* set surface vars */
  for(f=0; f<6; f++)
    if(CI->iSurf[f] > 0)
    {
      int dir = f/2;
      int pl = (n[dir]-1)*(f%2);
      int d1  = Dir1_norm(dir);
      int d2  = Dir2_norm(dir);
      double *sig = Vard(node, CI->iSurf[f]);
      double C[2], F;

      forplaneN(dir, i,j,k, n, pl)
      {
        int ijk = Ind_n(i,j,k, n);
        int ind = Ind_n_norm(i,j,k, n, dir);

        C[0] = pX[d1][ijk];
        C[1] = pX[d2][ijk];
        CI->FSurf[f](pat, f, C, &F);
        sig[ind] = F;
      }

      /* and their derivs */
      if(CI->idSurfdX[f][d1] > 0)
      {
        double *dsig1 = Vard(node, CI->idSurfdX[f][d1]);
        double *dsig2 = Vard(node, CI->idSurfdX[f][d2]);
        double dF[2];

        forplaneN(dir, i,j,k, n, pl)
        {
          int ijk = Ind_n(i,j,k, n);
          int ind = Ind_n_norm(i,j,k, n, dir);

          C[0] = pX[d1][ijk];
          C[1] = pX[d2][ijk];
          CI->dFSurfdC[f](pat, f, C, dF);
          dsig1[ind] = dF[0];
          dsig2[ind] = dF[1];
        }
      }
    }

  /* set extra vars to store stuff in the middle between grid points */
  if(Getb(Par("coordinates_midpoint_data")))
  {
    int iXm_dXdx = Ind("Xm_dXdx");
    int iYm_dXdx = Ind("Ym_dXdx");
    int iZm_dXdx = Ind("Zm_dXdx");
    int iXm_sqrtgdiagx = Ind("Xm_sqrtgdiagx");
    int iYm_sqrtgdiagx = Ind("Ym_sqrtgdiagx");
    int iZm_sqrtgdiagx = Ind("Zm_sqrtgdiagx");
    double *pXm_dXdx[3][3] =
     {{Vard(node,iXm_dXdx),   Vard(node,iXm_dXdx+1), Vard(node,iXm_dXdx+2)},
      {Vard(node,iXm_dXdx+3), Vard(node,iXm_dXdx+4), Vard(node,iXm_dXdx+5)},
      {Vard(node,iXm_dXdx+6), Vard(node,iXm_dXdx+7), Vard(node,iXm_dXdx+8)} };
    double *pYm_dXdx[3][3] =
     {{Vard(node,iYm_dXdx),   Vard(node,iYm_dXdx+1), Vard(node,iYm_dXdx+2)},
      {Vard(node,iYm_dXdx+3), Vard(node,iYm_dXdx+4), Vard(node,iYm_dXdx+5)},
      {Vard(node,iYm_dXdx+6), Vard(node,iYm_dXdx+7), Vard(node,iYm_dXdx+8)} };
    double *pZm_dXdx[3][3] =
     {{Vard(node,iZm_dXdx),   Vard(node,iZm_dXdx+1), Vard(node,iZm_dXdx+2)},
      {Vard(node,iZm_dXdx+3), Vard(node,iZm_dXdx+4), Vard(node,iZm_dXdx+5)},
      {Vard(node,iZm_dXdx+6), Vard(node,iZm_dXdx+7), Vard(node,iZm_dXdx+8)} };
    int *Xm_n = Arrn(VarA(node, iXm_dXdx));
    int *Ym_n = Arrn(VarA(node, iYm_dXdx));
    int *Zm_n = Arrn(VarA(node, iZm_dXdx));

    /* set Xm_dXdx, Ym_dXdx, Zm_dXdx */
    {
      forijk(i,j,k, n)
      {
        double Xmid[3];
        double Ymid[3];
        double Zmid[3];
        double X[3], x[3], dXdx[3][3];
        int gotXmid, gotYmid, gotZmid;
        int ijk;

        //printf("Xm_n[]=%d %d %d\n", Xm_n[0],Xm_n[1],Xm_n[2]);
        //printf("Ym_n[]=%d %d %d\n", Ym_n[0],Ym_n[1],Ym_n[2]);
        //printf("Zm_n[]=%d %d %d\n", Zm_n[0],Zm_n[1],Zm_n[2]);
        //exit(88);

        /* find mid points in X-dir */
        gotXmid = set_nodemidpoint_XbYbZb(node, i,j,k, 0, Xmid);
        gotYmid = set_nodemidpoint_XbYbZb(node, i,j,k, 1, Ymid);
        gotZmid = set_nodemidpoint_XbYbZb(node, i,j,k, 2, Zmid);

        /* now set x, dXdx */
        if(gotXmid)
        {
          XYZ_of_XbYbZb(node, Xmid, X);
          set_xyz_dXYZdxyz(pat, node, -1, X, x, dXdx);
          ijk = Ind_n(i,j,k, Xm_n);

          for(d=0; d<3; d++)
            for(e=0; e<3; e++) pXm_dXdx[d][e][ijk] = dXdx[d][e];
        }
        if(gotYmid)
        {
          XYZ_of_XbYbZb(node, Ymid, X);
          set_xyz_dXYZdxyz(pat, node, -1, X, x, dXdx);
          ijk = Ind_n(i,j,k, Ym_n);

          for(d=0; d<3; d++)
            for(e=0; e<3; e++) pYm_dXdx[d][e][ijk] = dXdx[d][e];
        }
        if(gotZmid)
        {
          XYZ_of_XbYbZb(node, Zmid, X);
          set_xyz_dXYZdxyz(pat, node, -1, X, x, dXdx);
          ijk = Ind_n(i,j,k, Zm_n);

          for(d=0; d<3; d++)
            for(e=0; e<3; e++) pZm_dXdx[d][e][ijk] = dXdx[d][e];
        }
      } /* end forijk */
    }

    /* set Xm_sqrtgdiag, Ym_sqrtgdiag, Zm_sqrtgdiag */
    {
      /* 3 arrays for dXdx on midpoints */
      tArray *AXm_dXdx[3][3] =
      {{VarA(node,iXm_dXdx),   VarA(node,iXm_dXdx+1), VarA(node,iXm_dXdx+2)},
       {VarA(node,iXm_dXdx+3), VarA(node,iXm_dXdx+4), VarA(node,iXm_dXdx+5)},
       {VarA(node,iXm_dXdx+6), VarA(node,iXm_dXdx+7), VarA(node,iXm_dXdx+8)}};
      tArray *AYm_dXdx[3][3] =
      {{VarA(node,iYm_dXdx),   VarA(node,iYm_dXdx+1), VarA(node,iYm_dXdx+2)},
       {VarA(node,iYm_dXdx+3), VarA(node,iYm_dXdx+4), VarA(node,iYm_dXdx+5)},
       {VarA(node,iYm_dXdx+6), VarA(node,iYm_dXdx+7), VarA(node,iYm_dXdx+8)}};
      tArray *AZm_dXdx[3][3] =
      {{VarA(node,iZm_dXdx),   VarA(node,iZm_dXdx+1), VarA(node,iZm_dXdx+2)},
       {VarA(node,iZm_dXdx+3), VarA(node,iZm_dXdx+4), VarA(node,iZm_dXdx+5)},
       {VarA(node,iZm_dXdx+6), VarA(node,iZm_dXdx+7), VarA(node,iZm_dXdx+8)}};
      /* 3 arrays for sqrtgdiag on midpoints */
      tArray *AXm_sqrtgdiag[3] = { VarA(node, iXm_sqrtgdiagx),
                                   VarA(node, iXm_sqrtgdiagx+1),
                                   VarA(node, iXm_sqrtgdiagx+2) };
      tArray *AYm_sqrtgdiag[3] = { VarA(node, iYm_sqrtgdiagx),
                                   VarA(node, iYm_sqrtgdiagx+1),
                                   VarA(node, iYm_sqrtgdiagx+2) };
      tArray *AZm_sqrtgdiag[3] = { VarA(node, iZm_sqrtgdiagx),
                                   VarA(node, iZm_sqrtgdiagx+1),
                                   VarA(node, iZm_sqrtgdiagx+2) };
      tArray *Am_g[6];

      /* NOTE: for now we only use a flat metric in DG */
      for(d=0; d<6; d++) Am_g[d] = NULL;
      if(i3metric>=0)
        errorexit("implement case where we use non-flat metric in DG");

      /* now set sqrtgdiag */
      coordinates_set_sqrtgdiag_array(node, AXm_dXdx, Am_g, AXm_sqrtgdiag);
      coordinates_set_sqrtgdiag_array(node, AYm_dXdx, Am_g, AYm_sqrtgdiag);
      coordinates_set_sqrtgdiag_array(node, AZm_dXdx, Am_g, AZm_sqrtgdiag);
    }
  }

  /* mark coords as set */
  dat->coords_set = 1;

  return 0;
}


/* initialize coordinates in each patch */
int coordinates_init(tMesh *mesh)
{
  PRF;printf(":\n");
  int surface_metric = Par("coordinates_surface_metric");

  /* global par values */
  coordinates->sqrtdet2g_o_det3gamma = Getv(surface_metric,
                                            "sqrtdet2g_o_det3gamma");
  /* set some global vars */
  coordinates->idXdx = Ind("dXdx");
  coordinates->itmp1 = Ind("coordinates_tmp1");
  coordinates->idet_dXbdx = Ind("det_dXbdx");
  coordinates->isqrtdet2g_o_det3gamma0 = Ind("sqrtdet2g_o_det3gamma0");
  coordinates->isqrtgdiagx = Ind("sqrtgdiagx");

  //PRF;printf(":  coordinates->idet_dXbdx=%d\n",  coordinates->idet_dXbdx);

  formylnodes(mesh)
  {
    tNode *node = MyLnode;
    coordinates_init_node(node);
  }

  return 0;
}


/********************************************************************/
/* some functions to set surface metric for DG-surface terms */
/********************************************************************/


/* Write sqrt(det(2g)/det(3gamma_ij)) on node faces into the 6 vars
   sqrtdet2g_o_det3gamma^i. We calculate sqrtdet2g_o_det3gamma from
   var dXdx and the symm. 3-metric igxx in x-coords.
   If igxx<0 we assume a flat 3-metric. */
void coordinates_set_sqrtdet2g_o_det3gamma_var(tNode *node, int igxx,
                                               int isqrtdet2g_o_det3gamma0)
{
  int *n = node->n;
  /* arrays to compute 2 metric sqrtdet2g_o_det3gamma on faces */
  tArray *a3gamT = alloc_empty_array2d(3,3); /* 3x3 for transp. of 3-metric */
  tArray *a2J = alloc_array2d(3,2);   /* 3x2 for 2-Jacobian */
  tArray *a2g = alloc_array2d(2,2);   /* 2x2 for induced metric on surface */
  tArray *tmp = alloc_array2d(3,2);
  tArray *Ag[6]; /* arrays for 3-metric */
  double *gxx, *gxy, *gxz, *gyy, *gyz, *gzz;
  int i,j,k, d, f;

  /* arrays for 3 metric */
  if(igxx>=0)
  {
    for(d=0; d<6; d++) Ag[d] = VarA(node, igxx + d);
    if(Ag[0]==NULL) errorexit("array for gxx is NULL!");
  }
  else
  {
    for(d=0; d<6; d++) Ag[d] = NULL;
  }
  /* pointers to 3 metric */
  gxx = Arrd(Ag[0]);
  gxy = Arrd(Ag[1]);
  gxz = Arrd(Ag[2]);
  gyy = Arrd(Ag[3]);
  gyz = Arrd(Ag[4]);
  gzz = Arrd(Ag[5]);

  /* set sqrtd2g_o_d3g on each face */
  for(f=0; f<6; f++)
  {
    int dir = f/2;
    int p = (f%2)*(n[dir] - 1);
    double *sqrtd2g_o_d3g = Vard(node, isqrtdet2g_o_det3gamma0 + f);

    forplaneN(dir, i,j,k, n, p)
    {
      double det2g, det3gam;
      int ijk = Ind_n(i,j,k, n);
      int JK = Ind_n_norm(i,j,k, n, dir);

      /* get 2-Jacobian of dx/dXb, which is a 3x2 matrix */
      array_2dxdXb(node, ijk, dir, a2J);

      /* if no 3-metric is specified, we assume it is flat */
      if(igxx<0)
      {
        /* compute 2-metric from a2J^T a2J */
        mm_array0_norestrict(a2J,a2J, a2g);
        det3gam = 1.; /* 3-metric is assumed flat */
      }
      else
      {
        /* compute 2-metric from a2J^T a3gam a2J, where a3gam is 3x3 */
        double M[3][3] = { { gxx[ijk], gxy[ijk], gxz[ijk] },
                           { gxy[ijk], gyy[ijk], gyz[ijk] },
                           { gxz[ijk], gyz[ijk], gzz[ijk] } };
        /* Put transpose of 3-metric into a3gamT.
           Trick here is that C-arrays are row major. So when we put this
           into the col. major a3gamT, we get the transpose automatically! */
        point_array_d_to_data(a3gamT, M, 1);
        mm_array0(a3gamT,a2J, tmp); /* a3gam a2J -> tmp */
        mm_array0(a2J,tmp, a2g);    /* a2J^T tmp -> a2g */
        det3gam = det_3Dmatrix(M);  /* det of 3-metric */
        if(det3gam <= 0.) det3gam = 1.; /* avoid NaNs */
      }
      det2g = det_2_2_array(a2g);
      //NOTE: there should be a faster way to get det2g. There should
      //      be an analogue to \det(g)=\alpha^2\det(\gamma)
      sqrtd2g_o_d3g[JK] = sqrt(det2g / det3gam);
    }
  }

  /* free arrays */
  free_array(tmp);
  free_array(a2g);
  free_array(a2J);
  free_array(a3gamT);
}


/* Write sqrt(g^{XX,YY,ZZ}) (in Xb-coords) into the 3 vars sqrtgdiag^i.
   We calculate sqrtgdiag from var dXdx and the symm. 3-metric igxx
   in x-coords. If igxx<0 we assume a flat 3-metric. */
void coordinates_set_sqrtgdiag_var(tNode *node, int idXdx, int igxx,
                                   int isqrtgdiagx)
{
  /* 3 arrays for dXdx on midpoints */
  tArray *AdXdx[3][3] =
    { { VarA(node,idXdx),   VarA(node,idXdx+1), VarA(node,idXdx+2) },
      { VarA(node,idXdx+3), VarA(node,idXdx+4), VarA(node,idXdx+5) },
      { VarA(node,idXdx+6), VarA(node,idXdx+7), VarA(node,idXdx+8) } };
  /* 3 arrays for sqrtgdiag on grid points */
  tArray *Asqrtgdiag[3] = { VarA(node, isqrtgdiagx),
                            VarA(node, isqrtgdiagx+1),
                            VarA(node, isqrtgdiagx+2) };
  tArray *Ag[6]; /* arrays for 3-metric */
  int d;

  /* set arrays for 3 metric */
  if(igxx>=0)
  {
    for(d=0; d<6; d++) Ag[d] = VarA(node, igxx + d);
    if(Ag[0]==NULL) errorexit("array for gxx is NULL!");
  }
  else
  {
    for(d=0; d<6; d++) Ag[d] = NULL;
  }

  /* now set sqrtgdiag */
  coordinates_set_sqrtgdiag_array(node, AdXdx, Ag, Asqrtgdiag);
}

/* Write sqrt(g^{XX,YY,ZZ}) (in Xb-coords) into the 3 arrays Asqrtgdiag[i].
   We calculate Asqrtgdiag from dXdx[3][3] and the symm. 3-metric in Ag[6]
   in x-coords. If Ag[0]=NULL we assume a flat 3-metric.*/
void coordinates_set_sqrtgdiag_array(tNode *node, tArray *AdXdx[3][3],
                                     tArray *Ag[6], tArray *Asqrtgdiag[3])
{
  /* arrays to compute sqrtgdiag */
  tArray *adXdxT = alloc_empty_array2d(3,3); /* 3x3 for coord. transf. */
  tArray *ainvM  = alloc_empty_array2d(3,3); /* 3x3 for inv. metric */
  tArray *agam = alloc_array2d(3,3); /* 3x3 for transf. inv. metric */
  tArray *tmp  = alloc_array2d(3,3);
  double dXdx[3][3]; /* coord. transf. */
  double invM[3][3]; /* inverse metric in x-coords  */
  double *gam = Arrd(agam);
  double *sqrtgdiagx = Arrd(Asqrtgdiag[0]);
  double *sqrtgdiagy = Arrd(Asqrtgdiag[1]);
  double *sqrtgdiagz = Arrd(Asqrtgdiag[2]);
  double *pdXdx[3][3];
  double *gxx = Arrd(Ag[0]);
  double *gxy = Arrd(Ag[1]);
  double *gxz = Arrd(Ag[2]);
  double *gyy = Arrd(Ag[3]);
  double *gyz = Arrd(Ag[4]);
  double *gzz = Arrd(Ag[5]);
  double dXbdX[3];
  int ijk, d,e;

  /* init pdXdx: point pdXdx to AdXdx data */
  for(d=0; d<3; d++)
    for(e=0; e<3; e++)
      pdXdx[d][e] = Arrd(AdXdx[d][e]);

  /* set dXbdX from node dimensions */
  dXbYbZb_dXYZ(node, dXbdX);

  /* Put coord. transf. into adXdx.
     Trick here is that C-arrays are row major. So when we put dXdx
     into the col. major adXdxT, we automatically take the transpose! */
  point_array_d_to_data(adXdxT, dXdx, 1);

  /* Put inv. of 3-metric into ainvM. */
  point_array_d_to_data(ainvM, invM, 1);

  /* loop over target arrays,
     we assume all arrays passed in have same dims */
  forarray(Asqrtgdiag[0], ijk)
  {
    /* set transpose of coord. transf. */
    for(d=0; d<3; d++)
      for(e=0; e<3; e++)
        dXdx[d][e] = pdXdx[d][e][ijk];

    /* if no 3-metric is specified, we assume it is flat */
    if(!gxx)
    {
      /* transform flat metric to X coords */
      mm_array0_norestrict(adXdxT,adXdxT, agam);
    }
    else
    {
      double M[3][3] = { { gxx[ijk], gxy[ijk], gxz[ijk] },
                         { gxy[ijk], gyy[ijk], gyz[ijk] },
                         { gxz[ijk], gyz[ijk], gzz[ijk] } };

      /* get inverse metric and multiply with dX/dx */
      inv3Dmat_from_3Dmat(M, invM); /* ainvM points to invM */
      /* gam = dX/dx invM (dX/dx)^T */
      mm_array0(ainvM,adXdxT, tmp); /* ainvM adXdxT -> tmp */
      mm_array0(adXdxT,tmp, agam);  /* adXdx tmp    -> agam */
    }

    /* go from X to Xb coords */
    sqrtgdiagx[ijk] = dXbdX[0] * sqrt(gam[0]);
    sqrtgdiagy[ijk] = dXbdX[1] * sqrt(gam[4]);
    sqrtgdiagz[ijk] = dXbdX[2] * sqrt(gam[8]);
  }
  free_array(tmp);
  free_array(agam);
  free_array(ainvM);
  free_array(adXdxT);
}


/********************************************************************/
/* some functions to convert between pat->xyz_of_XYZ and labels */
/********************************************************************/

/* helper typedef for coord trafo function */
typedef int (*x_of_X)(tPat *pat, tNode *node, int ind,
                      const double X[3], double x[3]);

/* return label of pat->xyz_of_XYZ */
int coordinates_get_label(tPat *pat)
{
  x_of_X xyz_of_XYZ = pat->xyz_of_XYZ;

  /* return label */
  if(xyz_of_XYZ == NULL)
    return Cartesian;
  else if(xyz_of_XYZ == xyz_of_lamAB_CubSph)
    return CubedSphere;
  else if(xyz_of_XYZ == xyz_of_rhoAB_CubSph)
    return CubedSphere_Stretch1;
  else if(xyz_of_XYZ == xyz_of_rh2AB_CubSph)
    return CubedSphere_Stretch2;
  else
    errorexit("implement label for this pat->xyz_of_XYZ");
}

/* return func. pointer xyz_of_XYZ depending on label */
x_of_X coordinates_get_xyz_of_XYZ(int label)
{
  switch(label)
  {
  case Cartesian:
    return NULL;
  case CubedSphere:
    return xyz_of_lamAB_CubSph;
  case CubedSphere_Stretch1:
    return xyz_of_rhoAB_CubSph;
  case CubedSphere_Stretch2:
    return xyz_of_rh2AB_CubSph;
  default:
    errorexiti("unknown or unimplemented label %i", label);
  }
}
