/* coordinates.c */
/* Wolfgang Tichy, 2/2019 */

#include "nmesh.h"
#include "coordinates.h"

#define PR 1



/* try to enable all coord vars, return 1 if we have dat */
int coordinates_coordvars_enabled(tNode *node)
{
  tPat *pat = node->pat;
  tMesh *mesh = pat->mesh;
  tDat *dat = node->dat;
  tCoordInfo *CI = pat->CI;
  int iX, ix, idXdx, idet_dXbdx;
  int isqrtdet2gamma0, isqrtgdiagx;
  int surface_metric, sqrtdet2gamma, sqrtgdiag;
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
  isqrtdet2gamma0 = Ind("sqrtdet2gamma0");
  isqrtgdiagx = Ind("sqrtgdiagx");

  /* which surface info do we set */
  surface_metric = Par("coordinates_surface_metric");
  sqrtdet2gamma  = Getv(surface_metric, "sqrtdet2gamma");
  sqrtgdiag      = Getv(surface_metric, "sqrtgdiag");

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
  if(1 || sqrtdet2gamma) enablevar_innode(node, isqrtdet2gamma0);
  if(1 || sqrtgdiag)     enablevar_innode(node, isqrtgdiagx);

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
  int surface_metric, sqrtdet2gamma, sqrtgdiag;
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
  double *gxx, *gxy, *gxz, *gyy, *gyz, *gzz;

  /* do nothing if coords are set already or if vars are off */
  if(!dat) return 0;
  if(dat->coords_set) return 0;
  if(!vars_on) return 0;

  PRF;printf(": nid%ld\n", node->nid);

  /* pointers to 3 metric */
  if(i3metric>=0)
  {
    gxx = Vard(node, i3metric);
    gxy = Vard(node, i3metric+1);
    gxz = Vard(node, i3metric+2);
    gyy = Vard(node, i3metric+3);
    gyz = Vard(node, i3metric+4);
    gzz = Vard(node, i3metric+5);
  }
  else
  { gxx = gxy = gxz = gyy = gyz = gzz = NULL; }

  /* which surface info do we set */
  surface_metric = Par("coordinates_surface_metric");
  sqrtdet2gamma  = Getv(surface_metric, "sqrtdet2gamma");
  sqrtgdiag      = Getv(surface_metric, "sqrtgdiag");

  /* get det of dXb/dX */
  dXbYbZb_dXYZ(node, dXbdX);
  det_dXbYbZb_dXYZ = dXbdX[0] * dXbdX[1] * dXbdX[2];

  /* set coords */
  forijk(i,j,k, n)
  {
    double Xb[] = { node->Xb[0]->d[i], node->Xb[1]->d[j], node->Xb[2]->d[k] };
    double X[3], x[3], dXdx[3][3];
    int ijk = Ind_n(i,j,k, n);

    /* get X from Xb */
    XYZ_of_XbYbZb(node, Xb, X);
    for(d=0; d<3; d++) pX[d][ijk] = X[d];

    /* now set x, dXdx, det(dXb/dx) */
    if(pat->dXYZ_dxyz)
    {
      pat->dXYZ_dxyz(pat, node, -1, X, x, dXdx);
      if(0) if( !isfinite(x[0]) || !isfinite(x[1]) || !isfinite(x[2]) )
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

  /* set sqrtdet2gamma on node faces */
  if(sqrtdet2gamma)
  {
    int isqrtdet2gamma0 = Ind("sqrtdet2gamma0");
    /* arrays to compute 2 metric sqrtdet2gamma on faces */
    tArray *a3gamT = alloc_empty_array2d(3,3); /* 3x3 for transp. of 3-metric */
    tArray *a2J = alloc_array2d(3,2);   /* 3x2 for 2-Jacobian */
    tArray *a2gam = alloc_array2d(2,2); /* 2x2 for induced metric on surface */
    tArray *tmp = alloc_array2d(3,2);

    for(f=0; f<6; f++)
    {
      int dir = f/2;
      int p = (f%2)*(n[dir] - 1);
      double *sqrtdet2gam = Vard(node, isqrtdet2gamma0 + f);

      forplaneN(dir, i,j,k, n, p)
      {
        double det2gam;
        int ijk = Ind_n(i,j,k, n);
        int JK = Ind_n_norm(i,j,k, n, dir);

        /* get 2-Jacobian of dx/dXb, which is a 3x2 matrix */
        array_2dxdXb(node, ijk, dir, a2J);

        /* if no 3-metric is specified, we assume it is flat */
        if(i3metric<0)
        {
          /* compute 2-metric from a2J^T a2J */
          mm_array0_norestrict(a2J,a2J, a2gam);
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
          point_array_a_to_data(a3gamT, M, 1);
          mm_array0(a3gamT,a2J, tmp); /* a3gam a2J -> tmp */
          mm_array0(a2J,tmp, a2gam);  /* a2J^T tmp -> a2gam */
        }
        det2gam = det_2_2_array(a2gam);
        //NOTE: there should be a faster way to get det2gam. There should
        //      be an analogue to \det(g)=\alpha^2\det(\gamma)
        sqrtdet2gam[JK] = sqrt(det2gam);
      }
    }
    /* free arrays */
    free_array(tmp);
    free_array(a2gam);
    free_array(a2J);
    free_array(a3gamT);
  }

  /* set sqrtgdiag */
  if(sqrtgdiag)
  {
    /* arrays to compute sqrtgdiag */
    tArray *adXdxT = alloc_empty_array2d(3,3); /* 3x3 for coord. transf. */
    tArray *ainvM  = alloc_empty_array2d(3,3); /* 3x3 for inv. metric */
    tArray *agam = alloc_array2d(3,3); /* 3x3 for transf. inv. metric */
    tArray *tmp  = alloc_array2d(3,3);
    double dXdx[3][3]; /* coord. transf. */
    double invM[3][3]; /* inverse metric in x-coords  */
    double *gam = Arrd(agam);
    int isqrtgdiagx = Ind("sqrtgdiagx");
    double *sqrtgdiagx = Vard(node, isqrtgdiagx);
    double *sqrtgdiagy = Vard(node, isqrtgdiagx+1);
    double *sqrtgdiagz = Vard(node, isqrtgdiagx+2);
    int ijk;

    /* Put coord. transf. into adXdx.
       Trick here is that C-arrays are row major. So when we put dXdx
       into the col. major adXdxT, we automatically take the transpose! */
    point_array_a_to_data(adXdxT, dXdx, 1);

    /* Put inv. of 3-metric into ainvM. */
    point_array_a_to_data(ainvM, invM, 1);

    forpoints(node, ijk)
    {
      /* set transpose of coord. transf. */
      for(d=0; d<3; d++)
        for(e=0; e<3; e++)
          dXdx[d][e] = pdXdx[d][e][ijk];

      /* if no 3-metric is specified, we assume it is flat */
      if(i3metric<0)
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

  /* mark coords as set */
  dat->coords_set = 1;

  return 0;
}

/* initialize coordinates in each patch */
int coordinates_init(tMesh *mesh)
{
  int myid;

  PRF;printf(":\n");

  formylnodes(mesh, myid)
  {
    tNode *node = MyNode(mesh, myid);
    coordinates_init_node(node);
  }

  return 0;
}
