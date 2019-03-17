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
  int iX   = Ind("X");
  int idXd = Ind("dXdx");
  int idet_dXbdx = Ind("det_dXbdx");
  int ix   = Ind("x");
  int f, d;

  /* do nothing if this is not my node */
  if(!dat) return 0;

  /* if coords are set already do nothing */
  if(dat->coords_set) return 1;

  /* give all these memory */
  enablevar_innode(node, iX);
  enablevar_innode(node, iX+1);
  enablevar_innode(node, iX+2);
  enablevar_innode(node, idXd);
  enablevar_innode(node, idXd+3);
  enablevar_innode(node, idXd+6);
  enablevar_innode(node, idet_dXbdx);
  enablevar_innode(node, ix);
  enablevar_innode(node, ix+1);
  enablevar_innode(node, ix+2);

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
  int iX   = Ind("X");
  int idXd = Ind("dXdx");
  int ix   = Ind("x");
  double *pX[] = { Vard(node,iX), Vard(node,iX+1), Vard(node,iX+2) };
  double *px[] = { Vard(node,ix), Vard(node,ix+1), Vard(node,ix+2) };
  double *pdXd[3][3]
            = { {Vard(node,idXd),   Vard(node,idXd+1), Vard(node,idXd+2)},
                {Vard(node,idXd+3), Vard(node,idXd+4), Vard(node,idXd+5)},
                {Vard(node,idXd+6), Vard(node,idXd+7), Vard(node,idXd+8)} };
  double *det_dXbdx = Vard(node, Ind("det_dXbdx"));
  double dXbdX[3];
  double det_dXbYbZb_dXYZ;

  /* do nothing if coords are set already or if vars are off */
  if(!dat) return 0;
  if(dat->coords_set) return 0;
  if(!vars_on) return 0;

  PRF;printf(": nid%ld\n", node->nid);

  /* get det of dXb/dX */
  dXbYbZb_dXYZ(node, dXbdX);
  det_dXbYbZb_dXYZ = dXbdX[0] * dXbdX[1] * dXbdX[2];

  /* set coords */
  forijk(i,j,k, n)
  {
    double Xb[] = { node->Xb[0]->d[i], node->Xb[1]->d[j], node->Xb[2]->d[k] };
    double X[3], x[3], dXd[3][3];
    int ijk = Ind_n(i,j,k, n);

    /* get X from Xb */
    XYZ_of_XbYbZb(node, Xb, X);
    for(d=0; d<3; d++) pX[d][ijk] = X[d];

    /* now set x, dXdx, det(dXb/dx) */
    if(pat->dXYZ_dxyz)
    {
      pat->dXYZ_dxyz(pat, node, -1, X, x, dXd);
      for(d=0; d<3; d++)
      {
        px[d][ijk] = x[d];
        for(e=0; e<3; e++) pdXd[d][e][ijk] = dXd[d][e];
      }
      det_dXbdx[ijk] = det_dXbYbZb_dXYZ * det_3Dmatrix(dXd);
    }
    else /* assume X,Y,Z are Cartesian*/
    {
      for(d=0; d<3; d++)
      {
        px[d][ijk] = pX[d][ijk];
        pdXd[d][d][ijk] = 1.;
      }
      det_dXbdx[ijk] = det_dXbYbZb_dXYZ;
    }
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
