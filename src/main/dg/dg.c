/* dg.c */
/* Wolfgang Tichy, April 2019 */


#include "nmesh.h"
#include "dg.h"


/* the dgs in here use the struct tINDIC within tDat to exchange data
   between neighbors. The exchange is done in evolve, which calls functions
   from main/amr/indicators.c to do the exchange. */


/* set data for MRS dg. This one uses the cons vars in vl */




/* use numerical flux FN^i to set F */
void dg_advection1_F(tMesh *mesh, tVarList *vlu)
{
  int iu = vlu->index[0];
  int ifx = Ind("advection1_fx");
  int iF  = Ind("advection1_F0");
  char *advdir = Gets(Par("advection1_direction"));
  double nx,ny,nz;
  int myid;

  /* prop. dir.*/
  sscanf(advdir, "%lg %lg %lg", &nx, &ny, &nz);

  /* compute boundary flux terms */
  formylnodes(mesh, myid)
  {
    tNode *node = MyNode(mesh, myid);
    int *n = node->n;
    double *fx = Vard(node, ifx);
    double *fy = Vard(node, ifx+1);
    double *fz = Vard(node, ifx+2);
    double FNx,FNy,FNz, norm[3];
    int face, dir, p, i,j,k, ijk, JK;

    /* set F on each face */
    for(face=0; face<6; face++)
    {
      double *F = Vard(node, iF+face);
      double *uaj = Varaj(node, iu, face);
      dir = face/2;
      p = (face%2)*(n[dir] - 1);
      forplaneN(dir, i,j,k, n, p)
      {
        ijk = Ind_n(i,j,k, n);
        JK = Ind_n_norm(i,j,k, n, dir);
        node_normal_at_ijk(node, face, ijk, norm);

        /* if stuff is coming in */
        if(norm[0]*nx + norm[1]*ny + norm[2]*nz < 0.)
        {
          /* if there is an adjacent surface */
          if(uaj)
          {
            FNx = uaj[JK] * nx;
            FNy = uaj[JK] * ny;
            FNz = uaj[JK] * nz;
          }
          else
          {
            FNx = FNy = FNz = 0;
          }
        }
        else
        {
          FNx = fx[ijk];
          FNy = fy[ijk];
          FNz = fz[ijk];
        }

        F[JK] = (FNx - fx[ijk])*norm[0] +
                (FNy - fy[ijk])*norm[1] +
                (FNz - fz[ijk])*norm[2];
      }
    }
  }
}


/* add surface flux terms */
int dg_surface_fluxes(tMesh *mesh, tVarList *vlr, tVarList *vlu)
{
  int ir = vlr->index[0];
  //int iu = vlu->index[0];
  int iooJ = Ind("det_dXbdx");
  int isqrtdet2gamma0 = Ind("sqrtdet2gamma0");
  int myid;

  /* get surfaces so that we can compute fluxes */
  get_all_myln_surfaces(mesh);

  /* get flux terms on boundary */
  dg_advection1_F(mesh, vlu);

int iF   = Ind("advection1_F0");

  /* add boundary flux terms */
  formylnodes(mesh, myid)
  {
    tNode *node = MyNode(mesh, myid);
    int *n = node->n;
    double *ooJ = Vard(node, iooJ);
    int face;

    for(face=0; face<6; face++)
    {
      int dir = face/2;
      int p = (face%2)*(n[dir] - 1);
      //double sig = 2*(face%2) - 1;
      double *sqrtdet2gam = Vard(node, isqrtdet2gamma0+face);
      double *F = Vard(node, iF+face);
      double *w = Wquad(node, dir);
      int i,j,k, ijk, JK, i0;

      forplaneN(dir, i,j,k, n, p)
      {
        ijk = Ind_n(i,j,k, n);
        JK = Ind_n_norm(i,j,k, n, dir);
        i0 = i0_norm(i,j,k, dir);

        forvl(vlr, ir)
        {
          double *r = Vard(node, ir);
          r[ijk] -= F[JK] * sqrtdet2gam[JK] * fabs(ooJ[ijk])/ w[i0];
        }
      }
    }
  }

  return 0;
}
