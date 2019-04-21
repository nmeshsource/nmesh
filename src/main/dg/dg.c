/* dg.c */
/* Wolfgang Tichy, April 2019 */


#include "nmesh.h"
#include "dg.h"


/* some funcions to add boundary fluxes of discontinous Galerkin (dg) */


/* add surface flux terms */
int dg_add_surface_fluxes(tMesh *mesh, tVarList *vlr, tVarList *vlu,
                          int (*u_f_lam)(tNode *node, int face,
                                         int i, int j, int k,
                                         tVarList *vlu,
                                         double *ui, double *ua,
                                         double *fi,  double *fa,
                                         double *lami, double *lama))
{
  int nvars = vlu->n;
  double *ui   = dmalloc(nvars); /* cons. vars inside this node */
  double *fi   = dmalloc(nvars);
  double *lami = dmalloc(nvars);
  double *ua   = dmalloc(nvars); /* cons. vars on adjacent side */
  double *fa   = dmalloc(nvars);
  double *lama = dmalloc(nvars);
  double *fnum = dmalloc(nvars);
  int iooJ = Ind("det_dXbdx");
  int isqrtdet2gamma0 = Ind("sqrtdet2gamma0");
  int myid;

  /* get surfaces so that we can compute fluxes */
  get_all_myln_surfaces(mesh);

  /* loop over nodes so we can add boundary flux terms */
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
      double *w = Wquad(node, dir);
      int i,j,k;

      forplaneN(dir, i,j,k, n, p)
      {
        int ijk = Ind_n(i,j,k, n);
        int JK = Ind_n_norm(i,j,k, n, dir);
        int i0 = i0_norm(i,j,k, dir);
        int l;

        /* set vars, fluxes and eigenvals on both sides */
        u_f_lam(node, face, i,j,k, vlu, ui,ua, fi,fa, lami,lama);

        /* compute numerical flux */
        numflux1d_LLF(mesh, nvars, fnum, ui,ua, fi,fa, lami,lama);
        //numflux1d_upwind(mesh, nvars, fnum, ui,ua, fi,fa, lami,lama);

if(0 && myid==4 && face==1)
{
printf("fnum=%g: ui=%g ua=%g fi=%g fa=%g lami=%g lama=%g\n",
*fnum, *ui,*ua, *fi,*fa, *lami,*lama);
}
        /* add boundary flux terms to RHS */
        forvl(vlr, l)
        {
          int ir = Vind(vlr,l);
          double *r = Vard_(node, ir);
          double F;
          //int iu = Vind(vlu,l);
          //double *uaj = Varaj(node, iu, face);

          /* do something special on outer boundary */
          //if(!uaj)
          //{
          //  if(lami[l] < 0.) fnum[l] = 0.;
          //  else             fnum[l] = fi[l];
          //}

          F = fnum[l] - fi[l];
          r[ijk] -= F * sqrtdet2gam[JK] * fabs(ooJ[ijk])/ w[i0];
        }
      }
    }
  }

  free(fnum);
  free(lama);
  free(fa);
  free(ua);
  free(lami);
  free(fi);
  free(ui);

  return 0;
}
