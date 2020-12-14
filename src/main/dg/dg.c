/* dg.c */
/* Wolfgang Tichy, April 2019 */


#include "nmesh.h"
#include "dg.h"


/* global vars for dg */
tDGglobals DGglobals[1];


/* some funcions to add boundary fluxes of discontinous Galerkin (dg) */

/* allocate DGinfo structure */
tDGinfo *alloc_DGinfo(tVarList *vlu, tVarList *vls)
{
  tDGinfo *dgi = calloc(1, sizeof(dgi[0]));
  int nvars = vlu->n;

  /* set varlists */
  dgi->vlu  = vlu; /* varlist with cons vars */
  dgi->vls  = vls; /* varlist with needed source terms, could be NULL */

  /* alloc mem for vars at point i,j,k */
  dgi->ui   = dmalloc(nvars); /* cons. vars inside this node */
  dgi->fi   = dmalloc(nvars);
  dgi->lami = dmalloc(nvars);
  dgi->ua   = dmalloc(nvars); /* cons. vars on adjacent side */
  dgi->fa   = dmalloc(nvars);
  dgi->lama = dmalloc(nvars);

  dgi->fnum = dmalloc(nvars);
  dgi->Ffac = 1.;

  return dgi;
}

/* free DGinfo structure */
void free_DGinfo(tDGinfo *dgi)
{
  /* free contents */
  free(dgi->ui);
  free(dgi->fi);
  free(dgi->lami);
  free(dgi->ua);
  free(dgi->fa);
  free(dgi->lama);
  free(dgi->fnum);

  /* free dgi */
  free(dgi);
}


/* add surface flux terms */
int dg_add_surface_fluxes(tMesh *mesh, tVarList *vlr, tVarList *vlu,
                          tVarList *vls,
                          void (*u_f_lam)(tDGinfo *d),
                          void (*numflux)(tDGinfo *d))
{
  int surface_metric = Par("coordinates_surface_metric");
  double det2gam     = Getv(surface_metric, "sqrtdet2gamma");
  double gdiag       = Getv(surface_metric, "sqrtgdiag");
  int isqrtdet2gamma0 = Ind("sqrtdet2gamma0");
  int isqrtgdiagx     = Ind("sqrtgdiagx");
  int iooJ = Ind("det_dXbdx");

  TIMER_START;

  /* we now call get_all_myln_surfaces in evolve_setrhs_mesh
     so we do not need to do it here */
  ///* get surfaces so that we can compute fluxes */
  //get_all_myln_surfaces(mesh);

  /* loop over nodes so we can add boundary flux terms */
  NODELEVEL_Pragma(omp parallel)
  {
    tDGinfo *dgi = alloc_DGinfo(vlu, vls); /* each thread gets its own dgi */

    /* Note the following leaf node loop cannot be a taskloop because we
       allocated one dgi per thread. But each task needs its own dgi!!!
       And one thread may do more than one task... */
    formylnodes_ompfor(mesh)
    {
      tNode *node = MyLnode;
      int *n = node->n;
      double *ooJ = Vard(node, iooJ);
      int face;
      int use_fv = node->dat->info->use_fv;
      double mod0 = (!use_fv);  /* set to 1 if we don't use fin. vol. */
      double mod1 = 1. - mod0;  /* set to 1 if we use fin. vol. */
      double distXb[6];

      /* find distance from faces to nearest midpoint */
      set_nodemidpoints_to_face_distXb(node, distXb);

      /* set DG node info */
      dgi->node = node;

      for(face=0; face<6; face++)
      {
        int dir = face/2;
        int p = (face%2)*(n[dir] - 1);
        //double sig = 2*(face%2) - 1;
        double *sqrtdet2gam = Vard(node, isqrtdet2gamma0+face);
        double *sqrtgdiag = Vard(node, isqrtgdiagx+dir);
        double *Wq = Wquad(node, dir);
        double Wqmod = fabs(distXb[face]);
        int i,j,k;

        /* set DG face info */
        dgi->face = face;

        forplaneN(dir, i,j,k, n, p)
        {
          int ijk = Ind_n(i,j,k, n);
          int JK = Ind_n_norm(i,j,k, n, dir);
          int i0 = i0_norm(i,j,k, dir);
          double oow = 1./(Wq[i0]*mod0 + Wqmod*mod1);
          double sdg_oJ_ow = sqrtdet2gam[JK] * fabs(ooJ[ijk]) * oow;
          double gd_ow = sqrtgdiag[ijk] * oow;
          double Ffac;
          int l;

          /* set DG i,j,k info */
          dgi->i = i;
          dgi->j = j;
          dgi->k = k;
          /* we do not use:  dgi->info = use_fv;
             because our regular grid points on the faces already are
             considered midpoints in fin. vol. approach */

          /* set vars, fluxes and eigenvals on both sides */
          u_f_lam(dgi);

          /* compute numerical flux */
          numflux(dgi);

          /* get Ffac, this can be set in u_f_lam or numflux */
          Ffac = dgi->Ffac; /* usually 1, set to 0 to turn off surface fluxes */

          /* get F from dgi and add boundary flux terms to RHS */
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
            //  if(dgi->lami[l] < 0.) dgi->fnum[l] = 0.;
            //  else                  dgi->fnum[l] = dgi->fi[l];
            //}

            F = (dgi->fnum[l] - dgi->fi[l]*mod0) * Ffac;
            r[ijk] -= F * (det2gam * sdg_oJ_ow + gdiag * gd_ow);
//if(l==0 && face==1)
//{
//printDGinfo(dgi);
//printf("nid%ld dgi->fnum[l]=%g gd_ow=%g dgi->fnum[l]*gd_ow=%g\n",
//node->nid, dgi->fnum[l], gd_ow, dgi->fnum[l]*gd_ow);
//}
          }
        }
      } /* end loop over faces */
    }
    free_DGinfo(dgi);
  }

  TIMER_STOP;

  return 0;
}

/* print */
void printDGinfo(tDGinfo *d)
{
  int k;

  PRFs(": ");pr_nodename(d->node);
  printf(": face=%d i,j,k=%d,%d,%d Ffac=%g info=%d\n",
         d->face, d->i,d->j,d->k, d->Ffac, d->info);
  if(d->vlu) prvarlist(d->vlu);
  if(d->vls) prvarlist(d->vls);
  forvl(d->vlu, k)
  {
    printf(" %3d: ui=%g fi=%g lami=%g\n", k, d->ui[k], d->fi[k], d->lami[k]);
    printf("      ua=%g fa=%g lama=%g", d->ua[k], d->fa[k], d->lama[k]);
    printf(" => fnum=%g\n", d->fnum[k]);
  }
}


/* init DGglobals struct */
int dg_set_DGglobals(tMesh *mesh)
{
  int fv_rec = Par("fv_rec");

  if(Getv(fv_rec, "1"))
    DGglobals->fv_rec_mode = FV_REC_1;
  else if(Getv(fv_rec, "WENOm3_2"))
    DGglobals->fv_rec_mode = FV_REC_WENOm3_2;
  else if(Getv(fv_rec, "WENO3if1away_1"))
    DGglobals->fv_rec_mode = FV_REC_WENO3if1away_1;
  else if(Getv(fv_rec, "WENO3if2away_1"))
    DGglobals->fv_rec_mode = FV_REC_WENO3if2away_1;
  else if(Getv(fv_rec, "WENO3_2"))
    DGglobals->fv_rec_mode = FV_REC_WENO3_2;
  else if(Getv(fv_rec, "WENO3_2g"))
    DGglobals->fv_rec_mode = FV_REC_WENO3_2g;
  else
    errorexits("unknown value %s in par fv_rec.", Gets(fv_rec));

  return 0;
}
