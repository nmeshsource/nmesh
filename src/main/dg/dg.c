/* dg.c */
/* Wolfgang Tichy, April 2019 */


#include "nmesh.h"
#include "dg.h"


/* get glabal amr and coordinate vars */
extern tAMR amr[1];
extern tcoordinates coordinates[1];

/* global vars for dg */
tDGglobals DGglobals[1];

/***************************************************************************/
/* some funcions to add boundary fluxes of discontinous Galerkin (dg) */
/***************************************************************************/

/* Add surface flux terms with a choice of sign (sign=+1 or sign=-1)
   to vldf. We compute the fluxes from vlu.
   If vldf contains the RHS of the DG eqns we need sign=-1, so that the
   surface flux terms get directly added to the RHS.
   If vldf contains div(flux) we need sign=+1, and div(flux) simply gets
   modified by the addition of the surface flux terms. Then something
   else (e.g. vladdto_onfaces_node) has to still subtract the now modified
   vldf from the RHS after this call.
   We also have the flag use_fv, to decide whether fv mode is active at
   all inside this function. */
int dg_add_surface_fluxes_sign_fvflag(tNode *node, double sign,
                          tVarList *vldf, tVarList *vlu, tVarList *vls,
                          void (*u_f_lam)(tDGinfo *d),
                          void (*numflux)(tDGinfo *d),
                          int use_fv)
{
  tDGinfo *dgi = alloc_DGinfo(vlu, vls);
  tMesh *mesh = vlu->mesh;
  double det2g, gdiag;
  int isqrtdet2g_o_det3gamma0 = coordinates->isqrtdet2g_o_det3gamma0;
  int isqrtgdiagx = coordinates->isqrtgdiagx;
  int iooJ        = coordinates->idet_dXbdx;
  int skip_fv = DGglobals->fv_divf_adds_surface_fluxes;
  int add_surface_fluxes = 1; /* by default we want to set fluxes here */
  double distXb[6] = {0};

  if(coordinates->sqrtdet2g_o_det3gamma) { det2g = 1.; gdiag = 0.; }
  else                                   { det2g = 0.; gdiag = 1.; }

  /* we now call get_all_myln_surfaces in evolve_setrhs_mesh
     so we do not need to do it here */
  ///* get surfaces so that we can compute fluxes */
  //get_all_myln_surfaces(mesh);

  /* set overall sign by multiplying det2g and gdiag */
  det2g *= sign;
  gdiag *= sign;

  /* special fv cases */
  if(use_fv)
  {
    /* do nothing if fv_divf has already taken care of surface fluxes */
    if(skip_fv)
    {
      add_surface_fluxes = 0;
    }
    else
    {
      /* find distance from faces to nearest midpoint */
      set_nodemidpoints_to_face_distXb(node, distXb);
    }
  }

  /* add boundary flux terms */
  if(add_surface_fluxes)
  {
    int *n = node->n;
    double *ooJ = Vard(node, iooJ);
    int face;
    int not_fv = !use_fv;
    int have_Xofx = ( node->pat->XYZ_of_xyz ? 1 : 0 );
    int subtr_fi  = DGglobals->fv_flux_is_fnum_minus_fi;
    double mod0 = not_fv;     /* set to 1 if we don't use fin. vol. */
    double mod1 = 1. - mod0;  /* set to 1 if we use fin. vol. */
    /* set s_fi=1 if we subtract inner flux dgi->fi in F */
    double s_fi = ( not_fv || (subtr_fi && have_Xofx) );

    /* set DG node info */
    dgi->node = node;

    for(face=0; face<6; face++)
    {
      int dir = face/2;
      int p = (face%2)*(n[dir] - 1);
      //double sig = 2*(face%2) - 1;
      double *sqrtd2g_o_d3g = Vard(node, isqrtdet2g_o_det3gamma0+face);
      double *sqrtgdiag = Vard(node, isqrtgdiagx+dir);
      double *Wq = Wquad(node, dir);
      double Wqmod = fabs(distXb[face]);
      int i,j,k;

      /* do nothing if dir is not active */
      if(!Getb(amr->dir_active[dir])) continue;

      /* set DG face info */
      dgi->face = face;

      forplaneN(dir, i,j,k, n, p)
      {
        int ijk = Ind_n(i,j,k, n);
        int JK = Ind_n_norm(i,j,k, n, dir);
        int i0 = i0_norm(i,j,k, dir);
        double oow = 1./(Wq[i0]*mod0 + Wqmod*mod1);
        double sdg_oJ_ow = sqrtd2g_o_d3g[JK] * fabs(ooJ[ijk]) * oow;
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

        /* compute numerical flux directly after u_f_lam,
           if not set already in u_f_lam */
        if(numflux) numflux(dgi);

        /* get Ffac, this can be set in u_f_lam or numflux */
        Ffac = dgi->Ffac; /* usually 1, set to 0 to turn off surface fluxes */

        /* get F from dgi and add boundary flux terms to vldf */
        forvl(vldf, l)
        {
          int idf = Vind(vldf,l);
          double *df = Vard_(node, idf);
          double F;

          F = (dgi->fnum[l] - dgi->fi[l]*s_fi) * Ffac;
          df[ijk] += F * (det2g * sdg_oJ_ow + gdiag * gd_ow);
        }
      }
    } /* end loop over faces */
  }

  free_DGinfo(dgi);
  return 0;
}

/* Add surface flux terms with a choice of sign (sign=+1 or sign=-1)
   to vldf. We compute the fluxes from vlu.
   If vldf contains the RHS of the DG eqns we need sign=-1, so that the
   surface flux terms get directly added to the RHS.
   If vldf contains div(flux) we need sign=+1, and div(flux) simply gets
   modified by the addition of the surface flux terms. Then something
   else (e.g. vladdto_onfaces_node) has to still subtract the now modified
   vldf from the RHS after this call. */
int dg_add_surface_fluxes_sign(tNode *node, double sign, tVarList *vldf,
                               tVarList *vlu, tVarList *vls,
                               void (*u_f_lam)(tDGinfo *d),
                               void (*numflux)(tDGinfo *d))
{
  int use_fv = node->dat->info->use_fv;
  return dg_add_surface_fluxes_sign_fvflag(node, sign, vldf, vlu, vls,
                                           u_f_lam, numflux, use_fv);
}

/* add surface flux terms of DG formulation to RHS */
int dg_add_surface_fluxes(tNode *node, tVarList *vlr, tVarList *vlu,
                          tVarList *vls,
                          void (*u_f_lam)(tDGinfo *d),
                          void (*numflux)(tDGinfo *d))
{
  return dg_add_surface_fluxes_sign(node, -1., vlr, vlu, vls,
                                    u_f_lam, numflux);
}


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

  /* extra space for source terms */
  if(vls)
  {
    dgi->si = dmalloc(vls->n); /* source vars inside this node */
    dgi->sa = dmalloc(vls->n); /* source vars on adjacent side */
  }
  else
  {
    dgi->si = dgi->sa = NULL;
  }

  return dgi;
}

/* free DGinfo structure */
void free_DGinfo(tDGinfo *dgi)
{
  /* free contents */
  free(dgi->sa);
  free(dgi->si);

  free(dgi->fnum);

  free(dgi->lama);
  free(dgi->fa);
  free(dgi->ua);

  free(dgi->lami);
  free(dgi->fi);
  free(dgi->ui);

  /* free dgi */
  free(dgi);
}


/* copy the parts of struct tDGINFO that are not allocated by alloc_DGinfo
   from dsrc to ddest */
void copy_nonallocd_DGinfo(tDGinfo *dsrc, tDGinfo *ddest)
{
  tDGinfo dbak[1];

  /* shallow backup copy of ddest */
  *dbak = *ddest;

  /* shallow copy of dsrc, to get everthing from dsrc into ddest */
  *ddest = *dsrc;

  /* now restore the allocd parts of ddest from dbak */
  ddest->ui   = dbak->ui;
  ddest->fi   = dbak->fi;
  ddest->lami = dbak->lami;
  ddest->ua   = dbak->ua;
  ddest->fa   = dbak->fa;
  ddest->lama = dbak->lama;

  ddest->fnum = dbak->fnum;

  ddest->si   = dbak->si;
  ddest->sa   = dbak->sa;
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
  int fv_divf_extrap = Par("fv_divf_extrap");
  int fv_surface_interp = Par("fv_surface_interp");
  double WENOm3_opt_weightratio = Getd(Par("fv_WENOm3_opt_weightratio"));
  char *list, *saveptr, *name;

  /* set reconstruction mode */
  if(Getv(fv_rec, "1"))
    DGglobals->fv_rec_mode = FV_REC_1;
  else if(Getv(fv_rec, "WENOm3_2"))
    DGglobals->fv_rec_mode = FV_REC_WENOm3_2;
  else if(Getv(fv_rec, "WENOm5_2"))
    DGglobals->fv_rec_mode = FV_REC_WENOm5_2;
  else if(Getv(fv_rec, "WENOmT_2"))
    DGglobals->fv_rec_mode = FV_REC_WENOmT_2;
  else if(Getv(fv_rec, "WENOmZ_2"))
    DGglobals->fv_rec_mode = FV_REC_WENOmZ_2;
  else if(Getv(fv_rec, "WENOm3_1"))
    DGglobals->fv_rec_mode = FV_REC_WENOm3_1;
  else if(Getv(fv_rec, "WENOm5_1"))
    DGglobals->fv_rec_mode = FV_REC_WENOm5_1;
  else if(Getv(fv_rec, "WENOmZ_1"))
    DGglobals->fv_rec_mode = FV_REC_WENOmZ_1;
  else if(Getv(fv_rec, "WENOm3"))
    DGglobals->fv_rec_mode = FV_REC_WENOm3;
  else if(Getv(fv_rec, "WENOm5"))
    DGglobals->fv_rec_mode = FV_REC_WENOm5;
  else if(Getv(fv_rec, "WENOmZ"))
    DGglobals->fv_rec_mode = FV_REC_WENOmZ;
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

  /* WENOm3,WENOm5,WENOmZ use 1st or 2nd order acc. in rec. near boundary */
  if( (DGglobals->fv_rec_mode == FV_REC_WENOm3) ||
      (DGglobals->fv_rec_mode == FV_REC_WENOm5) ||
      (DGglobals->fv_rec_mode == FV_REC_WENOmZ) )
    DGglobals->fv_rec_mode_WENOm = 1; //need to decide on 1st or 2nd order
  DGglobals->fv_rec_WENOm_s1  = Getd(Par("fv_rec_WENOm_s1"));
  DGglobals->fv_rec_WENOm_s2  = Getd(Par("fv_rec_WENOm_s2"));
  DGglobals->fv_rec_WENOm_opt = Getd(Par("fv_rec_WENOm_opt"));

  /* save how we get fv flux */
  if(Getv(Par("fv_flux"), "fnum_minus_fi"))
    DGglobals->fv_flux_is_fnum_minus_fi = 1;

  /* set extrapolation mode for div(flux) */
  if(Getv(fv_divf_extrap, "no"))
  {
    DGglobals->fv_divf_extrap_mode = FV_NO_EXTRAP;
    DGglobals->fv_divf_adds_surface_fluxes = 0;
  }
  else if(Getv(fv_divf_extrap, "divf_extrap1"))
  {
    DGglobals->fv_divf_extrap_mode = FV_DIVF_EXTRAP1;
    DGglobals->fv_divf_adds_surface_fluxes = 0;
  }
  else if(Getv(fv_divf_extrap, "dnfn_extrap1"))
  {
    DGglobals->fv_divf_extrap_mode = FV_DNFN_EXTRAP1;
    DGglobals->fv_divf_adds_surface_fluxes = 1;
  }
  else
  {
    errorexits("unknown value %s in par fv_divf_extrap.",
               Gets(fv_divf_extrap));
  }

  /* set surface interpolation mode */
  if(Getv(fv_surface_interp,      "linear"))
    DGglobals->fv_surface_interp_mode = FV_2DINTERP_LINEAR;
  else if(Getv(fv_surface_interp, "parabolic"))
    DGglobals->fv_surface_interp_mode = FV_2DINTERP_PARAB;
  else
    errorexits("unknown value %s in par fv_surface_interp.",
               Gets(fv_surface_interp));

  /* set optimal weights for WENOm3 */
  DGglobals->fv_WENOm3_optw[0] = 1.;
  DGglobals->fv_WENOm3_optw[1] = WENOm3_opt_weightratio;

  /* set flux factors for outer BCs */
  if( sscanf(Gets(Par("dg_outerBC_flux_fac")), "%lg %lg %lg",
             &(DGglobals->outerBC_flux_fac[0]),
             &(DGglobals->outerBC_flux_fac[1]),
             &(DGglobals->outerBC_flux_fac[2])) != 3 )
    errorexit("par dg_outerBC_flux_fac must contain 3 numbers");

  /* set some more par values */
  DGglobals->fv_divf_extrap_s1 = Getd(Par("fv_divf_extrap_s1"));
  DGglobals->fv_divf_extrap_s2 = Getd(Par("fv_divf_extrap_s2"));
  DGglobals->fv_divf_extrap_opt = Geti(Par("fv_divf_extrap_opt"));
  DGglobals->fv_divf_use_only_right_flux
    = Getb(Par("fv_divf_use_only_right_flux"));

  /* info from fv2dg_interp_use_extrap1 pars */
  DGglobals->fv2dg_interp_use_extrap1 = Getb(Par("fv2dg_interp_use_extrap1"));
  vlfree(DGglobals->fv2dg_interp_use_extrap1_vl);
  DGglobals->fv2dg_interp_use_extrap1_vl = vlalloc(mesh);
  list = strdup(Gets(Par("fv2dg_interp_use_extrap1_vars")));
  for(name=strtok_r(list, " ", &saveptr); name!=NULL;
      name=strtok_r(NULL, " ", &saveptr))
  {
    vlpush(DGglobals->fv2dg_interp_use_extrap1_vl, Ind(name));
  }
  free(list);

  return 0;
}

/* free any memory that was allocated in dg_set_DGglobals */
int dg_free_DGglobals(tMesh *mesh)
{
  vlfree(DGglobals->fv2dg_interp_use_extrap1_vl);
  DGglobals->fv2dg_interp_use_extrap1_vl = NULL;
  return 0;
}

/* print what we have in DGglobals */
int dg_print_DGglobals(tMesh *mesh)
{
  int d;

  PRFs(":\n");
  printf(" DGglobals->outerBC_flux_fac = {");
  for(d=0; d<3; d++) printf(" %.16g", DGglobals->outerBC_flux_fac[d]);
  printf(" }\n");
  printf(" DGglobals->fv_rec_mode = %d\n", DGglobals->fv_rec_mode);
  printf(" DGglobals->fv_rec_mode_WENOm = %d\n",
         DGglobals->fv_rec_mode_WENOm);
  printf(" DGglobals->fv_rec_WENOm_s1 = %g\n",
         DGglobals->fv_rec_WENOm_s1);
  printf(" DGglobals->fv_rec_WENOm_s2 = %g\n",
         DGglobals->fv_rec_WENOm_s2);
  printf(" DGglobals->fv_rec_WENOm_opt = %d\n",
         DGglobals->fv_rec_WENOm_opt);
  printf(" DGglobals->fv_flux_is_fnum_minus_fi = %d\n",
         DGglobals->fv_flux_is_fnum_minus_fi);
  printf(" DGglobals->fv_divf_extrap_mode = %d\n",
         DGglobals->fv_divf_extrap_mode);
  printf(" DGglobals->fv_divf_adds_surface_fluxes = %d\n",
         DGglobals->fv_divf_adds_surface_fluxes);
  printf(" DGglobals->fv_surface_interp_mode = %d\n",
         DGglobals->fv_surface_interp_mode);
  printf(" DGglobals->fv_WENOm3_optw = { %g %g }\n",
         DGglobals->fv_WENOm3_optw[0], DGglobals->fv_WENOm3_optw[1]);
  printf(" DGglobals->fv_divf_extrap_s1 = %g\n",
         DGglobals->fv_divf_extrap_s1);
  printf(" DGglobals->fv_divf_extrap_s2 = %g\n",
         DGglobals->fv_divf_extrap_s2);
  printf(" DGglobals->fv_divf_extrap_opt = %d\n",
         DGglobals->fv_divf_extrap_opt);
  printf(" DGglobals->fv_divf_use_only_right_flux = %d\n",
         DGglobals->fv_divf_use_only_right_flux);

  printf(" DGglobals->fv2dg_interp_use_extrap1 = %d\n",
         DGglobals->fv2dg_interp_use_extrap1);
  printf(" DGglobals->fv2dg_interp_use_extrap1_vl: ");
  prvarlist(DGglobals->fv2dg_interp_use_extrap1_vl);

  return 0;
}


/* get normal at gridpoint or midpoint to left or right of gridpoint ijk
   depending on info in tDGinfo *d */
double node_normal_from_DGinfo(tDGinfo *d, double nrm[3])
{
  tNode *node = d->node;
  int ijk = Ind_n(d->i,d->j,d->k, node->n);
  int f = d->face;

  /* get normal at midpoint next to ijk */
  if(d->info & DGINFO_MIDPTNORM)
    return node_normal_at_midpt_nextto_ijk(node, f, ijk, nrm);
  else /* or rather on grid point ijk */
    return node_normal_at_ijk(node, f, ijk, nrm);
}
