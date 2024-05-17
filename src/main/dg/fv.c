/* fv.c */
/* Wolfgang Tichy, Dec 2020 */


#include "nmesh.h"
#include "dg.h"

/* get glabal amr vars */
extern tAMR amr[1];

/* use DGglobals */
extern tDGglobals DGglobals[1];
extern tcoordinates coordinates[1];

/***********************************************************************/
/* funcs needed for finite volume method in nmesh */
/***********************************************************************/

/* Function that reconstructs q-vars at midpoint im,
   directly from arrays qc that contain q-vars at gridpoints.
   Note: fv->qc[l][i], qm_p[l], qm_m[l] are allocated in fv_divf.
   In: nq, qc, npts, im, q_scale, rec1d_p,rec1d_m. Out: qm_p, qm_m */
void fv_rec1d_q_midpt(tFVinfo *fv)
{
  int nq      = fv->nq;
  double **qc = fv->qc;   // qc[0..nvars-1][0..npts-1]
  int npts    = fv->npts;
  int im      = fv->im;   // im = 0..npts-2, im is midpt to right of grdpt im
  double q_scale = fv->q_scale;
  int l;

  /* interpolate fields qc towards the midpoint at im */
  for(l=0; l<nq; l++)
  {
    /* reconstruct from both sides of midpoint at im */
    fv->qm_p[l] = fv->rec1d_p(npts, qc[l], im, q_scale);
    fv->qm_m[l] = fv->rec1d_m(npts, qc[l], im, q_scale);
  }
}


/* compute d_i f^i with finite vol. methods on one node.
   In:
     vlq are vars we reconstruct (we can use cons. vars u here)
     vlu are cons. vars u
     vls are extra source vars we may need (can be NULL)
     rec1d_u_f_lam_midpt = rec. cons u, fluxes f & eigenvals at a midpoint
     numflux             = numerical flux we want
   Out:
     vldivf = div(f(u)) on all inner gridpoints and a piece of div(f(u)) on
              face points */
void fv_divf(tNode *node, tVarList *vldivf, tVarList *vlq,
             tVarList *vlu, tVarList *vls,
             void (*rec1d_u_f_lam_midpt)(tFVinfo *f, tDGinfo *d),
             void (*u_f_lam)(tDGinfo *d),
             void (*numflux)(tDGinfo *d))
{
  tMesh *mesh = vlq->mesh;
  int nqvars = vlq->n;
  int nfvars = vldivf->n;
  int have_XYZ_of_xyz = ( node->pat->XYZ_of_xyz ? 1 : 0 );

  int norms_and_sqrtgdiag_on_midpoints = Getb(coordinates->midpoint_data);
  int idet_dXbdx = coordinates->idet_dXbdx;
  int isqrtgdiagx = coordinates->isqrtgdiagx;
  int iXm_det_dXbdx = coordinates->iXm_det_dXbdx;
  int iYm_det_dXbdx = coordinates->iYm_det_dXbdx;
  int iZm_det_dXbdx = coordinates->iZm_det_dXbdx;
  int im_det_dXbdx[] = { iXm_det_dXbdx, iYm_det_dXbdx, iZm_det_dXbdx };
  int iXm_sqrtgdiagx = coordinates->iXm_sqrtgdiagx;
  int iYm_sqrtgdiagy = coordinates->iYm_sqrtgdiagy;
  int iZm_sqrtgdiagz = coordinates->iZm_sqrtgdiagz;
  int im_sqrtgdiag[] = { iXm_sqrtgdiagx, iYm_sqrtgdiagy, iZm_sqrtgdiagz };

  /* func ptrs for reconstruction */
  double (*rec1d_p)(int n, const double *u, int im, double u_scale);
  double (*rec1d_m)(int n, const double *u, int im, double u_scale);
  double q_scale = 1.; /* typical order of magnitude of fields */
  int nghosts;         /* number of ghost points on each end */
  int add_surface_fluxes; /* whether we set all of divf on faces */
  int use_left_flux;   /* whether we set and use the left fluxes in fnumL */
  int subtract_fi = DGglobals->fv_flux_is_fnum_minus_fi;
  int extrap_mode = DGglobals->fv_divf_extrap_mode;
  double extrap_s1 = DGglobals->fv_divf_extrap_s1;
  double extrap_s2 = DGglobals->fv_divf_extrap_s2;
  int extrap_opt = DGglobals->fv_divf_extrap_opt;

  /* set func ptrs for rec. */
  switch(DGglobals->fv_rec_mode)
  {
  case FV_REC_1:
    /* reconstruct from both sides of midpoint at i0m */
    rec1d_p = rec1d_p_1;
    rec1d_m = rec1d_m_1;
    nghosts = 0;
    break;
  /* use WENO3 from both sides of midpoint at i0m */
  case FV_REC_WENOm3_2:
    rec1d_p = rec1d_p_WENOm3_2;
    rec1d_m = rec1d_m_WENOm3_2;
    nghosts = 0;
    break;
  case FV_REC_WENOm5_2:
    rec1d_p = rec1d_p_WENOm5_2;
    rec1d_m = rec1d_m_WENOm5_2;
    nghosts = 0;
    break;
  case FV_REC_WENOmT_2:
    rec1d_p = rec1d_p_WENOmT_2;
    rec1d_m = rec1d_m_WENOmT_2;
    nghosts = 0;
    break;
  case FV_REC_WENOmZ_2:
    rec1d_p = rec1d_p_WENOmZ_2;
    rec1d_m = rec1d_m_WENOmZ_2;
    nghosts = 0;
    break;
  /* use WENO3 from both sides of midpoint, but copy near boundary */
  case FV_REC_WENOm3_1:
    rec1d_p = rec1d_p_WENOm3_1;
    rec1d_m = rec1d_m_WENOm3_1;
    nghosts = 0;
    break;
  case FV_REC_WENOm5_1:
    rec1d_p = rec1d_p_WENOm5_1;
    rec1d_m = rec1d_m_WENOm5_1;
    nghosts = 0;
    break;
  case FV_REC_WENOmZ_1:
    rec1d_p = rec1d_p_WENOmZ_1;
    rec1d_m = rec1d_m_WENOmZ_1;
    nghosts = 0;
    break;
  /* WENO3 experiments that didn't show good convergence: */
  case FV_REC_WENO3if2away_1:
    rec1d_p = rec1d_p_WENO3_if2away;
    rec1d_m = rec1d_m_WENO3_if2away;
    nghosts = 0;
    break;
  case FV_REC_WENO3if1away_1:
    rec1d_p = rec1d_p_WENO3_if1away;
    rec1d_m = rec1d_m_WENO3_if1away;
    nghosts = 0;
    break;
  case FV_REC_WENO3_2:
    rec1d_p = rec1d_p_WENO3_2;
    rec1d_m = rec1d_m_WENO3_2;
    nghosts = 0;
    break;
  case FV_REC_WENO3_2g:
    rec1d_p = rec1d_p_WENO3_2g;
    rec1d_m = rec1d_m_WENO3_2g;
    nghosts = 1;
    break;
  default:
    errorexit("unknown DGglobals->fv_rec_mode");
  }

  /* do we add in the surface fluxes here already? */
  add_surface_fluxes = DGglobals->fv_divf_adds_surface_fluxes;

  /* do we use left fluxes? */
  use_left_flux = !(DGglobals->fv_divf_use_only_right_flux);

  /* set var list for div of fluxes to zero */
  vlsetconstant_node(node, vldivf, 0.);

  /* RHS */
  {
    tDGinfo *d = alloc_DGinfo(vlu, vls);
    int *n = node->n;
    int maxn = max3(n[0],n[1],n[2]);
    double *Xbm = dmalloc(maxn);
    double *dXb = dmalloc(maxn);
    double *qc[nqvars];          //pointers to data of the q-fields
    double fnumR[nfvars];        //the right fluxes
    double fnumL[nfvars];        //the left fluxes
    double fiR[nfvars];          //the inner fluxes with normal on the right
    double fiL[nfvars];          //the inner fluxes with normal on the left
    int npg = maxn + 2*nghosts;  //number of points in qcg[l]
    double (*qcg)[npg] = dtensor(nqvars*npg);     //array for the q-fields
    double (*di0fi0)[maxn] = dtensor(nfvars*maxn); //array for d_i J*sgd*flux^i
    double *qm_p = dmalloc(nqvars); // array for rec u at one point
    double *qm_m = dmalloc(nqvars);
    int l; /* field index */
    int d_info_midnorm;
    int dir;

    /* set qc to part of qcg without ghosts */
    for(l=0; l<nqvars; l++) qc[l] = &(qcg[l][nghosts]);
    /* NOTE: now qc[l][-1] = qcg[l][0] i.e. ghost on left */

    /* init some vars with 0 */
    forvl(vldivf, l)
      fiR[l] = fiL[l] = fnumR[l] = fnumL[l] = 0.;

    /* write node into d because numflux needs this */
    d->node = node;
    if(norms_and_sqrtgdiag_on_midpoints)
      d_info_midnorm = DGINFO_MIDPTNORM;
    else
      d_info_midnorm = DGINFO_NULL;

    /* get nbsurf and ajsurf already */
    if(nghosts || add_surface_fluxes) get_all_surfaces(node);

    /* add fluxes in each direction to RHS */
    for(dir=0; dir<3; dir++)
    {
      int dir_active = Getb(amr->dir_active[dir]);
      double *sqrtgdiagm = Vard(node, im_sqrtgdiag[dir]);
      double *ooJm = Vard(node, im_det_dXbdx[dir]);
      double *sqrtgdiag = Vard(node, isqrtgdiagx+dir);
      double *ooJ = Vard(node, idet_dXbdx);
      int i,j,k;

      /* do nothing if dir is not active */
      if(!dir_active) continue;

      /* get midpoints */
      set_nm_nodemidpoints_Xb_dir(node, n[dir]-1,0, dir, Xbm);
      set_nm_nodemidpoint_distsXb_dir(node, dir, Xbm, dXb);

      /* loop over plane */
      forplaneN(dir, i,j,k, n, 0)
      {
        int i1 = i1_norm(i,j,k, dir); /* 1st and 2nd index in plane */
        int i2 = i2_norm(i,j,k, dir);
        int i0;                       /* index orthogonal to plane */

        /* fill field arrays qc, i0 runs orth. to plane */
        for(i0=0; i0<n[dir]; i0++)
        {
          int ic,jc,kc, ccc;          /* index of gridpoints */
          ijk_inplaneN(dir, ic,jc,kc, i1,i2, i0);
          ccc = Ind_n(ic,jc,kc, n);

          forvl(vlq, l)
          {
            double *q = Vard(node, Vind(vlq, l));
            qc[l][i0] = q[ccc];
          }
        }

        /* fill in ghost points in qc with value from adjacent node */
        if(nghosts)
        {
          int JK = Ind_n_norm(i,j,k, n, dir);
          forvl(vlq, l)
          {
            int vi = Vind(vlq, l);
            double *qaj;
            /* put adj. val on left face in left ghost: qc[l][-1] */
            qaj = Varaj(node, vi, dir*2);
            if(qaj) qc[l][-1] = qaj[JK];
            else    qc[l][-1] = 1e30; //large value that WENO should ignore
            /* get adj. val on right face in right ghost: qc[l][n] */
            qaj = Varaj(node, vi, dir*2+1);
            if(qaj) qc[l][n[dir]] = qaj[JK];
            else    qc[l][n[dir]] = 1e30;
          }
        }

        /* loop over points in dir */
        for(i0=0; i0<n[dir]; i0++)
        {
          int ic,jc,kc, ccc;          /* index of gridpoints */
          int im0,   cccR;            /* index of right midpoints */
          int im0m1, cccL;            /* index of left midpoints */
          int im,jm,km;
          int i0g0, i0lN;
          double wc;
          double Jgdow_R, Jgdow_L;

          /* set 1d index of left and right midpoint and some flags if we
             are at endpoints */
          if(i0>0) { i0g0=1; im0m1 = i0-1; }
          else     { i0g0=0; im0m1 = i0; /* safe value */ }
          if(i0<n[dir]-1) { i0lN=1; im0 = i0; }
          else            { i0lN=0; im0 = i0-1; /* safe value */ }

          /* gridpoint index and weight */
          ijk_inplaneN(dir, ic,jc,kc, i1,i2, i0);
          ccc = Ind_n(ic,jc,kc, n);
          wc = dXb[i0];

          /* write gridpoint index into d */
          d->i = ic;
          d->j = jc;
          d->k = kc;

          /* set right midpoint index */
          ijk_inplaneN(dir, im,jm,km, i1,i2,im0);
          cccR = Ind_n(im,jm,km, n);

          /* set left midpoint index */
          ijk_inplaneN(dir, im,jm,km, i1,i2,im0m1);
          cccL = Ind_n(im,jm,km, n);

          /* set all factors in flux once */
          if(i0g0 && i0lN) /* in middle */
          {
            Jgdow_R = sqrtgdiagm[cccR] / (ooJm[cccR] * wc);
            Jgdow_L = sqrtgdiagm[cccL] / (ooJm[cccL] * wc);
          }
          else if(i0g0==0 && i0lN) /* left end */
          {
            Jgdow_R = sqrtgdiagm[cccR] / (ooJm[cccR] * wc);
            Jgdow_L = sqrtgdiag[ccc] / (ooJ[ccc] * wc);
          }
          else if(i0g0 && i0lN==0) /* right end */
          {
            Jgdow_R = sqrtgdiag[ccc] / (ooJ[ccc] * wc);
            Jgdow_L = sqrtgdiagm[cccL] / (ooJm[cccL] * wc);
          }
          else /* there is only 1 point */
          {
            Jgdow_R = Jgdow_L = sqrtgdiag[ccc] / (ooJ[ccc] * wc);
          }

          /* Set factors in flux on faces to zero, if we don't add surface
             terms.
             Note: i0g0=0 on left face, i0lN=0 on right face */
          if(!add_surface_fluxes)
          {
            Jgdow_R = i0lN * Jgdow_R;
            Jgdow_L = i0g0 * Jgdow_L;
          }

          /* if i0 has a midpoint to its left */
          if(i0g0)
          {
            if(use_left_flux)
            {
              tFVinfo fv[1];

              /* set fv  */
              fv->nq = nqvars;
              fv->qc = qc;
              fv->npts = n[dir];
              fv->im = im0m1;
              fv->q_scale = q_scale;
              fv->rec1d_p = rec1d_p;
              fv->rec1d_m = rec1d_m;
              fv->qm_p = qm_p;
              fv->qm_m = qm_m;

              d->face = dir*2;
              d->info = d_info_midnorm;

              /* reconstruct q,u and then set fluxes and eigenvalues in d */
              rec1d_u_f_lam_midpt(fv, d);

              /* compute numerical flux directly after rec1d_u_f_lam_midpt,
                 if not set already in rec1d_u_f_lam_midpt */
              if(numflux) numflux(d);

              /* save d->fnum in fnumL for each field and point */
              forvl(vldivf, l) fnumL[l] = d->fnum[l];
              //printDGinfo(d);
            }
            else /* do not separately compute the left flux */
            {
              /* here we set fnumL = -fnumR_{previous point} */
              forvl(vldivf, l) fnumL[l] = -fnumR[l];
              /* for the first pt (i0=0) we do not get here and thus
                 fnumR[l] has already been compouted */
            }
          }
          else /* left end: we need fnumL for sure */
          {
            if(add_surface_fluxes)
            {
              /* compute numerical fluxes on the left side of node */
              d->face = dir*2;        /* normal points to the left */
              d->info = DGINFO_NULL;  /* facepoint is gridpoint */
              u_f_lam(d);
              numflux(d);

              /* save d->fnum in fnumL for each field and point */
              forvl(vldivf, l) fnumL[l] = d->fnum[l];
            }
            /* if add_surface_fluxes=0 we do not set fnumL at all */
            //else
            //{
            //  forvl(vldivf, l) fnumL[l] = 0.;
            //}
          }

          /* if i0 has a midpoint to its right */
          if(i0lN)
          {
            tFVinfo fv[1];

            /* set fv  */
            fv->nq = nqvars;
            fv->qc = qc;
            fv->npts = n[dir];
            fv->im = im0;
            fv->q_scale = q_scale;
            fv->rec1d_p = rec1d_p;
            fv->rec1d_m = rec1d_m;
            fv->qm_p = qm_p;
            fv->qm_m = qm_m;

            d->face = dir*2 + 1;
            d->info = d_info_midnorm;

            /* reconstruct q,u and then set fluxes and eigenvalues in d */
            rec1d_u_f_lam_midpt(fv, d);

            /* compute numerical flux directly after rec1d_u_f_lam_midpt,
               if not set already in rec1d_u_f_lam_midpt */
            if(numflux) numflux(d);

            /* save d->fnum in fnumR for each field */
            forvl(vldivf, l) fnumR[l] = d->fnum[l];
            /* above we have a case for fnumL (with normL=-normR), but I think
               this results in fnumL = -fnumR. So it should be enough to only
               use fnumR. */
            //printDGinfo(d);
            //if(d->fnum[0]>0.00000001 && i0==4) errorexit("STOP");
          }
          else /* right end */
          {
            if(add_surface_fluxes)
            {
              /* compute numerical fluxes on the right side of node */
              d->face = dir*2 + 1;    /* normal points to the right */
              d->info = DGINFO_NULL;  /* facepoint is gridpoint */
              u_f_lam(d);
              numflux(d);

              /* save fluxes of right face in fnumR */
              forvl(vldivf, l) fnumR[l] = d->fnum[l];
            }
            /* if add_surface_fluxes=0 we do not set fnumR at all */
            //else
            //{
            //  forvl(vldivf, l) fnumR[l] = 0.;
            //}
          }


          /* we still need the inner flux computed with normals from both
             the left and right midpoints */
          if(subtract_fi && have_XYZ_of_xyz)
          {
            /* now get d->fi on gridpoint for normal pointing to the right */
            d->info  = d_info_midnorm;
            d->info |= DGINFO_INNONLY; /* it's enough to get ui,fi only */

            d->face = dir*2 + 1; /* ==> normal points to the right */
            u_f_lam(d);
            forvl(vldivf, l) fiR[l] = d->fi[l];

            d->face = dir*2;     /* ==> normal points to the left */
            u_f_lam(d);
            forvl(vldivf, l) fiL[l] = d->fi[l];
          }

          /* get piece of div(flux) in direction dir with FV method */
          forvl(vldivf, l)
          {
            double *df = di0fi0[l];
            double fR  = fnumR[l] - fiR[l]; // fiR can be zero
            double fL  = fnumL[l] - fiL[l]; // fiL can be zero
            df[i0] = Jgdow_R * fR + Jgdow_L * fL;
          }
        } /* end i0 loop */

        /* extrapolate df = d_i0 f^i0 to face */
        if(extrap_mode == FV_DNFN_EXTRAP1)
          forvl(vldivf, l)
          {
            double *df = di0fi0[l];
            rec1d_uface_to_uin_1_Carray(n[dir], df, 0, q_scale,
                                        extrap_s1, extrap_s2, extrap_opt);
          }

        /* final loop over points in dir */
        for(i0=0; i0<n[dir]; i0++)
        {
          int ic,jc,kc, ccc;          /* index of gridpoints */

          /* set points and their index */
          ijk_inplaneN(dir, ic,jc,kc, i1,i2, i0);
          ccc = Ind_n(ic,jc,kc, n);

          /* add to divf */
          forvl(vldivf, l)
          {
            int idivf = Vind(vldivf,l);
            double *divf = Vard_(node, idivf);
            double *df = di0fi0[l];

            /* add (d_i0 f^i0 J) / J term to div(flux) */
            divf[ccc] += df[i0] * ooJ[ccc];
          }
        }
      } /* end plane loop */
    } /* end dir-loop*/

    /* release mem */
    free(qm_m);
    free(qm_p);
    free(di0fi0);
    free(qcg);
    free(dXb);
    free(Xbm);
    free_DGinfo(d);
  }
}



/**************************************************************************/
/* old versions that are not used any longer */
/**************************************************************************/

/* compute d_i f^i with finite vol. methods on one node.
   In:
     vlq are vars we reconstruct (we can use cons. vars u here)
     vlu are cons. vars u
     vls are extra source vars we may need (can be NULL)
     rec1d_u_f_lam_midpt = rec. cons u, fluxes f & eigenvals at a midpoint
     numflux             = numerical flux we want
   Out:
     vldivf = div(f(u)) on all inner gridpoints and a piece of div(f(u)) on
              face points */
void fv_divf__old(tNode *node, tVarList *vldivf, tVarList *vlq,
             tVarList *vlu, tVarList *vls,
             void (*rec1d_u_f_lam_midpt)(tFVinfo *f, tDGinfo *d),
             void (*u_f_lam)(tDGinfo *d),
             void (*numflux)(tDGinfo *d))
{
  tMesh *mesh = vlq->mesh;
  int norms_and_sqrtgdiag_on_midpoints = Getb(coordinates->midpoint_data);
  int nqvars = vlq->n;
  int nfvars = vldivf->n;
  int isqrtgdiagx = coordinates->isqrtgdiagx;
  int idet_dXbdx  = coordinates->idet_dXbdx;
  int iXm_sqrtgdiagx, iYm_sqrtgdiagy, iZm_sqrtgdiagz;
  int iXm_det_dXbdx, iYm_det_dXbdx, iZm_det_dXbdx;
  /* func ptrs for reconstruction */
  double (*rec1d_p)(int n, const double *u, int im, double u_scale);
  double (*rec1d_m)(int n, const double *u, int im, double u_scale);
  double q_scale = 1.; /* typical order of magnitude of fields */
  int nghosts;         /* number of ghost points on each end */
  int add_surface_fluxes; /* whether we set all of divf on faces */
  int use_left_flux;   /* whether we set and use the left fluxes in fnumL */
  int extrap_mode = DGglobals->fv_divf_extrap_mode;
  double extrap_s1 = DGglobals->fv_divf_extrap_s1;
  double extrap_s2 = DGglobals->fv_divf_extrap_s2;
  int extrap_opt = DGglobals->fv_divf_extrap_opt;

  if(norms_and_sqrtgdiag_on_midpoints)
  {
    iXm_sqrtgdiagx = coordinates->iXm_sqrtgdiagx;
    iYm_sqrtgdiagy = coordinates->iYm_sqrtgdiagy;
    iZm_sqrtgdiagz = coordinates->iZm_sqrtgdiagz;
    iXm_det_dXbdx = coordinates->iXm_det_dXbdx;
    iYm_det_dXbdx = coordinates->iYm_det_dXbdx;
    iZm_det_dXbdx = coordinates->iZm_det_dXbdx;
  }
  else
  {
    iXm_sqrtgdiagx = isqrtgdiagx;
    iYm_sqrtgdiagy = isqrtgdiagx + 1;
    iZm_sqrtgdiagz = isqrtgdiagx + 2;
    iXm_det_dXbdx = iYm_det_dXbdx = iZm_det_dXbdx = idet_dXbdx;
  }

  /* set func ptrs for rec. */
  switch(DGglobals->fv_rec_mode)
  {
  case FV_REC_1:
    /* reconstruct from both sides of midpoint at i0m */
    rec1d_p = rec1d_p_1;
    rec1d_m = rec1d_m_1;
    nghosts = 0;
    break;
  /* use WENO3 from both sides of midpoint at i0m */
  case FV_REC_WENOm3_2:
    rec1d_p = rec1d_p_WENOm3_2;
    rec1d_m = rec1d_m_WENOm3_2;
    nghosts = 0;
    break;
  case FV_REC_WENOm5_2:
    rec1d_p = rec1d_p_WENOm5_2;
    rec1d_m = rec1d_m_WENOm5_2;
    nghosts = 0;
    break;
  case FV_REC_WENOmZ_2:
    rec1d_p = rec1d_p_WENOmZ_2;
    rec1d_m = rec1d_m_WENOmZ_2;
    nghosts = 0;
    break;
  /* use WENO3 from both sides of midpoint, but copy near boundary */
  case FV_REC_WENOm3_1:
    rec1d_p = rec1d_p_WENOm3_1;
    rec1d_m = rec1d_m_WENOm3_1;
    nghosts = 0;
    break;
  case FV_REC_WENOm5_1:
    rec1d_p = rec1d_p_WENOm5_1;
    rec1d_m = rec1d_m_WENOm5_1;
    nghosts = 0;
    break;
  case FV_REC_WENOmZ_1:
    rec1d_p = rec1d_p_WENOmZ_1;
    rec1d_m = rec1d_m_WENOmZ_1;
    nghosts = 0;
    break;
  /* WENO3 experiments that didn't show good convergence: */
  case FV_REC_WENO3if2away_1:
    rec1d_p = rec1d_p_WENO3_if2away;
    rec1d_m = rec1d_m_WENO3_if2away;
    nghosts = 0;
    break;
  case FV_REC_WENO3if1away_1:
    rec1d_p = rec1d_p_WENO3_if1away;
    rec1d_m = rec1d_m_WENO3_if1away;
    nghosts = 0;
    break;
  case FV_REC_WENO3_2:
    rec1d_p = rec1d_p_WENO3_2;
    rec1d_m = rec1d_m_WENO3_2;
    nghosts = 0;
    break;
  case FV_REC_WENO3_2g:
    rec1d_p = rec1d_p_WENO3_2g;
    rec1d_m = rec1d_m_WENO3_2g;
    nghosts = 1;
    break;
  default:
    errorexit("unknown DGglobals->fv_rec_mode");
  }

  /* do we add in the surface fluxes here already? */
  add_surface_fluxes = DGglobals->fv_divf_adds_surface_fluxes;

  /* do we use left fluxes? */
  use_left_flux = !(DGglobals->fv_divf_use_only_right_flux);

  /* set var list for div of fluxes to zero */
  vlsetconstant_node(node, vldivf, 0.);

  /* RHS */
  {
    tDGinfo *d = alloc_DGinfo(vlu, vls);
    double *g_sqrtgdiag[3] = { Vard(node, isqrtgdiagx), //sqrtgdiag on gridpts
                               Vard(node, isqrtgdiagx+1),
                               Vard(node, isqrtgdiagx+2) };
    double *m_sqrtgdiag[3] = { Vard(node, iXm_sqrtgdiagx),//sqrtgdiag on midpts
                               Vard(node, iYm_sqrtgdiagy),
                               Vard(node, iZm_sqrtgdiagz) };
    double *g_det_dXbdx = Vard(node, idet_dXbdx); /* 1/J at gridpoint */
    double *m_det_dXbdx[3] = { Vard(node, iXm_det_dXbdx),  /* 1/J at Xmid */
                               Vard(node, iYm_det_dXbdx),  /* 1/J at Ymid */
                               Vard(node, iZm_det_dXbdx) };
    int *n = node->n;
    int maxn = max3(n[0],n[1],n[2]);
    double *Xbm = dmalloc(maxn);
    double *dXb = dmalloc(maxn);
    double *qc[nqvars];          //pointers to data of the q-fields
    double *fnumR[nfvars];       //pointers to data of the right fluxes
    double *fnumL[nfvars];       //pointers to data of the left fluxes
    int npg = maxn + 2*nghosts;  //number of points in qcg[l]
    int npe = maxn + 2;          //number of points in fnumRe[l]
    double (*qcg)[npg] = dtensor(nqvars*npg);     //array for the q-fields
    double (*fnumRe)[npe] = calloc(nfvars, sizeof *fnumRe); //for fluxes
    double (*fnumLe)[npe] = calloc(nfvars, sizeof *fnumLe); //for fluxes
    double (*di0fi0J)[maxn] = dtensor(nfvars*maxn); //array for d_i flux^i*J
    double *qm_p = dmalloc(nqvars); // array for rec u at one point
    double *qm_m = dmalloc(nqvars);
    int l; /* field index */
    int dir;

    /* set qc to part of qcg without ghosts */
    for(l=0; l<nqvars; l++) qc[l] = &(qcg[l][nghosts]);
    /* NOTE: now qc[l][-1] = qcg[l][0] i.e. ghost on left */
    for(l=0; l<nfvars; l++) fnumR[l] = &(fnumRe[l][1]);
    for(l=0; l<nfvars; l++) fnumL[l] = &(fnumLe[l][1]);

    /* write node into d because numflux needs this */
    d->node = node;
    if(norms_and_sqrtgdiag_on_midpoints)
      d->info = DGINFO_MIDPTNORM; // was DGINFO_MIDPT in old version !!!
    else
      d->info = DGINFO_NULL;

    /* get nbsurf and ajsurf already */
    if(nghosts || add_surface_fluxes) get_all_surfaces(node);

    /* add fluxes in each direction to RHS */
    for(dir=0; dir<3; dir++)
    {
      int dir_active = Getb(amr->dir_active[dir]);
      double *sqrtgdiag  = g_sqrtgdiag[dir];
      double *sqrtgdiagm = m_sqrtgdiag[dir];
      double *ooJ  = g_det_dXbdx;
      double *ooJm = m_det_dXbdx[dir];
      int i,j,k;

      /* get midpoints */
      set_nm_nodemidpoints_Xb_dir(node, n[dir]-1,0, dir, Xbm);
      set_nm_nodemidpoint_distsXb_dir(node, dir, Xbm, dXb);

      /* loop over plane */
      forplaneN(dir, i,j,k, n, 0)
      {
        int i1 = i1_norm(i,j,k, dir); /* 1st and 2nd index in plane */
        int i2 = i2_norm(i,j,k, dir);
        int i0;                       /* index orthogonal to plane */
        int ic,jc,kc, ccc;            /* index of gridpoints */

        /* fill field arrays qc, i0 runs orth. to plane */
        for(i0=0; i0<n[dir]; i0++)
        {
          ijk_inplaneN(dir, ic,jc,kc, i1,i2, i0);
          ccc = Ind_n(ic,jc,kc, n);

          forvl(vlq, l)
          {
            double *q = Vard(node, Vind(vlq, l));
            qc[l][i0] = q[ccc];
          }
        }
        /* fill in ghost points in qc with value from adjacent node */
        if(nghosts)
        {
          int JK = Ind_n_norm(i,j,k, n, dir);
          forvl(vlq, l)
          {
            int vi = Vind(vlq, l);
            double *qaj;
            /* put adj. val on left face in left ghost: qc[l][-1] */
            qaj = Varaj(node, vi, dir*2);
            if(qaj) qc[l][-1] = qaj[JK];
            else    qc[l][-1] = 1e30; //large value that WENO should ignore
            /* get adj. val on right face in right ghost: qc[l][n] */
            qaj = Varaj(node, vi, dir*2+1);
            if(qaj) qc[l][n[dir]] = qaj[JK];
            else    qc[l][n[dir]] = 1e30;
          }
        }

        /* loop over points in dir */
        for(i0=0; i0<n[dir]; i0++)
        {
          int im0, im0m1, im,jm,km, cccm1;
          double i0g0, i0lN, Jgd_ow_m, Jgd_ow_m1, wm;

          /* if we have only 1 point do nothing, as there are no midpoints */
//          if(n[dir]<=1) break;

          /* set points and their index */
          ijk_inplaneN(dir, ic,jc,kc, i1,i2, i0);
          ccc = Ind_n(ic,jc,kc, n);

          /* set 1d index of right and left midpoint and some flags if we
             are at endpoints */
          if(i0>0) { i0g0=1; im0m1 = i0-1; }
          else     { i0g0=0; im0m1 = i0; /* safe value */ }
          if(i0<n[dir]-1) { i0lN=1; im0 = i0; }
          else            { i0lN=0; im0 = i0-1; /* safe value */ }

          /* set left midpoint index */
          ijk_inplaneN(dir, im,jm,km, i1,i2,im0m1);
          cccm1 = Ind_n(im,jm,km, n);

          /* set index and face of the left midpoint */
          d->i = im;
          d->j = jm;
          d->k = km;
          d->face = dir*2;

          /* if we want to construct left fluxes AND
             if i0 has a midpoint to its left */
          if(use_left_flux && i0g0)
          {
            tFVinfo fv[1];

            /* set fv  */
            fv->nq = nqvars;
            fv->qc = qc;
            fv->npts = n[dir];
            fv->im = im0m1;
            fv->q_scale = q_scale;
            fv->rec1d_p = rec1d_p;
            fv->rec1d_m = rec1d_m;
            fv->qm_p = qm_p;
            fv->qm_m = qm_m;

            /* reconstruct q,u and then set fluxes and eigenvalues in d */
            rec1d_u_f_lam_midpt(fv, d);

            /* compute numerical flux directly after rec1d_u_f_lam_midpt,
               if not set already in rec1d_u_f_lam_midpt */
            if(numflux) numflux(d);

            /* save d->fnum in fnumL for each field and point */
            forvl(vldivf, l)
              fnumL[l][im0m1] = d->fnum[l];
            //printDGinfo(d);
          }

          /* set index and face of the right midpoint */
          d->i = ic;
          d->j = jc;
          d->k = kc;
          d->face = dir*2 + 1;

          /* if i0 has a midpoint to its right */
          if(i0lN)
          {
            tFVinfo fv[1];

            /* set fv  */
            fv->nq = nqvars;
            fv->qc = qc;
            fv->npts = n[dir];
            fv->im = im0;
            fv->q_scale = q_scale;
            fv->rec1d_p = rec1d_p;
            fv->rec1d_m = rec1d_m;
            fv->qm_p = qm_p;
            fv->qm_m = qm_m;

            /* reconstruct q,u and then set fluxes and eigenvalues in d */
            rec1d_u_f_lam_midpt(fv, d);

            /* compute numerical flux directly after rec1d_u_f_lam_midpt,
               if not set already in rec1d_u_f_lam_midpt */
            if(numflux) numflux(d);

            /* save d->fnum in fnumR for each field and point */
            forvl(vldivf, l)
              fnumR[l][im0] = d->fnum[l];
            /* above we have a case for fnumL (with normL=-normR), but I think
               this results in fnumL = -fnumR. So it should be enough to only
               use fnumR. */
            //printDGinfo(d);
            //if(d->fnum[0]>0.00000001 && i0==4) errorexit("STOP");
          }

          /* factors in flux terms on RHS at right and left midpoint */
          wm  = dXb[i0];
          if(norms_and_sqrtgdiag_on_midpoints && i0g0 && i0lN)
          {
            Jgd_ow_m  = sqrtgdiagm[ccc]  /(ooJm[ccc] * wm);
            Jgd_ow_m1 = sqrtgdiagm[cccm1]/(ooJm[cccm1] * wm);
          }
          else /* get sqrtgdiag on grid points */
          {
            Jgd_ow_m  = sqrtgdiag[ccc]/(ooJ[ccc] * wm);
            Jgd_ow_m1 = sqrtgdiag[ccc]/(ooJ[ccc] * wm);
          }

          /* include flux terms on facepoints */
          if(add_surface_fluxes && dir_active)
          {
            int d_face_sav = d->face; /* save parts of d we may alter */
            int d_info_sav = d->info;

            if(i0 == 0)
            {
              /* set non-zero weight on left boundary */
              Jgd_ow_m1 = sqrtgdiag[ccc]/(ooJ[ccc] * wm);

              /* compute numerical fluxes on the left side of node */
              d->face = dir*2;        /* normal points to the left */
              d->info = DGINFO_NULL;  /* facepoint is grid point */
              u_f_lam(d);
              numflux(d);

              /* save fluxes of left face in fnumR, use fnumR = -fnumL */
              forvl(vldivf, l) fnumR[l][-1] = -( d->fnum[l] );
              if(use_left_flux)
                forvl(vldivf, l) fnumL[l][-1] = d->fnum[l];
            }
            if(i0 == n[dir]-1)
            {
              /* set non-zero weight on right boundary */
              Jgd_ow_m = sqrtgdiag[ccc]/(ooJ[ccc] * wm);

              /* compute numerical fluxes on the right side of node */
              d->face = dir*2 + 1;   /* normal points to the right */
              d->info = DGINFO_NULL; /* facepoint is grid point */
              u_f_lam(d);
              numflux(d);

              /* save fluxes of right face in fnumR */
              forvl(vldivf, l) fnumR[l][i0] = d->fnum[l];
              if(use_left_flux)
                forvl(vldivf, l) fnumL[l][i0] = -( d->fnum[l] );
                                 /* ^-here we used fnumL = -fnumR */
            }

            /* restore altered parts of d */
            d->face = d_face_sav;
            d->info = d_info_sav;
          }
          else
          {
            /* Set factors in flux on faces to zero.
               Note: i0g0=0 on left face, i0lN=0 on right face */
            Jgd_ow_m  = i0lN * Jgd_ow_m;
            Jgd_ow_m1 = i0g0 * Jgd_ow_m1;
          }

          //printf("i0=%d im0=%d im0m1=%d: wm=%g Jgd_ow_m=%g Jgd_ow_m1=%g\n",
          //i0, im0, im0m1, wm, Jgd_ow_m, Jgd_ow_m1);

          /* get piece of div(flux) in direction dir with FV method */
          if(use_left_flux)
            forvl(vldivf, l)
            {
              double *fR  = fnumR[l];
              double *fL  = fnumL[l];
              double *dfJ = di0fi0J[l];
              dfJ[i0] = fR[i0]*Jgd_ow_m + fL[i0-1]*Jgd_ow_m1;
            }
          else /* fv like in BAM, where we only need right flux fnumR */
            forvl(vldivf, l)
            {
              double *fnum = fnumR[l];
              double *dfJ = di0fi0J[l];
              dfJ[i0] = fnum[i0]*Jgd_ow_m - fnum[i0-1]*Jgd_ow_m1;
            }
        } /* end i0 loop */

        /* extrapolate dfJ = d_i0 f^i0 J to face */
        if(extrap_mode == FV_DNFN_EXTRAP1)
          forvl(vldivf, l)
          {
            double *dfJ = di0fi0J[l];
            rec1d_uface_to_uin_1_Carray(n[dir], dfJ, 0, q_scale,
                                        extrap_s1, extrap_s2, extrap_opt);
          }

        /* final loop over points in dir */
        for(i0=0; i0<n[dir]; i0++)
        {
          /* set points and their index */
          ijk_inplaneN(dir, ic,jc,kc, i1,i2, i0);
          ccc = Ind_n(ic,jc,kc, n);

          /* add to divf */
          forvl(vldivf, l)
          {
            int idivf = Vind(vldivf,l);
            double *divf = Vard_(node, idivf);
            double *dfJ = di0fi0J[l];

            /* add (d_i0 f^i0 J) / J term to div(flux) */
            divf[ccc] += dfJ[i0] * ooJ[ccc];
          }
        }
      } /* end plane loop */
    } /* end dir-loop*/

    /* release mem */
    free(qm_m);
    free(qm_p);
    free(di0fi0J);
    free(fnumLe);
    free(fnumRe);
    free(qcg);
    free(dXb);
    free(Xbm);
    free_DGinfo(d);
  }
}


/* print */
void printFVinfo(tFVinfo *fv)
{
  int l, i;
  PRFs(": ");
  printf("nq=%d npts=%d q_scale=%g\n", fv->nq, fv->npts, fv->q_scale);

  for(l=0; l<=fv->nq; l++)
  {
    printf("qc[%d] =", l);
    for(i=0; i<fv->npts; i++) printf(" %g", fv->qc[l][i]);
    printf("\n");
  }

  printf("im=%d\n", fv->im);
  printf("qm_p =");
  for(l=0; l<=fv->nq; l++) printf(" %g", fv->qm_p[l]);
  printf("\n");
  printf("qm_m =");
  for(l=0; l<=fv->nq; l++) printf(" %g", fv->qm_m[l]);
  printf("\n");
}
